#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Local SR↔Primary Community Builder (PPR → Conductance → Leiden) + Metrics

What it does (end-to-end):
1) Seeded, typed graph expansion using Personalized PageRank (local push) on a bipartite CIT backbone:
   - If node is SR/MA → neighbors = referenced primaries (from refs)
   - If node is Primary → neighbors = SR citers (from citers)
   (We fetch both refs+citers at first touch; NIH iCite is fast; we batch to avoid 413.)
2) Conductance sweep on the discovered nodes (CIT backbone) to pick a compact local community S*.
3) Leiden community detection on S* with semantically-weighted edges for interpretability.
4) Build SR→Primary edges inside S* (both: SR refs→Primary, and Primary citers→SR) for completeness.
5) Compute metrics per SR→Primary edge:
   - cosine (title+abstract embeddings via LM Studio)
   - co-citation (count + Jaccard)
   - bibliographic coupling (count + Jaccard)
   - age difference (years)
6) Coverage & gaps: primaries with zero SR coverage; for each gap, propose top SRs by a composite score.

Outputs:
- out/local_group.compact.json
- out/local_group.edges.csv
- out/local_group.graph.graphml (optional)
"""

import os, sys, re, json, math, time, argparse, hashlib, csv, logging
from collections import defaultdict, deque
from typing import Dict, List, Set, Tuple, Any, Optional

import requests
import numpy as np
import pandas as pd
from tqdm import tqdm
from xml.etree import ElementTree as ET

import igraph as ig
import leidenalg as la

# ---------------------------- Config defaults ---------------------------- #

ICITE_BASE = "https://icite.od.nih.gov/api/pubs"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

OUT_DIR = "out"
CACHE_DIR = "cache"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)

ICITE_BATCH = 150              # keep URL size safe
EFETCH_BATCH = 200
EMB_MAX_CHARS = 6000

# PPR knobs (local push)
ALPHA_DEFAULT = 0.15           # teleport
EPS_DEFAULT = 1e-4             # residual threshold

# Edge weight mixing (for interpretation / Leiden)
THETA1_CIT = 1.0               # backbone if any citation between u,v
THETA2_COCIT = 0.5
THETA3_BC = 0.5
THETA4_SEM = 0.3
THETA5_TIME = 0.1
LAMBDA_TIME = 0.05             # decay per year difference

# Leiden
LEIDEN_RESOLUTION = 0.5        # CPM-like resolution for sensible subtopics

# PubType logic
SRMA_TAGS = {"Systematic Review", "Meta-Analysis"}
PRIMARY_TAGS = {
    "Randomized Controlled Trial",
    "Clinical Trial", "Clinical Trial, Phase I", "Clinical Trial, Phase II",
    "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial",
    "Observational Study", "Cohort Studies", "Case-Control Studies",
    "Cross-Sectional Studies",
    "Equivalence Trial",
    "Comparative Study"  # kept as requested
}

# Gaps ranking weights
BETA_COS = 0.6
BETA_COCIT = 0.2
BETA_BC = 0.2

# ---------------------------- Logging ---------------------------- #

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)]
    )

# ---------------------------- Simple JSONL cache ---------------------------- #

class JsonlCache:
    def __init__(self, path: str, key_field: str):
        self.path = path
        self.key_field = key_field
        self.mem: Dict[str, Any] = {}
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as f:
                for line in f:
                    try:
                        rec = json.loads(line)
                        k = str(rec.get(key_field))
                        if k != "None":
                            self.mem[k] = rec
                    except Exception:
                        continue

    def get(self, key: Any) -> Optional[Dict[str, Any]]:
        return self.mem.get(str(key))

    def set(self, key: Any, value: Dict[str, Any]):
        self.mem[str(key)] = value
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(json.dumps(value, ensure_ascii=False) + "\n")

ICACHE = JsonlCache(os.path.join(CACHE_DIR, "icite.jsonl"), "pmid")
TCACHE = JsonlCache(os.path.join(CACHE_DIR, "efetch.jsonl"), "pmid")
ECACHE = JsonlCache(os.path.join(CACHE_DIR, "embeddings.jsonl"), "key")

# ---------------------------- Helpers ---------------------------- #

def safe_int(x) -> Optional[int]:
    try:
        return int(str(x))
    except Exception:
        return None

def year_from_pubdate(art: ET.Element) -> Optional[int]:
    y = None
    y_node = art.find(".//JournalIssue/PubDate/Year")
    if y_node is not None and y_node.text and y_node.text.strip().isdigit():
        y = int(y_node.text.strip())
    else:
        md = art.find(".//JournalIssue/PubDate/MedlineDate")
        if md is not None and md.text:
            m = re.search(r"(\d{4})", md.text)
            if m:
                y = int(m.group(1))
    return y

def chunk_text(s: str, n: int) -> List[str]:
    if not s:
        return []
    return [s[i:i+n] for i in range(0, len(s), n)]

def l2_norm(v: List[float]) -> float:
    return float(np.sqrt(np.dot(v, v)))

def l2_normalize(v: List[float]) -> List[float]:
    v = np.array(v, dtype=np.float32)
    n = np.linalg.norm(v)
    return (v / n).tolist() if n > 0 else v.tolist()

def cosine(a: List[float], b: List[float]) -> float:
    if not a or not b:
        return 0.0
    return float(np.dot(a, b))

def jaccard(a: Set[int], b: Set[int]) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / float(len(a | b)) if (a or b) else 0.0

def exp_time_decay(y1: Optional[int], y2: Optional[int], lam: float) -> float:
    if not y1 or not y2:
        return 1.0
    return math.exp(-lam * abs(y1 - y2))

def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

# ---------------------------- iCite (refs+citers) ---------------------------- #

def icite_fetch(pmids: List[int]) -> Dict[int, Dict[str, Any]]:
    """Fetch refs+citers for PMIDs (batched), store in ICACHE, return dict."""
    pmids = [int(p) for p in pmids if p is not None]
    missing = [p for p in pmids if ICACHE.get(p) is None]
    out: Dict[int, Dict[str, Any]] = {}
    for i in range(0, len(missing), ICITE_BATCH):
        sub = missing[i:i+ICITE_BATCH]
        params = {
            "pmids": ",".join(str(x) for x in sub),
            "format": "json",
            "refs": "true",
            "legacy": "false"
        }
        r = requests.get(ICITE_BASE, params=params, timeout=90)
        r.raise_for_status()
        data = r.json()
        recs = data.get("data", [])
        for rec in recs:
            pmid = safe_int(rec.get("pmid"))
            if not pmid:
                continue
            # normalize fields
            refs = rec.get("citedPmids")
            citers = rec.get("citedByPmids")
            if refs is None and "references" in rec:
                refs = rec.get("references")
            if citers is None and "cited_by" in rec:
                citers = rec.get("cited_by")
            refs = [int(x) for x in (refs or []) if str(x).isdigit()]
            citers = [int(x) for x in (citers or []) if str(x).isdigit()]
            ICACHE.set(pmid, {"pmid": pmid, "refs": refs, "citers": citers})
    # Assemble output
    for p in pmids:
        rec = ICACHE.get(p)
        if rec:
            out[p] = {"pmid": p, "refs": rec.get("refs", []), "citers": rec.get("citers", [])}
    return out

# ---------------------------- EFetch (title/abstract/types/year) ---------------------------- #

def efetch_meta(pmids: List[int]) -> Dict[int, Dict[str, Any]]:
    pmids = [int(p) for p in pmids if p is not None]
    missing = [p for p in pmids if TCACHE.get(p) is None]
    for i in range(0, len(missing), EFETCH_BATCH):
        sub = missing[i:i+EFETCH_BATCH]
        params = {
            "db": "pubmed",
            "id": ",".join(str(x) for x in sub),
            "retmode": "xml",
        }
        r = requests.get(f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=90)
        r.raise_for_status()
        root = ET.fromstring(r.text)
        for art in root.findall(".//PubmedArticle"):
            pmid_node = art.find(".//PMID")
            if pmid_node is None or not pmid_node.text:
                continue
            pmid = safe_int(pmid_node.text)
            if not pmid:
                continue
            title_node = art.find(".//ArticleTitle")
            title = "".join(title_node.itertext()).strip() if title_node is not None else ""
            abs_parts = [("".join(a.itertext())).strip() for a in art.findall(".//Abstract/AbstractText")]
            abstract = "\n".join([p for p in abs_parts if p])
            pts = [("".join(p.itertext())).strip() for p in art.findall(".//PublicationTypeList/PublicationType")]
            y = year_from_pubdate(art)
            TCACHE.set(pmid, {"pmid": pmid, "title": title, "abstract": abstract, "pub_types": pts, "year": y})
    out: Dict[int, Dict[str, Any]] = {}
    for p in pmids:
        rec = TCACHE.get(p)
        if rec:
            out[p] = {
                "pmid": p,
                "title": rec.get("title", ""),
                "abstract": rec.get("abstract", ""),
                "pub_types": rec.get("pub_types", []),
                "year": rec.get("year", None)
            }
    return out

def is_sr(rec: Dict[str, Any]) -> bool:
    pts = set(rec.get("pub_types", []))
    return bool(pts & SRMA_TAGS)

def is_primary(rec: Dict[str, Any]) -> bool:
    pts = set(rec.get("pub_types", []))
    # If it's SR/MA, do NOT treat as primary even if "Comparative Study" is present.
    if pts & SRMA_TAGS:
        return False
    return bool(pts & PRIMARY_TAGS)


# ---------------------------- Embeddings via LM Studio ---------------------------- #

def lmstudio_embed(text: str, url: str, model: str, timeout: float=90.0) -> List[float]:
    def one(chunk: str) -> List[float]:
        payload = {"model": model, "input": chunk}
        r = requests.post(url, json=payload, timeout=timeout)
        r.raise_for_status()
        data = r.json()
        return [float(x) for x in data["data"][0]["embedding"]]
    chunks = chunk_text(text, EMB_MAX_CHARS)
    if not chunks:
        return []
    if len(chunks) == 1:
        return one(chunks[0])
    acc = None
    for ch in chunks:
        v = one(ch)
        if acc is None:
            acc = np.array(v, dtype=np.float32)
        else:
            acc += np.array(v, dtype=np.float32)
    acc = (acc / len(chunks)).tolist()
    return acc

def get_embedding_for(pmid: int, meta: Dict[int, Dict[str, Any]], lm_url: str, lm_model: str) -> List[float]:
    m = meta.get(pmid, {})
    text = (m.get("title","") + "\n\n" + m.get("abstract","")).strip()
    key = hashlib.md5((str(pmid) + "|" + text).encode("utf-8")).hexdigest()
    cached = ECACHE.get(key)
    if cached and cached.get("vec"):
        return cached["vec"]
    vec = lmstudio_embed(text, lm_url, lm_model)
    unit = l2_normalize(vec) if vec else []
    ECACHE.set(key, {"key": key, "pmid": pmid, "vec": unit})
    return unit

# ---------------------------- Typed neighborhood (CIT backbone) ---------------------------- #

class LocalGraph:
    """
    Maintains:
      nodes: pmid -> node data (refs, citers, type flags, year, emb vec on demand)
      adj: undirected, CIT-only backbone for diffusion (neighbors sets)
    """
    def __init__(self):
        self.nodes: Dict[int, Dict[str, Any]] = {}
        self.adj: Dict[int, Set[int]] = defaultdict(set)

    def ensure_nodes(self, pmids: List[int]):
        # fetch iCite and meta, add to nodes
        if not pmids:
            return
        ic = icite_fetch(pmids)
        md = efetch_meta(pmids)
        for p in pmids:
            if p not in self.nodes:
                self.nodes[p] = {
                    "pmid": p,
                    "refs": set(ic.get(p, {}).get("refs", [])),
                    "citers": set(ic.get(p, {}).get("citers", [])),
                    "title": md.get(p, {}).get("title", ""),
                    "abstract": md.get(p, {}).get("abstract", ""),
                    "pub_types": md.get(p, {}).get("pub_types", []),
                    "year": md.get(p, {}).get("year", None),
                    "is_sr": is_sr(md.get(p, {})),
                    "is_primary": is_primary(md.get(p, {})),
                    "emb": None  # lazy
                }

    def get_neighbors_typed(self, pmid: int) -> List[int]:
        """
        Typed rule:
          - if SR -> neighbors = referenced primaries
          - if Primary -> neighbors = SR citers
        """
        n = self.nodes.get(pmid)
        if not n:
            return []
        out: List[int] = []
        if n["is_sr"]:
            cand = list(n["refs"])
            self.ensure_nodes(cand)
            out = [q for q in cand if self.nodes.get(q, {}).get("is_primary")]
        elif n["is_primary"]:
            cand = list(n["citers"])
            self.ensure_nodes(cand)
            out = [q for q in cand if self.nodes.get(q, {}).get("is_sr")]
        else:
            # non-target types: allow mild connectivity both ways to avoid dead-ends
            # SR citers + refs primaries if they exist
            cand = list(n["refs"] | n["citers"])
            self.ensure_nodes(cand)
            out = [q for q in cand if self.nodes.get(q, {}).get("is_primary") or self.nodes.get(q, {}).get("is_sr")]
        return out

    def ensure_adj(self, pmid: int):
        """
        Always attempt to add typed neighbors.
        Previously we returned early if the node already had neighbors,
        which prevented SR nodes from ever expanding to their referenced primaries.
        """
        self.ensure_nodes([pmid])
        neigh = self.get_neighbors_typed(pmid)
        for q in neigh:
            self.adj[pmid].add(q)
            self.adj[q].add(pmid)  # undirected backbone

    def degree(self, pmid: int) -> int:
        return max(1, len(self.adj.get(pmid, set())))

# ---------------------------- Local PPR (ACL push) on CIT backbone ---------------------------- #

def local_ppr_push(graph: LocalGraph, seeds: List[int], alpha: float, eps: float, max_nodes: int = 5000) -> Dict[int, float]:
    p: Dict[int, float] = defaultdict(float)
    r: Dict[int, float] = defaultdict(float)
    # initialize
    total_seeds = len(seeds)
    for s in seeds:
        graph.ensure_adj(s)  # ensure neighbors for seeds
        r[s] += 1.0 / total_seeds

    # simple queue: nodes with potential push
    queue = deque([s for s in seeds])
    in_q = set(queue)

    def enqueue(v: int):
        if v not in in_q:
            queue.append(v); in_q.add(v)

    discovered: Set[int] = set(seeds)

    while queue and len(discovered) < max_nodes:
        v = queue.popleft(); in_q.discard(v)
        graph.ensure_adj(v)
        deg_v = graph.degree(v)
        if r[v] <= eps * deg_v:
            continue
        # push
        push = alpha * r[v]
        p[v] += push
        rem = (1 - alpha) * r[v]
        r[v] = 0.0
        # neighbors
        neigh = list(graph.adj[v])
        if not neigh:
            continue
        w = 1.0 / len(neigh)  # unweighted CIT backbone
        for u in neigh:
            r[u] += rem * w
            discovered.add(u)
            enqueue(u)
        # also re-enqueue self if residual still big (rare)
        if r[v] > eps * deg_v:
            enqueue(v)

    return dict(p)

# ---------------------------- Conductance on CIT backbone ---------------------------- #

def compute_conductance(graph: LocalGraph, S: Set[int]) -> float:
    if not S:
        return 1.0
    cut = 0
    volS = 0
    # undirected, unit weights for backbone
    for v in S:
        neigh = graph.adj.get(v, set())
        volS += len(neigh)
        for u in neigh:
            if u not in S:
                cut += 1
    # each boundary edge counted once from S
    volOut = sum(len(graph.adj.get(u, set())) for u in set(graph.adj.keys()) - S)
    denom = min(volS, volOut) if min(volS, volOut) > 0 else volS if volS > 0 else 1
    return cut / float(denom) if denom > 0 else 1.0

def sweep_cut(graph: LocalGraph, p_scores: Dict[int, float]) -> Set[int]:
    # order by p/deg
    items = [(v, p_scores.get(v,0.0)/graph.degree(v)) for v in p_scores.keys()]
    items.sort(key=lambda x: x[1], reverse=True)
    best_phi = 1.0
    best_k = 0
    S: Set[int] = set()
    # incremental sweep
    for k, (v, _) in enumerate(items, start=1):
        S.add(v)
        phi = compute_conductance(graph, S)
        if phi < best_phi:
            best_phi = phi
            best_k = k
    return set([items[i][0] for i in range(best_k)])

# ---------------------------- Full weights (for Leiden & metrics) ---------------------------- #

def ensure_embeddings(graph: LocalGraph, pmids: List[int], meta: Dict[int, Dict[str, Any]], lm_url: str, lm_model: str):
    for p in pmids:
        n = graph.nodes.get(p)
        if n is None:
            continue
        if n["emb"] is None:
            vec = get_embedding_for(p, meta, lm_url, lm_model)
            n["emb"] = vec

def blended_weight(u: int, v: int, graph: LocalGraph) -> float:
    nu = graph.nodes[u]; nv = graph.nodes[v]
    # components
    has_cit = (v in nu["refs"]) or (u in nv["refs"])  # either cites either
    w_cit = THETA1_CIT * (1.0 if has_cit else 0.0)
    w_cocit = THETA2_COCIT * jaccard(nu["citers"], nv["citers"])
    w_bc = THETA3_BC * jaccard(nu["refs"], nv["refs"])
    # semantics
    w_sem = THETA4_SEM * cosine(nu["emb"] or [], nv["emb"] or [])
    # time
    w_time = THETA5_TIME * exp_time_decay(nu.get("year"), nv.get("year"), LAMBDA_TIME)
    w = w_cit + w_cocit + w_bc + w_sem + w_time
    return float(max(0.0, w))

# ---------------------------- SR→Primary edges inside S* ---------------------------- #

def sr_primary_edges_in_S(graph: LocalGraph, S: Set[int]) -> List[Tuple[int,int]]:
    S_sr = [p for p in S if graph.nodes[p]["is_sr"]]
    S_pr = [p for p in S if graph.nodes[p]["is_primary"]]
    Sset = set(S)

    edges = set()

    # SR refs → Primary
    for s in S_sr:
        for r in graph.nodes[s]["refs"]:
            if r in Sset and graph.nodes.get(r, {}).get("is_primary"):
                edges.add((s, r))

    # Primary citers → SR
    for p in S_pr:
        for c in graph.nodes[p]["citers"]:
            if c in Sset and graph.nodes.get(c, {}).get("is_sr"):
                edges.add((c, p))

    return sorted(list(edges))

# ---------------------------- Metrics ---------------------------- #

def pearson(xs: List[float], ys: List[float]) -> Optional[float]:
    if len(xs) != len(ys) or len(xs) < 2:
        return None
    x = np.array(xs); y = np.array(ys)
    vx = x - x.mean(); vy = y - y.mean()
    denom = float(np.linalg.norm(vx) * np.linalg.norm(vy))
    return float((vx @ vy) / denom) if denom > 0 else None

def spearman(xs: List[float], ys: List[float]) -> Optional[float]:
    if len(xs) != len(ys) or len(xs) < 2:
        return None
    rx = pd.Series(xs).rank(method="average").to_numpy()
    ry = pd.Series(ys).rank(method="average").to_numpy()
    return pearson(rx.tolist(), ry.tolist())

# ---------------------------- Leiden on S* ---------------------------- #

def leiden_on_S(graph: LocalGraph, S: Set[int], meta: Dict[int, Dict[str, Any]], lm_url: str, lm_model: str) -> Tuple[List[int], List[int], ig.Graph]:
    Slist = sorted(list(S))
    idx_of = {pmid:i for i,pmid in enumerate(Slist)}
    # ensure embeddings for semantic component
    ensure_embeddings(graph, Slist, meta, lm_url, lm_model)

    # Build edges with blended weights; we only add edges where w > 0 and nodes connected via backbone or share structure
    esrc = []
    etgt = []
    ew = []
    # limit to pairs with any structural relation to avoid O(n^2)
    for u in Slist:
        # neighbors we already have (backbone)
        neigh = graph.adj.get(u, set()) & S
        # also include extra candidates: SR→Primary via refs and Primary→SR via citers are already in backbone by construction
        for v in neigh:
            if u < v:
                w = blended_weight(u, v, graph)
                if w > 0:
                    esrc.append(idx_of[u]); etgt.append(idx_of[v]); ew.append(w)

    g = ig.Graph(n=len(Slist), edges=list(zip(esrc, etgt)), directed=False)
    g.vs["pmid"] = Slist
    g.es["weight"] = ew

    if len(g.es) == 0:
        # trivial: one community
        parts = [0]*len(Slist)
        return Slist, parts, g

    # Leiden with weights
    part = la.find_partition(g, la.CPMVertexPartition, weights="weight", resolution_parameter=LEIDEN_RESOLUTION, seed=42)
    # part.membership is a list aligned to vertices
    return Slist, part.membership, g

# ---------------------------- Gaps ranking ---------------------------- #

def rank_gap_sr_for_primary(p: int, SRs: List[int], graph: LocalGraph) -> List[Dict[str, Any]]:
    res = []
    p_node = graph.nodes[p]
    p_citers = p_node["citers"]
    p_refs = p_node["refs"]
    for s in SRs:
        s_node = graph.nodes[s]
        # components
        cosv = cosine(p_node["emb"] or [], s_node["emb"] or [])
        cocit = jaccard(p_citers, s_node["citers"])
        bc = jaccard(p_refs, s_node["refs"])
        score = BETA_COS*cosv + BETA_COCIT*cocit + BETA_BC*bc
        res.append({"sr_pmid": s, "cosine": float(cosv), "cocit_j": float(cocit), "bc_j": float(bc), "score": float(score)})
    res.sort(key=lambda d: d["score"], reverse=True)
    return res[:10]

# ---------------------------- Main ---------------------------- #

def main():
    setup_logging()

    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, required=True, help="Seed PMID")
    ap.add_argument("--alpha", type=float, default=ALPHA_DEFAULT)
    ap.add_argument("--eps", type=float, default=EPS_DEFAULT)
    ap.add_argument("--lm-url", type=str, required=True, help="LM Studio embeddings endpoint URL")
    ap.add_argument("--lm-model", type=str, required=True, help="Embedding model name")
    ap.add_argument("--max-nodes", type=int, default=3000, help="Max nodes to explore in diffusion")
    ap.add_argument("--graphml", action="store_true", help="Write GraphML of S*")
    args = ap.parse_args()

    seed = args.seed
    logging.info(f"Seed PMID: {seed}")

    # Seed load
    lg = LocalGraph()
    lg.ensure_nodes([seed])   # iCite + EFetch
    # Pre-load seed neighbors for backbone
    lg.ensure_adj(seed)

    # PPR local push on backbone
    logging.info("Running local PPR push (CIT backbone)…")
    p_scores = local_ppr_push(lg, [seed], alpha=args.alpha, eps=args.eps, max_nodes=args.max_nodes)
    logging.info(f"Discovered nodes: {len(p_scores)}")

    # Sweep cut to pick S*
    logging.info("Sweep cut (conductance) to pick S*…")
    S_star = sweep_cut(lg, p_scores)
    logging.info(f"Selected community S*: n={len(S_star)}")

    # Metadata map for embeddings
    meta_pmids = sorted(list(S_star))
    meta = efetch_meta(meta_pmids)

    # Ensure embeddings for S* (once)
    logging.info("Embedding S* nodes (cached)…")
    ensure_embeddings(lg, meta_pmids, meta, args.lm_url, args.lm_model)

    # Leiden on S* for interpretability
    logging.info("Leiden clustering on S*…")
    vlist, membership, g_leiden = leiden_on_S(lg, S_star, meta, args.lm_url, args.lm_model)
    # map pmid -> community id
    comm_of = {pmid: cid for pmid, cid in zip(vlist, membership)}
    n_comms = max(membership)+1 if membership else 1
    logging.info(f"Leiden communities: {n_comms}")

    # SR→Primary edges inside S*
    logging.info("Building SR→Primary edges inside S*…")
    edges_sr_p = sr_primary_edges_in_S(lg, S_star)
    logging.info(f"SR→Primary edges: {len(edges_sr_p)}")

    # Metrics per edge
    logging.info("Computing edge metrics (cosine, co-citation, coupling)…")
    rows = []
    cos_vals, cocit_cnts, bc_cnts = [], [], []
    for s, p in tqdm(edges_sr_p, ncols=80):
        ns = lg.nodes[s]; np_ = lg.nodes[p]
        cosv = cosine(ns["emb"] or [], np_["emb"] or [])
        cocit_count = len(ns["citers"] & np_["citers"])
        cocit_j = jaccard(ns["citers"], np_["citers"])
        bc_count = len(ns["refs"] & np_["refs"])
        bc_j = jaccard(ns["refs"], np_["refs"])
        age = None
        if ns.get("year") and np_.get("year"):
            age = abs(ns["year"] - np_["year"])
        rows.append({
            "sr_pmid": s, "primary_pmid": p,
            "sr_title": ns["title"][:300], "primary_title": np_["title"][:300],
            "cosine": round(cosv, 6),
            "cocit_count": cocit_count, "cocit_j": round(cocit_j, 6),
            "bc_count": bc_count, "bc_j": round(bc_j, 6),
            "age_abs_years": age,
            "sr_comm": comm_of.get(s), "primary_comm": comm_of.get(p)
        })
        cos_vals.append(cosv); cocit_cnts.append(cocit_count); bc_cnts.append(bc_count)

    # Coverage & gaps
    S_primaries = [p for p in S_star if lg.nodes[p]["is_primary"]]
    S_srs = [s for s in S_star if lg.nodes[s]["is_sr"]]

    covered = set([p for (_, p) in edges_sr_p])
    gaps = [p for p in S_primaries if p not in covered]

    # rank SRs for each gap
    logging.info(f"Gaps (primaries with no SR in S*): {len(gaps)} — ranking candidates…")
    gap_cards = []
    # ensure SR embeddings present (already ensured)
    for p in gaps:
        recs = rank_gap_sr_for_primary(p, S_srs, lg)
        gp = lg.nodes[p]
        gap_cards.append({
            "primary_pmid": p,
            "title": gp["title"],
            "top_sr_candidates": recs
        })

    # Correlations
    pear_cocit = pearson(cos_vals, cocit_cnts) if cos_vals else None
    spear_cocit = spearman(cos_vals, cocit_cnts) if cos_vals else None
    pear_bc = pearson(cos_vals, bc_cnts) if cos_vals else None
    spear_bc = spearman(cos_vals, bc_cnts) if cos_vals else None

    # CSV edges
    edges_csv = os.path.join(OUT_DIR, "local_group.edges.csv")
    with open(edges_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else
                           ["sr_pmid","primary_pmid","sr_title","primary_title","cosine","cocit_count","cocit_j","bc_count","bc_j","age_abs_years","sr_comm","primary_comm"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # Compact JSON
    compact = {
        "timestamp_utc": now_iso(),
        "seed": seed,
        "params": {
            "alpha": args.alpha,
            "eps": args.eps,
            "max_nodes": args.max_nodes,
            "theta": {"cit": THETA1_CIT, "cocit": THETA2_COCIT, "bc": THETA3_BC, "sem": THETA4_SEM, "time": THETA5_TIME, "lambda_time": LAMBDA_TIME},
            "leiden_resolution": LEIDEN_RESOLUTION,
            "lm_model": args.lm_model
        },
        "community": {
            "size": len(S_star),
            "n_sr": sum(1 for x in S_star if lg.nodes[x]["is_sr"]),
            "n_primary": sum(1 for x in S_star if lg.nodes[x]["is_primary"]),
            "conductance_method": "sweep on CIT backbone",
        },
        "leiden": {
            "n_communities": int(max(membership)+1 if membership else 1),
            "membership_preview": [{"pmid": vlist[i], "comm": int(membership[i])} for i in range(min(len(vlist), 50))]
        },
        "edges_summary": {
            "n_edges": len(rows),
            "cosine_summary": {
                "n": len(cos_vals),
                "mean": float(np.mean(cos_vals)) if cos_vals else None,
                "median": float(np.median(cos_vals)) if cos_vals else None,
                "min": float(np.min(cos_vals)) if cos_vals else None,
                "max": float(np.max(cos_vals)) if cos_vals else None,
                "p25": float(np.percentile(cos_vals, 25)) if cos_vals else None,
                "p75": float(np.percentile(cos_vals, 75)) if cos_vals else None
            },
            "correlations": {
                "pearson_cos_vs_cocit_count": pear_cocit,
                "spearman_cos_vs_cocit_count": spear_cocit,
                "pearson_cos_vs_bc_count": pear_bc,
                "spearman_cos_vs_bc_count": spear_bc
            }
        },
        "coverage": {
            "n_primaries_in_S": len(S_primaries),
            "n_covered": len(covered),
            "n_gaps": len(gaps)
        },
        "gaps": gap_cards[:50],  # keep compact
        "files": {
            "edges_csv": edges_csv
        }
    }
    out_json = os.path.join(OUT_DIR, "local_group.compact.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(compact, f, ensure_ascii=False, indent=2)

    logging.info(f"Wrote: {out_json}")
    logging.info(f"Wrote: {edges_csv}")

    # GraphML (optional)
    if args.graphml:
        graphml_path = os.path.join(OUT_DIR, "local_group.graph.graphml")
        try:
            g_leiden.write_graphml(graphml_path)
            logging.info(f"Wrote: {graphml_path}")
        except Exception as e:
            logging.warning(f"GraphML write failed: {e}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error("Fatal error:", exc_info=True)
        sys.exit(1)
