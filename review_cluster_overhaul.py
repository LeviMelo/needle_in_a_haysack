#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Local Review Cluster Builder + Structural & Semantic Analytics
- Depth-limited, ranked SR/MA <-> Primary ping-pong expansion until saturation
- Embedding cache (JSONL) to avoid recomputing
- Structural metrics: co-citation & bibliographic coupling
- Semantic: SR vs cited Primary cosine (LM Studio embeddings)
- Outputs: FULL (debug) + COMPACT (audit)

Run: python review_cluster_overhaul.py
"""

import os, sys, time, json, math, hashlib, itertools, logging, statistics
from typing import Dict, List, Set, Tuple, Any, Optional
from xml.etree import ElementTree as ET
import requests

# ---------------------------- Config ---------------------------- #

SEED_PMID = 33591115

# Controlled exploration parameters
MAX_DEPTH_SR = 2            # SR waves
MAX_DEPTH_PRIMARY = 2       # Primary waves (increase from 1 to allow ping-pong)
MAX_PASSES = 6              # safety bound on SR<->Primary outer passes

CAP_SR_TOTAL = 200          # overall cap on SR/MA nodes
CAP_PRIMARY_TOTAL = 800     # overall cap on primary nodes

# Per-node caps (when exceeded, rank & keep top-K)
PER_NODE_CAP_SR_CITERS = 60        # SR ← citers (other SRs)
PER_NODE_CAP_SR_REFS_PRIMARY = 160 # SR → refs (primaries)
PER_NODE_CAP_PRIMARY_CITERS_SR = 80# Primary → citers (SRs)

# Ranking mixture weights (structural only)
W_COREFS_SEED = 2.0    # overlap(candidate.refs, seed.refs)
W_COCIT_SEED  = 1.0    # overlap(candidate.citers, seed.citers)
W_LINK_TO_CUR = 3.0    # any link into current cluster

# iCite
ICITE_BASE = "https://icite.od.nih.gov/api/pubs"
ICITE_TIMEOUT = 40

# E-utilities
EUTILS_BASE  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EUTILS_TOOL  = os.environ.get("EUTILS_TOOL", "review_cluster_overhaul")
EUTILS_EMAIL = os.environ.get("EUTILS_EMAIL", "your_email@example.com")
EUTILS_KEY   = os.environ.get("NCBI_API_KEY")  # optional
EUTILS_TIMEOUT = 40
EUTILS_SLEEP   = 0.35  # polite

# Embeddings (LM Studio)
LMSTUDIO_EMB_URL = "http://127.0.0.1:1234/v1/embeddings"
EMB_MODEL = "text-embedding-qwen3-embedding-0.6b@f16"
EMB_TIMEOUT = 120
EMB_MAX_CHARS = 6000

# Embedding cache
OUT_DIR = "out"
CACHE_DIR = OUT_DIR
CACHE_FILE = os.path.join(CACHE_DIR, f"emb_cache_{EMB_MODEL.replace('/','_').replace('@','_at_')}.jsonl")

# Output
FULL_JSON    = os.path.join(OUT_DIR, "review_cluster.full.json")
COMPACT_JSON = os.path.join(OUT_DIR, "review_cluster.compact.json")
LOG_LEVEL = logging.INFO

# iCite polite batching & sleep
ICITE_PMIDS_PER_CALL = int(os.environ.get("ICITE_PMIDS_PER_CALL", "120"))
ICITE_SLEEP = float(os.environ.get("ICITE_SLEEP", "0.25"))

# On-disk caches (links + pubmed)
ICITE_CACHE_FILE  = os.path.join(OUT_DIR, "icite_links.cache.jsonl")
PUBMED_CACHE_FILE = os.path.join(OUT_DIR, "pubmed_meta.cache.jsonl")

# Publication-type filters
SRMA_PT = {"Meta-Analysis", "Systematic Review"}

# Tight primary whitelist (no “Comparative Study”, etc.)
PRIMARY_WHITELIST = {
    "Randomized Controlled Trial", "Clinical Trial",
    "Clinical Trial, Phase I", "Clinical Trial, Phase II",
    "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial", "Observational Study",
    "Cohort Studies", "Case-Control Studies", "Cross-Sectional Studies",
    "Equivalence Trial"
}


# ---------------------------- Utils ---------------------------- #

def _jsonl_load(path: str, key: str) -> Dict[int, dict]:
    d = {}
    if not os.path.exists(path): return d
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            try:
                rec = json.loads(line)
                if key not in rec: continue
                d[int(rec[key])] = rec
            except Exception:
                continue
    return d

def _jsonl_append(path: str, rec: dict):
    with open(path, "a", encoding="utf-8") as f:
        json.dump(rec, f)
        f.write("\n")

def setup_logging():
    os.makedirs(OUT_DIR, exist_ok=True)
    logging.basicConfig(
        level=LOG_LEVEL,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)]
    )

def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

def chunk_text(s: str, max_chars: int) -> List[str]:
    if not s: return []
    return [s[i:i+max_chars] for i in range(0, len(s), max_chars)]

def l2_norm(v: List[float]) -> float:
    return math.sqrt(sum(x*x for x in v))

def l2_normalize(v: List[float]) -> List[float]:
    n = l2_norm(v)
    return [x/n for x in v] if n > 0 else v

def cosine_unit(a: List[float], b: List[float]) -> float:
    return sum(x*y for x,y in zip(a,b))

def angle_deg_from_cos(c: float) -> float:
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))

def jaccard(a: Set[int], b: Set[int]) -> float:
    if not a and not b: return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0

def summarize(vals: List[float]) -> Dict[str, Any]:
    if not vals: return {"n": 0}
    v = sorted(vals); n = len(v)
    def q(p):
        if n==1: return v[0]
        k=(n-1)*p; f=math.floor(k); c=math.ceil(k)
        return v[f] if f==c else v[f]+(v[c]-v[f])*(k-f)
    out = {"n": n, "mean": sum(v)/n, "median": q(0.5),
           "p10": q(0.10), "p25": q(0.25), "p75": q(0.75), "p90": q(0.90),
           "min": v[0], "max": v[-1]}
    return {k: (round(v, 6) if isinstance(v, float) else v) for k, v in out.items()}

def pearson(xs: List[float], ys: List[float]) -> Optional[float]:
    if len(xs) != len(ys) or len(xs) < 2: return None
    mx = statistics.mean(xs); my = statistics.mean(ys)
    num = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    denx = math.sqrt(sum((x-mx)**2 for x in xs))
    deny = math.sqrt(sum((y-my)**2 for y in ys))
    if denx==0 or deny==0: return None
    return num/(denx*deny)

def rankdata(vals: List[float]) -> List[float]:
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0]*len(vals)
    i = 0
    while i < len(vals):
        j = i
        while j+1 < len(vals) and vals[order[j+1]] == vals[order[i]]:
            j += 1
        r = (i + j + 2) / 2.0
        for k in range(i, j+1):
            ranks[order[k]] = r
        i = j + 1
    return ranks

def spearman(xs: List[float], ys: List[float]) -> Optional[float]:
    if len(xs) != len(ys) or len(xs) < 2: return None
    rx = rankdata(xs); ry = rankdata(ys)
    return pearson(rx, ry)

def md5_text(s: str) -> str:
    return hashlib.md5(s.encode("utf-8")).hexdigest()

# ---------------------------- Classifiers ---------------------------- #

def is_srma(pub_types: List[str]) -> bool:
    s = set(pub_types or [])
    return bool(s & SRMA_PT)

def is_primary(pub_types: List[str]) -> bool:
    s = set(pub_types or [])
    return bool(s & PRIMARY_WHITELIST)

# ---------------------------- iCite ---------------------------- #

def _icite_call(pmids: List[int], want_refs: bool, fields: str) -> List[Dict[str, Any]]:
    """
    Low-level iCite call with small retry for 5xx and dynamic 413 handling in caller.
    """
    params = {
        "pmids": ",".join(str(p) for p in pmids),
        "legacy": "false",
        "format": "json",
        "fl": fields
    }
    if want_refs:
        params["refs"] = "true"
    # light retry on 5xx
    for attempt in range(3):
        r = requests.get(ICITE_BASE, params=params, timeout=ICITE_TIMEOUT)
        if r.status_code >= 500:
            time.sleep(0.5 * (attempt + 1))
            continue
        r.raise_for_status()
        data = r.json()
        if isinstance(data, dict) and "data" in data:
            return data["data"]
        if isinstance(data, list):
            return data
        if isinstance(data, dict) and "pmid" in data:
            return [data]
        return []
    r.raise_for_status()

def icite_fetch_links_light(pmids: List[int]) -> Dict[int, Dict[str, Any]]:
    """
    Fetch only 'citers' cheaply (no refs payload):
      fields = pmid,citedByPmids
    """
    uniq = sorted(set(int(p) for p in pmids if str(p).isdigit()))
    out: Dict[int, Dict[str, Any]] = {}
    if not uniq: return out

    step = max(40, min(ICITE_PMIDS_PER_CALL, 200))
    i = 0
    while i < len(uniq):
        group = uniq[i:i+step]
        try:
            recs = _icite_call(group, want_refs=False, fields="pmid,citedByPmids")
            for d in recs:
                if "pmid" in d:
                    out[int(d["pmid"])] = d
            i += step
            logging.debug(f"iCite light fetched {i}/{len(uniq)} (step={step})")
            time.sleep(ICITE_SLEEP)
        except requests.exceptions.HTTPError as e:
            if getattr(e, "response", None) is not None and e.response.status_code == 413:
                new_step = max(20, step // 2)
                logging.warning(f"iCite 413 (light): shrink batch {step}->{new_step}")
                if new_step == step: raise
                step = new_step
                continue
            raise
    return out

def icite_fetch_links_heavy(pmids: List[int]) -> Dict[int, Dict[str, Any]]:
    """
    Fetch refs (+ citers if present):
      fields = pmid,citedPmids,citedByPmids with refs=true
    """
    uniq = sorted(set(int(p) for p in pmids if str(p).isdigit()))
    out: Dict[int, Dict[str, Any]] = {}
    if not uniq: return out

    step = max(40, min(ICITE_PMIDS_PER_CALL, 200))
    i = 0
    while i < len(uniq):
        group = uniq[i:i+step]
        try:
            recs = _icite_call(group, want_refs=True, fields="pmid,citedPmids,citedByPmids")
            for d in recs:
                if "pmid" in d:
                    out[int(d["pmid"])] = d
            i += step
            logging.debug(f"iCite heavy fetched {i}/{len(uniq)} (step={step})")
            time.sleep(ICITE_SLEEP)
        except requests.exceptions.HTTPError as e:
            if getattr(e, "response", None) is not None and e.response.status_code == 413:
                new_step = max(20, step // 2)
                logging.warning(f"iCite 413 (heavy): shrink batch {step}->{new_step}")
                if new_step == step: raise
                step = new_step
                continue
            raise
    return out

def _ids(x) -> List[int]:
    if not x: return []
    if isinstance(x, list):
        return [int(v) for v in x if str(v).isdigit()]
    return []

_ICITE_MEM = {}  # pmid -> {"pmid":..., "citedPmids":[...], "citedByPmids":[...]}

def icite_cache_load():
    _ICITE_MEM.clear()
    _ICITE_MEM.update(_jsonl_load(ICITE_CACHE_FILE, key="pmid"))

def icite_cache_flush(records: Dict[int, Dict[str, Any]]):
    for pm, rec in records.items():
        # normalize shape
        pmid = int(rec["pmid"])
        out = {"pmid": pmid}
        if "citedPmids" in rec and rec["citedPmids"] is not None:
            out["citedPmids"] = rec["citedPmids"]
        if "citedByPmids" in rec and rec["citedByPmids"] is not None:
            out["citedByPmids"] = rec["citedByPmids"]
        _ICITE_MEM[pmid] = {**_ICITE_MEM.get(pmid, {}), **out}
        _jsonl_append(ICITE_CACHE_FILE, _ICITE_MEM[pmid])

def icite_get_links(pmids: List[int], need_refs: bool) -> Dict[int, Dict[str, Any]]:
    """Return iCite records from cache or fetch missing; heavy-call only if need_refs=True."""
    want = sorted(set(int(p) for p in pmids if str(p).isdigit()))
    missing = [p for p in want if p not in _ICITE_MEM or (need_refs and "citedPmids" not in _ICITE_MEM[p])]
    got: Dict[int, Dict[str, Any]] = {p: _ICITE_MEM[p] for p in want if p in _ICITE_MEM and (not need_refs or "citedPmids" in _ICITE_MEM[p])}

    if missing:
        logging.info(f"iCite fetch ({'heavy' if need_refs else 'light'}): {len(missing)} missing")
        fetched = icite_fetch_links_heavy(missing) if need_refs else icite_fetch_links_light(missing)
        icite_cache_flush(fetched)
        got.update({int(k): v for k, v in fetched.items()})
    return got

def icite_links(doc: Dict[str, Any]) -> Tuple[Set[int], Set[int]]:
    if not doc: return set(), set()
    refs = set(int(x) for x in (doc.get("citedPmids") or []) if str(x).isdigit())
    citers = set(int(x) for x in (doc.get("citedByPmids") or []) if str(x).isdigit())
    return refs, citers

# ---------------------------- PubMed (EFetch) ---------------------------- #

def efetch_pubmed_xml(pmids: List[int]) -> ET.Element:
    ids = ",".join(str(x) for x in pmids)
    params = {
        "db": "pubmed", "id": ids, "retmode": "xml",
        "tool": EUTILS_TOOL, "email": EUTILS_EMAIL,
    }
    if EUTILS_KEY:
        params["api_key"] = EUTILS_KEY
    r = requests.get(f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=EUTILS_TIMEOUT)
    r.raise_for_status()
    return ET.fromstring(r.text)

def xml_text(node: Optional[ET.Element]) -> str:
    return "".join(node.itertext()).strip() if node is not None else ""

def parse_pubmed(root: ET.Element) -> Dict[int, Dict[str, Any]]:
    out: Dict[int, Dict[str, Any]] = {}
    for art in root.findall(".//PubmedArticle"):
        pmid_node = art.find(".//MedlineCitation/PMID")
        if pmid_node is None or not pmid_node.text: continue
        pmid = int(pmid_node.text)
        title = xml_text(art.find(".//MedlineCitation/Article/ArticleTitle"))
        abs_parts = [xml_text(x) for x in art.findall(".//MedlineCitation/Article/Abstract/AbstractText")]
        abstract = "\n".join([p for p in abs_parts if p])
        pts = [xml_text(x) for x in art.findall(".//MedlineCitation/Article/PublicationTypeList/PublicationType")]
        out[pmid] = {"title": title, "abstract": abstract, "pub_types": pts}
    return out

_PUBMED_MEM: Dict[int, Dict[str, Any]] = {}

def pubmed_cache_load():
    _PUBMED_MEM.clear()
    _PUBMED_MEM.update(_jsonl_load(PUBMED_CACHE_FILE, key="pmid"))

def pubmed_cache_flush(records: Dict[int, Dict[str, Any]]):
    for pmid, rec in records.items():
        out = {"pmid": int(pmid), "title": rec.get("title",""), "abstract": rec.get("abstract",""), "pub_types": rec.get("pub_types", [])}
        _PUBMED_MEM[int(pmid)] = out
        _jsonl_append(PUBMED_CACHE_FILE, out)

def get_pubmed_meta(pmids: List[int]) -> Dict[int, Dict[str, Any]]:
    want = sorted(set(int(p) for p in pmids if str(p).isdigit()))
    missing = [p for p in want if p not in _PUBMED_MEM]
    got = {p: _PUBMED_MEM[p] for p in want if p in _PUBMED_MEM}
    if not missing:
        return got

    BATCH = 200
    logging.info(f"EFetch: fetching metadata for {len(missing)} PMIDs in batches of {BATCH}")
    fetched: Dict[int, Dict[str, Any]] = {}
    for i in range(0, len(missing), BATCH):
        sub = missing[i:i+BATCH]
        try:
            root = efetch_pubmed_xml(sub)
            res = parse_pubmed(root)
            fetched.update(res)
        except Exception as e:
            logging.warning(f"EFetch batch failed ({sub[0]}..{sub[-1]}): {e}")
        time.sleep(EUTILS_SLEEP)
    # flush to disk cache
    pubmed_cache_flush(fetched)
    got.update({int(k): v for k, v in fetched.items()})
    return got

# ---------------------------- Embeddings + Cache ---------------------------- #

def lmstudio_embed_single(text: str) -> List[float]:
    payload = {"model": EMB_MODEL, "input": text}
    r = requests.post(LMSTUDIO_EMB_URL, json=payload, timeout=EMB_TIMEOUT)
    if r.status_code != 200:
        raise RuntimeError(f"Embeddings HTTP {r.status_code}: {r.text[:200]}")
    data = r.json()
    return [float(x) for x in data["data"][0]["embedding"]]

def lmstudio_embed_text(text: str) -> List[float]:
    chunks = chunk_text(text, EMB_MAX_CHARS)
    if not chunks: return []
    if len(chunks) == 1:
        return lmstudio_embed_single(chunks[0])
    acc = None
    for ch in chunks:
        v = lmstudio_embed_single(ch)
        if acc is None: acc = [0.0]*len(v)
        for i,x in enumerate(v): acc[i] += x
    return [x/len(chunks) for x in acc]

# Simple JSONL cache: one line per PMID
def load_emb_cache(path: str) -> Dict[int, Dict[str, Any]]:
    cache = {}
    if not os.path.exists(path): return cache
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            try:
                rec = json.loads(line)
                pmid = int(rec["pmid"])
                cache[pmid] = rec
            except Exception:
                continue
    return cache

def append_emb_cache(path: str, pmid: int, text_hash: str, vec: List[float], norm: float, hsh: str):
    with open(path, "a", encoding="utf-8") as f:
        json.dump({"pmid": pmid, "text_hash": text_hash, "vec": vec, "norm": norm, "hash": hsh}, f)
        f.write("\n")

def embed_with_cache(pmid: int, text: str, cache: Dict[int, Dict[str, Any]]) -> Dict[str, Any]:
    th = md5_text(text)
    if pmid in cache and cache[pmid].get("text_hash") == th:
        vec = cache[pmid]["vec"]
        unit = l2_normalize(vec) if vec else []
        return {"dim": len(vec), "norm": cache[pmid].get("norm", l2_norm(vec)), "hash": cache[pmid].get("hash"), "vec": unit}

    vec = lmstudio_embed_text(text) if text else []
    unit = l2_normalize(vec) if vec else []
    h = hashlib.md5((",".join(f"{x:.6f}" for x in unit)).encode("utf-8")).hexdigest() if unit else None
    rec = {"dim": len(vec), "norm": l2_norm(vec) if vec else 0.0, "hash": h, "vec": unit}
    # persist raw (non-normalized) for cache fidelity
    try:
        append_emb_cache(CACHE_FILE, pmid, th, vec, rec["norm"], h or "")
        cache[pmid] = {"pmid": pmid, "text_hash": th, "vec": vec, "norm": rec["norm"], "hash": h or ""}
    except Exception as e:
        logging.warning(f"Cache write failed for {pmid}: {e}")
    return rec

# ---------------------------- Ranking ---------------------------- #

def score_candidate(candidate_pmid: int,
                    cand_refs: Set[int], cand_citers: Set[int],
                    seed_refs: Set[int], seed_citers: Set[int],
                    cluster_nodes: Set[int]) -> float:
    corefs_seed = len(cand_refs & seed_refs)
    cocit_seed  = len(cand_citers & seed_citers)
    link_to_cur = len((cand_refs | cand_citers) & cluster_nodes)
    return (W_COREFS_SEED*corefs_seed) + (W_COCIT_SEED*cocit_seed) + (W_LINK_TO_CUR*link_to_cur)

def pick_top_k(candidates: List[int],
               k: int,
               links: Dict[int, Tuple[Set[int], Set[int]]],
               seed_refs: Set[int], seed_citers: Set[int],
               cluster_nodes: Set[int]) -> List[int]:
    if len(candidates) <= k: return candidates
    scored = []
    for pmid in candidates:
        refs, citers = links.get(pmid, (set(), set()))
        s = score_candidate(pmid, refs, citers, seed_refs, seed_citers, cluster_nodes)
        scored.append((s, pmid))
    scored.sort(reverse=True)
    return [p for _, p in scored[:k]]

# ---------------------------- Cluster Builder (ping-pong) ---------------------------- #

def join_title_abstract(m: Dict[str, Any]) -> str:
    return f"{m.get('title','').strip()}\n\n{m.get('abstract','').strip()}".strip()

def build_cluster(seed_pmid: int) -> Dict[str, Any]:
    # Seed links (heavy only for the seed to get both refs & citers)
    seed_doc = icite_get_links([seed_pmid], need_refs=True).get(seed_pmid, {})
    seed_refs, seed_citers = icite_links(seed_doc)
    logging.info(f"Seed refs: {len(seed_refs)} | Seed citers: {len(seed_citers)}")

    # Initial frontier links: light fetch (citers only) — cheap
    initial_frontier = sorted(list(seed_refs | seed_citers))
    _ = icite_get_links(initial_frontier + [seed_pmid], need_refs=False)

    # Identify SRs among seed citers; primaries among seed refs
    logging.info("Classifying seed neighborhood (SRs among citers; Primaries among refs)…")
    meta_frontier = get_pubmed_meta(initial_frontier) if initial_frontier else {}
    seed_citer_srmas = [p for p in seed_citers if is_srma(meta_frontier.get(p, {}).get("pub_types", []))]
    seed_ref_primaries = [p for p in seed_refs if is_primary(meta_frontier.get(p, {}).get("pub_types", []))]

    SR_set: Set[int] = set()
    PRIM_set: Set[int] = set()
    sr_depth: Dict[int, int] = {}
    pri_depth: Dict[int, int] = {}
    cluster_nodes: Set[int] = {seed_pmid}

    sr_queue: List[Tuple[int,int]] = []
    prim_queue: List[Tuple[int,int]] = []

    for s in seed_citer_srmas:
        if len(SR_set) < CAP_SR_TOTAL and s not in SR_set:
            SR_set.add(s); cluster_nodes.add(s); sr_depth[s] = 0
            sr_queue.append((s, 0))
    for p in seed_ref_primaries:
        if len(PRIM_set) < CAP_PRIMARY_TOTAL and p not in PRIM_set:
            PRIM_set.add(p); cluster_nodes.add(p); pri_depth[p] = 0
            prim_queue.append((p, 0))

    passes = 0
    while passes < MAX_PASSES:
        passes += 1
        logging.info(f"--- Expansion pass {passes} ---")
        added_any = False

        # ---- SR wave ----
        logging.info(f"SR wave: queue={len(sr_queue)} (SR total={len(SR_set)})")
        while sr_queue:
            spmid, d = sr_queue.pop(0)
            if d >= MAX_DEPTH_SR: continue
            # heavy fetch only for THIS SR to get its refs set
            sdoc = icite_get_links([spmid], need_refs=True).get(spmid, {})
            srefs, sciters = icite_links(sdoc)

            cand_prim = list(srefs)
            # classify primaries (batch metadata)
            meta = get_pubmed_meta(cand_prim) if cand_prim else {}
            cand_prim = [x for x in cand_prim if is_primary(meta.get(x, {}).get("pub_types", []))]
            if len(cand_prim) > PER_NODE_CAP_SR_REFS_PRIMARY:
                logging.debug(f"cap prim refs for SR {spmid}: {len(cand_prim)} -> {PER_NODE_CAP_SR_REFS_PRIMARY}")
                # need links (light) to rank
                lnk = icite_get_links(cand_prim, need_refs=False)
                def score(pm):
                    r = set(lnk.get(pm, {}).get("citedPmids", []) or [])
                    c = set(lnk.get(pm, {}).get("citedByPmids", []) or [])
                    return (W_COREFS_SEED*len(r & seed_refs)
                            + W_COCIT_SEED*len(c & seed_citers)
                            + W_LINK_TO_CUR*len((r|c) & cluster_nodes))
                cand_prim = sorted(cand_prim, key=score, reverse=True)[:PER_NODE_CAP_SR_REFS_PRIMARY]
            for p in cand_prim:
                if len(PRIM_set) >= CAP_PRIMARY_TOTAL: break
                if p not in PRIM_set:
                    PRIM_set.add(p); cluster_nodes.add(p); pri_depth[p] = 0
                    prim_queue.append((p, 0))
                    added_any = True

            # SR citers (other SRs) – light links + classify
            cands = list(sciters)
            _ = icite_get_links(cands, need_refs=False)
            meta_sr = get_pubmed_meta(cands) if cands else {}
            cand_sr = [x for x in cands if is_srma(meta_sr.get(x, {}).get("pub_types", []))]
            if len(cand_sr) > PER_NODE_CAP_SR_CITERS:
                logging.debug(f"cap SR citers for SR {spmid}: {len(cand_sr)} -> {PER_NODE_CAP_SR_CITERS}")
                # rank
                lnk = icite_get_links(cand_sr, need_refs=False)
                def score2(pm):
                    r = set(lnk.get(pm, {}).get("citedPmids", []) or [])
                    c = set(lnk.get(pm, {}).get("citedByPmids", []) or [])
                    return (W_COREFS_SEED*len(r & seed_refs)
                            + W_COCIT_SEED*len(c & seed_citers)
                            + W_LINK_TO_CUR*len((r|c) & cluster_nodes))
                cand_sr = sorted(cand_sr, key=score2, reverse=True)[:PER_NODE_CAP_SR_CITERS]
            for s in cand_sr:
                if len(SR_set) >= CAP_SR_TOTAL: break
                if s not in SR_set:
                    SR_set.add(s); cluster_nodes.add(s); sr_depth[s] = d+1
                    if sr_depth[s] < MAX_DEPTH_SR:
                        sr_queue.append((s, d+1))
                    added_any = True

        # ---- Primary wave ----
        logging.info(f"Primary wave: queue={len(prim_queue)} (Primary total={len(PRIM_set)})")
        while prim_queue:
            pmid, d = prim_queue.pop(0)
            if d >= MAX_DEPTH_PRIMARY: continue
            # For primary nodes we only need citers (SRs). Use light.
            pdoc = icite_get_links([pmid], need_refs=False).get(pmid, {})
            prefs, pciters = icite_links(pdoc)
            cands = list(pciters)
            _ = icite_get_links(cands, need_refs=False)
            meta_sr = get_pubmed_meta(cands) if cands else {}
            cand_sr = [x for x in cands if is_srma(meta_sr.get(x, {}).get("pub_types", []))]
            if len(cand_sr) > PER_NODE_CAP_PRIMARY_CITERS_SR:
                logging.debug(f"cap SR citers for Primary {pmid}: {len(cand_sr)} -> {PER_NODE_CAP_PRIMARY_CITERS_SR}")
                lnk = icite_get_links(cand_sr, need_refs=False)
                def score3(pm):
                    r = set(lnk.get(pm, {}).get("citedPmids", []) or [])
                    c = set(lnk.get(pm, {}).get("citedByPmids", []) or [])
                    return (W_COREFS_SEED*len(r & seed_refs)
                            + W_COCIT_SEED*len(c & seed_citers)
                            + W_LINK_TO_CUR*len((r|c) & cluster_nodes))
                cand_sr = sorted(cand_sr, key=score3, reverse=True)[:PER_NODE_CAP_PRIMARY_CITERS_SR]
            for s in cand_sr:
                if len(SR_set) >= CAP_SR_TOTAL: break
                if s not in SR_set:
                    SR_set.add(s); cluster_nodes.add(s); sr_depth[s] = 0
                    sr_queue.append((s, 0))
                    added_any = True
            pri_depth[pmid] = d + 1

        logging.info(f"Pass {passes} result: SR={len(SR_set)} | Primary={len(PRIM_set)} | nodes={len(SR_set)+len(PRIM_set)+1}")
        if not added_any:
            logging.info("No new nodes — cluster saturated or capped.")
            break

    # Now that we’ll compute coupling later, make sure refs are available for all SRs (heavy)
    _ = icite_get_links(list(SR_set), need_refs=True)
    # We’ll fetch primary refs later in analytics if needed (for coupling).
    # Build link view from cache:
    links: Dict[int, Dict[str, Any]] = {}
    for pm in (SR_set | PRIM_set | {seed_pmid}):
        rec = _ICITE_MEM.get(pm, {})
        refs = list(rec.get("citedPmids", []) or [])
        citers = list(rec.get("citedByPmids", []) or [])
        links[pm] = {"refs": sorted(int(x) for x in refs if str(x).isdigit()),
                     "citers": sorted(int(x) for x in citers if str(x).isdigit())}

    return {
        "seed": seed_pmid,
        "seed_refs": seed_refs,
        "seed_citers": seed_citers,
        "SR_set": SR_set,
        "PRIM_set": PRIM_set,
        "links": links,
        "depths": {"sr_depth": sr_depth, "pri_depth": pri_depth}
    }

# ---------------------------- Analytics ---------------------------- #

def compute_structural_and_semantic(
    cluster: Dict[str, Any],
    meta_all: Dict[int, Dict[str, Any]],
) -> Dict[str, Any]:
    SRs = cluster["SR_set"]; PRIM = cluster["PRIM_set"]; seed = cluster["seed"]
    links = cluster["links"]

    def text_for(pmid: int) -> str:
        m = meta_all.get(pmid, {})
        return join_title_abstract(m)

    # Embedding cache warm
    emb_cache = load_emb_cache(CACHE_FILE)

    embeddings: Dict[int, Dict[str, Any]] = {}
    for pmid in (SRs | PRIM | {seed}):
        t = text_for(pmid)
        try:
            rec = embed_with_cache(pmid, t, emb_cache) if t else {"dim":0,"norm":0.0,"hash":None,"vec":[]}
        except Exception as e:
            logging.error(f"Embedding failed for PMID {pmid}: {e}")
            rec = {"dim":0,"norm":0.0,"hash":None,"vec":[]}
        embeddings[pmid] = rec
        
    # Ensure we have refs for primaries as we need coupling; fetch heavy only for missing
    need_refs_for = [p for p in PRIM if not _ICITE_MEM.get(p, {}).get("citedPmids")]
    if need_refs_for:
        logging.info(f"Fetching refs for {len(need_refs_for)} primaries (heavy) for coupling…")
        _ = icite_get_links(need_refs_for, need_refs=True)
        # Refresh local links map for those primaries
        for p in need_refs_for:
            rec = _ICITE_MEM.get(p, {})
            links[p] = {
                "refs": sorted(int(x) for x in (rec.get("citedPmids") or []) if str(x).isdigit()),
                "citers": sorted(int(x) for x in (rec.get("citedByPmids") or []) if str(x).isdigit())
            }

    # SR ↔ PRIMARY pairs where SR cites the primary
    pairs = []
    for s in SRs:
        s_refs = set(links.get(s, {}).get("refs", []))
        for p in (PRIM & s_refs):
            sv = embeddings.get(s, {}).get("vec") or []
            pv = embeddings.get(p, {}).get("vec") or []
            if not sv or not pv:
                continue
            cos = cosine_unit(sv, pv)
            angle = angle_deg_from_cos(cos)

            # Structural metrics
            s_citers = set(links.get(s, {}).get("citers", []))
            p_citers = set(links.get(p, {}).get("citers", []))
            s_refs_full = set(links.get(s, {}).get("refs", []))
            p_refs_full = set(links.get(p, {}).get("refs", []))

            cocit_count = len(s_citers & p_citers)
            cocit_j = jaccard(s_citers, p_citers)
            bc_count = len(s_refs_full & p_refs_full)
            bc_j = jaccard(s_refs_full, p_refs_full)

            pairs.append({
                "sr_pmid": s,
                "primary_pmid": p,
                "cosine": float(cos),
                "angle_deg": float(angle),
                "co_citation_count": cocit_count,
                "co_citation_jaccard": round(cocit_j, 6),
                "biblio_coupling_count": bc_count,
                "biblio_coupling_jaccard": round(bc_j, 6),
            })

    # Coverage per primary
    coverage = {}
    for p in PRIM:
        citing_srs = [s for s in SRs if p in set(cluster["links"].get(s, {}).get("refs", []))]
        coverage[p] = {
            "coverage_count": len(citing_srs),
            "citing_srs": citing_srs  # keep in FULL
        }

    # Gaps: primaries with zero coverage → nearest SRs by cosine
    gaps = []
    for p in PRIM:
        if coverage[p]["coverage_count"] > 0:
            continue
        pv = embeddings.get(p, {}).get("vec") or []
        if not pv:
            gaps.append({"primary_pmid": p, "issue": "no_embedding", "title": meta_all.get(p, {}).get("title","")})
            continue
        sims = []
        for s in SRs:
            sv = embeddings.get(s, {}).get("vec") or []
            if not sv: continue
            sims.append((cosine_unit(sv, pv), s))
        sims.sort(reverse=True)
        top = [{"sr_pmid": s, "cosine": float(c)} for c, s in sims[:5]]
        gaps.append({"primary_pmid": p, "nearest_srs": top, "title": meta_all.get(p, {}).get("title","")})

    # Correlations
    cos_vals = [e["cosine"] for e in pairs]
    co_cit_counts = [e["co_citation_count"] for e in pairs]
    bc_counts = [e["biblio_coupling_count"] for e in pairs]
    corr = {
        "pearson_cos_vs_cocit_count": pearson(cos_vals, co_cit_counts),
        "spearman_cos_vs_cocit_count": spearman(cos_vals, co_cit_counts),
        "pearson_cos_vs_bc_count": pearson(cos_vals, bc_counts),
        "spearman_cos_vs_bc_count": spearman(cos_vals, bc_counts),
    }
    for k,v in list(corr.items()):
        if v is not None:
            corr[k] = round(v, 6)

    return {
        "embeddings": embeddings,  # FULL keeps vec; COMPACT strips them
        "pairs": pairs,
        "coverage": coverage,
        "gaps": gaps,
        "corr": corr,
        "pairwise_cosine_summary": summarize(cos_vals)
    }

# ---------------------------- Compact JSON ---------------------------- #

def make_compact(
    seed: int,
    model: str,
    cluster: Dict[str, Any],
    meta_all: Dict[int, Dict[str, Any]],
    analytics: Dict[str, Any]
) -> Dict[str, Any]:
    SRs = sorted(list(cluster["SR_set"]))
    PRIM = sorted(list(cluster["PRIM_set"]))
    links = cluster["links"]

    def card(pmid: int) -> Dict[str, Any]:
        m = meta_all.get(pmid, {})
        pts = m.get("pub_types", [])
        return {
            "pmid": pmid,
            "title": m.get("title",""),
            "pub_types": pts,
            "refs_count": len(links.get(pmid, {}).get("refs", [])),
            "citers_count": len(links.get(pmid, {}).get("citers", [])),
        }

    sr_cards = [card(p) for p in SRs]
    prim_cards = [card(p) for p in PRIM]

    pairs_view = []
    for e in sorted(analytics["pairs"], key=lambda d: d["cosine"], reverse=True):
        pairs_view.append({
            "sr_pmid": e["sr_pmid"],
            "primary_pmid": e["primary_pmid"],
            "primary_title": meta_all.get(e["primary_pmid"], {}).get("title", "")[:300],
            "cosine": round(e["cosine"], 6),
            "angle_deg": round(e["angle_deg"], 6),
            "co_citation_count": e["co_citation_count"],
            "co_citation_jaccard": e["co_citation_jaccard"],
            "biblio_coupling_count": e["biblio_coupling_count"],
            "biblio_coupling_jaccard": e["biblio_coupling_jaccard"],
        })

    gaps_view = []
    for g in analytics["gaps"][:25]:
        gaps_view.append({
            "primary_pmid": g.get("primary_pmid"),
            "title": g.get("title",""),
            "nearest_srs": [{"sr_pmid": r["sr_pmid"], "cosine": round(r["cosine"], 6)} for r in g.get("nearest_srs", [])]
        })

    cov_counts = [v["coverage_count"] for v in analytics["coverage"].values()]
    coverage_summary = {
        "n_primaries": len(PRIM),
        "n_reviewed": sum(1 for c in cov_counts if c > 0),
        "n_gaps": sum(1 for c in cov_counts if c == 0),
        "coverage_count_summary": summarize(cov_counts)
    }

    return {
        "params": {
            "seed": seed,
            "model": model,
            "max_depth_sr": MAX_DEPTH_SR,
            "max_depth_primary": MAX_DEPTH_PRIMARY,
            "caps": {
                "cap_sr_total": CAP_SR_TOTAL,
                "cap_primary_total": CAP_PRIMARY_TOTAL,
                "per_node_cap_sr_citers": PER_NODE_CAP_SR_CITERS,
                "per_node_cap_sr_refs_primary": PER_NODE_CAP_SR_REFS_PRIMARY,
                "per_node_cap_primary_citers_sr": PER_NODE_CAP_PRIMARY_CITERS_SR
            }
        },
        "cluster_summary": {
            "sr_count": len(SRs),
            "primary_count": len(PRIM),
            "universe_size": len(SRs) + len(PRIM) + 1
        },
        "srmas": sr_cards[:100],
        "primaries": prim_cards[:400],
        "sr_primary_pairs": pairs_view[:500],
        "coverage": coverage_summary,
        "gaps_candidates": gaps_view,
        "cosine_summary": analytics["pairwise_cosine_summary"],
        "correlations": analytics["corr"]
    }

# ---------------------------- Helpers: JSON-safe cluster ---------------------------- #

def sanitize_cluster_for_json(cluster: Dict[str, Any]) -> Dict[str, Any]:
    """Convert sets to sorted lists so json.dump won't crash."""
    return {
        "seed": cluster["seed"],
        "seed_refs": sorted(list(cluster.get("seed_refs", []))),
        "seed_citers": sorted(list(cluster.get("seed_citers", []))),
        "SR_set": sorted(list(cluster.get("SR_set", []))),
        "PRIM_set": sorted(list(cluster.get("PRIM_set", []))),
        "links": cluster.get("links", {}),
        "depths": {
            "sr_depth": cluster.get("depths", {}).get("sr_depth", {}),
            "pri_depth": cluster.get("depths", {}).get("pri_depth", {}),
        }
    }

# ---------------------------- Main ---------------------------- #

def main():
    setup_logging()
    t0 = time.time()
    logging.info(f"Seed PMID: {SEED_PMID}")
    
    icite_cache_load()
    pubmed_cache_load()

    # 1) Build cluster (ranked, depth-controlled ping-pong)
    cluster = build_cluster(SEED_PMID)
    SRs = cluster["SR_set"]; PRIM = cluster["PRIM_set"]
    logging.info(f"Cluster sizes → SR={len(SRs)} | PRIMARY={len(PRIM)} | total={len(SRs)+len(PRIM)+1}")

    # 2) PubMed metadata for all nodes
    all_pmids = sorted(list(SRs | PRIM | {SEED_PMID}))
    meta_all = get_pubmed_meta(all_pmids)
    # Re-filter (safety)
    SRs = {p for p in SRs if is_srma(meta_all.get(p, {}).get("pub_types", []))}
    PRIM = {p for p in PRIM if is_primary(meta_all.get(p, {}).get("pub_types", []))}
    cluster["SR_set"] = SRs
    cluster["PRIM_set"] = PRIM

    # 3) Analytics
    analytics = compute_structural_and_semantic(cluster, meta_all)

    # 4) Persist FULL (debug) with JSON-safe cluster (no sets)
    full = {
        "timestamp_utc": now_iso(),
        "seed": SEED_PMID,
        "params": {
            "max_depth_sr": MAX_DEPTH_SR,
            "max_depth_primary": MAX_DEPTH_PRIMARY,
            "max_passes": MAX_PASSES,
            "caps": {
                "cap_sr_total": CAP_SR_TOTAL,
                "cap_primary_total": CAP_PRIMARY_TOTAL,
                "per_node_cap_sr_citers": PER_NODE_CAP_SR_CITERS,
                "per_node_cap_sr_refs_primary": PER_NODE_CAP_SR_REFS_PRIMARY,
                "per_node_cap_primary_citers_sr": PER_NODE_CAP_PRIMARY_CITERS_SR
            }
        },
        "model": {"endpoint": LMSTUDIO_EMB_URL, "name": EMB_MODEL},
        "cluster": sanitize_cluster_for_json(cluster),
        "pubmed_meta": meta_all,
        "analytics": analytics
    }
    with open(FULL_JSON, "w", encoding="utf-8") as f:
        json.dump(full, f, ensure_ascii=False, indent=2)

    # 5) Persist COMPACT (no vectors)
    compact = make_compact(SEED_PMID, EMB_MODEL, cluster, meta_all, analytics)
    with open(COMPACT_JSON, "w", encoding="utf-8") as f:
        json.dump(compact, f, ensure_ascii=False, indent=2)

    logging.info(f"Saved FULL   → {FULL_JSON}")
    logging.info(f"Saved COMPACT→ {COMPACT_JSON}")

    elapsed = time.time() - t0
    logging.info(f"Elapsed {elapsed:.1f}s")

    # === Topline ===
    logging.info("=== Topline ===")
    cos_sum = analytics.get("pairwise_cosine_summary", {})
    logging.info(f"SRs: {len(SRs)} | Primaries: {len(PRIM)} | SR–Primary pairs: {cos_sum.get('n', 0)}")
    if cos_sum.get("n", 0):
        logging.info(
            "Cosine "
            f"mean={cos_sum.get('mean', 0):.3f} | "
            f"median={cos_sum.get('median', 0):.3f} | "
            f"p25={cos_sum.get('p25', 0):.3f} | "
            f"p75={cos_sum.get('p75', 0):.3f} | "
            f"min={cos_sum.get('min', 0):.3f} | "
            f"max={cos_sum.get('max', 0):.3f}"
        )

    cov_counts = [v["coverage_count"] for v in analytics["coverage"].values()] if analytics.get("coverage") else []
    if cov_counts:
        logging.info(
            f"Coverage: primaries={len(PRIM)} | "
            f"reviewed={sum(1 for c in cov_counts if c > 0)} | "
            f"gaps={sum(1 for c in cov_counts if c == 0)}"
        )

def _main():
    try:
        main()
    except KeyboardInterrupt:
        logging.error("Interrupted by user (Ctrl+C).")
        sys.exit(130)
    except Exception:
        logging.exception("Fatal error:")
        sys.exit(1)

if __name__ == "__main__":
    _main()
