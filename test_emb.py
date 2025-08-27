#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Seed → SR/MA graphlet + semantic alignment

Pipeline
1) Seed PMID (default: 33591115).
2) iCite: get citers and references of seed.
3) PubMed EFetch: fetch title, abstract, publication types for citers; keep SR/MA (pub_types include 'Meta-Analysis' or 'Systematic Review').
4) For each SR/MA:
   - iCite: get its citers+references.
   - Compute COMMON neighbors with seed: (refs_seed ∪ citers_seed) ∩ (refs_SR ∪ citers_SR).
5) Universe U = {seed} ∪ SRMAs ∪ COMMON neighbors ∪ all primary evidence that each SR/MA cites.
6) EFetch: fetch title/abstract/pub_types for U.
7) Embed all (title+abstract) via LM Studio embeddings endpoint (model: text-embedding-qwen3-embedding-0.6b@f16).
8) For each SR/MA, find its referenced PRIMARY evidence within U, compute cosine(SR, primary).
9) Save full JSON (data + metrics) and print verbose report.

Notes
- iCite fields used: 'citedByPmids' (citers), 'citedPmids' (references) with legacy=false API. (Also fall back to legacy names if needed.)
- E-utilities etiquette: include 'tool' and 'email' parameters; throttle to ≤3 req/s without API key.
- Publication types: use PubMed EFetch XML PublicationTypeList/PublicationType.
"""

import os
import sys
import time
import json
import math
import queue
import hashlib
import logging
import textwrap
import itertools
from typing import Dict, List, Set, Tuple, Any, Optional
from xml.etree import ElementTree as ET

import requests

# ---------------------------- Config ---------------------------- #

SEED_PMID = 33591115

ICITE_BASE = "https://icite.od.nih.gov/api"
ICITE_TIMEOUT = 30

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EUTILS_TOOL = os.environ.get("EUTILS_TOOL", "sr_semantic_eval")
EUTILS_EMAIL = os.environ.get("EUTILS_EMAIL", "your_email@example.com")
EUTILS_TIMEOUT = 30
EUTILS_SLEEP = 0.35  # polite throttling (~ <=3 req/s)

LMSTUDIO_EMB_URL = "http://127.0.0.1:1234/v1/embeddings"
EMB_MODEL = "text-embedding-qwen3-embedding-0.6b@f16"
EMB_TIMEOUT = 120
EMB_MAX_CHARS = 6000  # crude chunking safeguard for very long abstracts

OUT_DIR = "out"
OUT_JSON = os.path.join(OUT_DIR, "sr_semantic_eval.json")
LOG_LEVEL = logging.INFO

# Primary evidence publication types (exact PubMed pt strings)
PRIMARY_PT = {
    "Randomized Controlled Trial",
    "Controlled Clinical Trial",
    "Clinical Trial",
    "Clinical Trial, Phase I",
    "Clinical Trial, Phase II",
    "Clinical Trial, Phase III",
    "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial",
    "Observational Study",
    "Case Reports",
    "Comparative Study",
    "Evaluation Study",
    "Multicenter Study",
    "Validation Study",
    "Prospective Studies",
    "Retrospective Studies",
    "Cohort Studies",
    "Cross-Sectional Studies",
    "Case-Control Studies",
    "Twin Study",
    "Meta-Analysis as Topic",  # not primary, but leave for completeness in parsing
}

# SR/MA types
SRMA_PT = {"Meta-Analysis", "Systematic Review"}

# ---------------------------- Utils ---------------------------- #

def setup_logging():
    os.makedirs(OUT_DIR, exist_ok=True)
    logging.basicConfig(
        level=LOG_LEVEL,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)]
    )

def chunk_text(s: str, max_chars: int) -> List[str]:
    if not s:
        return []
    return [s[i:i+max_chars] for i in range(0, len(s), max_chars)]

def l2_norm(vec: List[float]) -> float:
    return math.sqrt(sum(x * x for x in vec))

def l2_normalize(vec: List[float]) -> List[float]:
    n = l2_norm(vec)
    return [x / n for x in vec] if n > 0 else vec

def cosine_unit(a: List[float], b: List[float]) -> float:
    # expects both pre-normalized
    return sum(x * y for x, y in zip(a, b))

def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

# ---------------------------- iCite ---------------------------- #

def icite_get_pub(pmid: int, legacy: bool = False) -> Dict[str, Any]:
    """
    GET /api/pubs/{pmid}; when legacy=false, field names include:
      - citedByPmids (citers)
      - citedPmids   (references)
    """
    url = f"{ICITE_BASE}/pubs/{pmid}"
    params = {}
    if legacy is False:
        # iCite doc indicates using legacy=false via the /api/pubs?pmids=... form.
        # /pubs/{pmid} returns modern keys by default today, but we handle both.
        pass
    r = requests.get(url, timeout=ICITE_TIMEOUT)
    r.raise_for_status()
    return r.json()

def icite_get_many(pmids: List[int], legacy: bool = False) -> Dict[int, Dict[str, Any]]:
    """
    GET /api/pubs?pmids=a,b,c (up to 1000). If legacy=false, modern keys:
      citedByPmids, citedPmids
    """
    if not pmids:
        return {}
    batch = ",".join(str(p) for p in pmids)
    url = f"{ICITE_BASE}/pubs"
    params = {"pmids": batch}
    if legacy is False:
        params["legacy"] = "false"
    r = requests.get(url, params=params, timeout=ICITE_TIMEOUT)
    r.raise_for_status()
    data = r.json()
    out = {}
    # Response may be list or dict with 'data'; normalize
    if isinstance(data, list):
        for d in data:
            if "pmid" in d:
                out[int(d["pmid"])] = d
    elif isinstance(data, dict) and "data" in data and isinstance(data["data"], list):
        for d in data["data"]:
            if "pmid" in d:
                out[int(d["pmid"])] = d
    else:
        # /pubs?pmids=... may also return a single pub JSON
        if "pmid" in data:
            out[int(data["pmid"])] = data
    return out

def extract_refs_citers(icite_doc: Dict[str, Any]) -> Tuple[Set[int], Set[int]]:
    """
    Return (references, citers) sets, handling legacy and modern keys.
    Modern (legacy=false mapping): citedPmids (refs), citedByPmids (citers).
    Legacy: 'references' (refs), 'cited_by' (citers).
    """
    refs = set()
    citers = set()
    if "citedPmids" in icite_doc or "citedByPmids" in icite_doc:
        refs = set(int(x) for x in icite_doc.get("citedPmids", []) if isinstance(x, int) or (isinstance(x, str) and x.isdigit()))
        citers = set(int(x) for x in icite_doc.get("citedByPmids", []) if isinstance(x, int) or (isinstance(x, str) and x.isdigit()))
    else:
        # legacy names
        refs = set(int(x) for x in icite_doc.get("references", []) if isinstance(x, int) or (isinstance(x, str) and x.isdigit()))
        citers = set(int(x) for x in icite_doc.get("cited_by", []) if isinstance(x, int) or (isinstance(x, str) and x.isdigit()))
    return refs, citers

# ---------------------------- E-utilities (PubMed) ---------------------------- #

def efetch_pubmed_xml(pmids: List[int]) -> ET.Element:
    """
    Fetch PubMed XML (EFetch) for a batch of PMIDs.
    Parse-able tags include:
      - MedlineCitation/Article/ArticleTitle
      - MedlineCitation/Article/Abstract/AbstractText
      - MedlineCitation/Article/PublicationTypeList/PublicationType
    """
    ids = ",".join(str(x) for x in pmids)
    params = {
        "db": "pubmed",
        "id": ids,
        "retmode": "xml",
        "tool": EUTILS_TOOL,
        "email": EUTILS_EMAIL,
    }
    r = requests.get(f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=EUTILS_TIMEOUT)
    r.raise_for_status()
    return ET.fromstring(r.text)

def xml_text(node: Optional[ET.Element]) -> str:
    return "".join(node.itertext()).strip() if node is not None else ""

def parse_pubmed_records(root: ET.Element) -> Dict[int, Dict[str, Any]]:
    """
    Returns dict[pmid] = {
      'title': str, 'abstract': str,
      'pub_types': [str, ...]
    }
    """
    ns = {}  # PubMed XML has no special namespaces for these tags
    out: Dict[int, Dict[str, Any]] = {}
    for art in root.findall(".//PubmedArticle", ns):
        pmid_node = art.find(".//MedlineCitation/PMID", ns)
        if pmid_node is None or not pmid_node.text:
            continue
        pmid = int(pmid_node.text)
        title = xml_text(art.find(".//MedlineCitation/Article/ArticleTitle", ns))
        abs_parts = [xml_text(x) for x in art.findall(".//MedlineCitation/Article/Abstract/AbstractText", ns)]
        abstract = "\n".join([p for p in abs_parts if p])
        pts = [xml_text(x) for x in art.findall(".//MedlineCitation/Article/PublicationTypeList/PublicationType", ns)]
        out[pmid] = {
            "title": title,
            "abstract": abstract,
            "pub_types": pts
        }
    return out

def get_pubmed_metadata(pmids: List[int]) -> Dict[int, Dict[str, Any]]:
    """
    Batched EFetch with politeness delay. Returns minimal dict for all pmids.
    """
    res: Dict[int, Dict[str, Any]] = {}
    BATCH = 150
    for i in range(0, len(pmids), BATCH):
        batch = pmids[i:i+BATCH]
        root = efetch_pubmed_xml(batch)
        res.update(parse_pubmed_records(root))
        time.sleep(EUTILS_SLEEP)
    return res

# ---------------------------- Embedding (LM Studio) ---------------------------- #

def lmstudio_embed(text: str) -> List[float]:
    """
    Call LM Studio OpenAI-compatible embeddings endpoint.
    For very long text, chunk and mean-pool.
    """
    def embed_once(chunk: str) -> List[float]:
        payload = {"model": EMB_MODEL, "input": chunk}
        r = requests.post(LMSTUDIO_EMB_URL, json=payload, timeout=EMB_TIMEOUT)
        if r.status_code != 200:
            raise RuntimeError(f"Embeddings HTTP {r.status_code}: {r.text[:200]}")
        data = r.json()
        emb = data["data"][0]["embedding"]
        return [float(x) for x in emb]

    chunks = chunk_text(text, EMB_MAX_CHARS)
    if not chunks:
        return []
    if len(chunks) == 1:
        return embed_once(chunks[0])

    # mean-pool chunk embeddings
    acc = None
    for ch in chunks:
        vec = embed_once(ch)
        if acc is None:
            acc = [0.0] * len(vec)
        for i, v in enumerate(vec):
            acc[i] += v
    return [v / len(chunks) for v in acc]

# ---------------------------- Core Logic ---------------------------- #

def is_srma(pub_types: List[str]) -> bool:
    s = set(pub_types or [])
    return bool(s & SRMA_PT)

def is_primary(pub_types: List[str]) -> bool:
    s = set(pub_types or [])
    if s & PRIMARY_PT:
        return True
    # Allow a simple heuristic on "Trial" to avoid missing specific subtypes
    return any("Trial" in pt for pt in s)

def join_title_abstract(meta: Dict[str, Any]) -> str:
    return f"{meta.get('title','').strip()}\n\n{meta.get('abstract','').strip()}".strip()

def summarize(vals: List[float]) -> Dict[str, Any]:
    if not vals:
        return {"n": 0}
    v = sorted(vals)
    n = len(v)
    def q(p):
        if n == 1:
            return v[0]
        k = (n-1)*p
        f = math.floor(k)
        c = math.ceil(k)
        if f == c: return v[f]
        return v[f] + (v[c]-v[f])*(k-f)
    return {
        "n": n,
        "mean": sum(v)/n,
        "median": q(0.5),
        "p10": q(0.10),
        "p25": q(0.25),
        "p75": q(0.75),
        "p90": q(0.90),
        "min": v[0],
        "max": v[-1]
    }

def main():
    setup_logging()
    t0 = time.time()
    logging.info(f"Seed PMID: {SEED_PMID}")

    # Seed iCite
    seed_doc = icite_get_pub(SEED_PMID)
    seed_refs, seed_citers = extract_refs_citers(seed_doc)
    logging.info(f"Seed refs: {len(seed_refs)} | Seed citers: {len(seed_citers)}")

    # Candidate SR/MAs = citers of seed filtered by pub types
    # First fetch PubMed metadata for all citers
    citer_list = sorted(seed_citers)
    meta_citers = get_pubmed_metadata(citer_list)
    srma_pmids = [p for p, m in meta_citers.items() if is_srma(m.get("pub_types", []))]
    logging.info(f"SR/MAs citing seed: {len(srma_pmids)}")

    # For each SR/MA, fetch iCite neighborhood
    sr_icite = icite_get_many(srma_pmids, legacy=False)
    sr_neighbors = {}  # pmid -> dict(refs, citers, common_with_seed)
    seed_neighbors = seed_refs | seed_citers
    for spmid, doc in sr_icite.items():
        srefs, sciters = extract_refs_citers(doc)
        common = (seed_neighbors) & (srefs | sciters)
        sr_neighbors[spmid] = {
            "refs": sorted(srefs),
            "citers": sorted(sciters),
            "common_with_seed": sorted(common),
        }

    # Universe U
    universe: Set[int] = {SEED_PMID} | set(srma_pmids)
    for s in sr_neighbors.values():
        universe.update(s["common_with_seed"])
    # Also ensure: include all primary evidence referenced BY each SR/MA
    all_sr_refs = set(itertools.chain.from_iterable(s["refs"] for s in sr_neighbors.values()))
    universe.update(all_sr_refs)

    logging.info(f"Universe size (seed + SR/MAs + commons + SR refs): n={len(universe)}")

    # Fetch PubMed metadata for universe
    meta_all = get_pubmed_metadata(sorted(universe))
    # Add any missing SR/MA metadata (if not already fetched via universe batch)
    for spmid in srma_pmids:
        if spmid not in meta_all and spmid in meta_citers:
            meta_all[spmid] = meta_citers[spmid]

    # Build role flags
    roles = {}
    for pmid, m in meta_all.items():
        pts = m.get("pub_types", [])
        roles[pmid] = {
            "is_seed": (pmid == SEED_PMID),
            "is_srma": is_srma(pts),
            "is_primary": is_primary(pts),
        }

    # Embed all (title+abstract)
    embeddings: Dict[int, Dict[str, Any]] = {}
    for pmid, m in meta_all.items():
        text = join_title_abstract(m)
        try:
            vec = lmstudio_embed(text) if text else []
        except Exception as e:
            logging.error(f"Embedding failed for PMID {pmid}: {e}")
            vec = []
        dim = len(vec)
        norm = l2_norm(vec) if dim > 0 else 0.0
        unit = l2_normalize(vec) if dim > 0 else []
        embeddings[pmid] = {
            "dim": dim,
            "norm": norm,
            "hash": hashlib.md5((",".join(f"{x:.6f}" for x in unit)).encode("utf-8")).hexdigest() if unit else None,
            "vec": unit,  # normalized vector so cosine == dot
        }

    # For each SR/MA, compute cosine with its referenced PRIMARY evidence
    pairwise = []  # list of dicts {sr_pmid, primary_pmid, cosine, in_common_with_seed: bool}
    for spmid in srma_pmids:
        s_meta = embeddings.get(spmid, {})
        s_vec = s_meta.get("vec") or []
        if not s_vec:
            continue
        s_ref_list = set(sr_neighbors.get(spmid, {}).get("refs", []))
        # keep only primary ones we actually fetched
        prims = [p for p in s_ref_list if roles.get(p, {}).get("is_primary", False) and p in embeddings]
        for p in prims:
            p_vec = embeddings[p]["vec"]
            if not p_vec:
                continue
            cos = cosine_unit(s_vec, p_vec)
            pairwise.append({
                "sr_pmid": spmid,
                "primary_pmid": p,
                "cosine": float(cos),
                "in_common_with_seed": bool(p in set(sr_neighbors.get(spmid, {}).get("common_with_seed", [])))
            })

    # Metrics
    cos_vals = [x["cosine"] for x in pairwise]
    overall = summarize(cos_vals)
    by_sr = {}
    for spmid, group in itertools.groupby(sorted(pairwise, key=lambda d: d["sr_pmid"]), key=lambda d: d["sr_pmid"]):
        vals = [g["cosine"] for g in group]
        by_sr[spmid] = summarize(vals)

    # Build JSON
    data = {
        "seed_pmid": SEED_PMID,
        "timestamp_utc": now_iso(),
        "icite_note": "iCite /api/pubs with legacy=false exposes 'citedByPmids' and 'citedPmids' (citers/references).",
        "lmstudio": {
            "endpoint": LMSTUDIO_EMB_URL,
            "model": EMB_MODEL
        },
        "srma_pmids": sorted(srma_pmids),
        "seed_neighbors": {
            "refs": sorted(seed_refs),
            "citers": sorted(seed_citers),
        },
        "sr_neighbors": sr_neighbors,   # per SR: refs, citers, common_with_seed
        "universe_pmids": sorted(list(universe)),
        "pubmed_meta": meta_all,        # title, abstract, pub_types
        "roles": roles,                 # per pmid flags
        "embeddings": embeddings,       # per pmid: dim, norm, hash, vec (normalized)
        "sr_primary_pairs": pairwise,   # list of {sr_pmid, primary_pmid, cosine, in_common_with_seed}
        "metrics": {
            "overall_pairwise_cosine": overall,
            "by_sr": by_sr
        }
    }

    # Save
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)

    t1 = time.time()
    logging.info(f"Saved JSON → {OUT_JSON} | elapsed {t1-t0:.1f}s")

    # --------- Verbose console summary --------- #
    def label(ptypes: List[str]) -> str:
        s = set(ptypes or [])
        if is_srma(ptypes): return "SR/MA"
        if is_primary(ptypes): return "PRIMARY"
        return "OTHER"

    # Topline
    logging.info("=== Topline ===")
    logging.info(f"SR/MAs citing seed: {len(srma_pmids)} | Universe size: {len(universe)} | SR–primary pairs: {overall.get('n',0)}")
    logging.info("Cosine (SR vs its PRIMARY refs) → "
                 f"mean={overall.get('mean','NA'):.3f} | median={overall.get('median','NA'):.3f} "
                 f"| p25={overall.get('p25','NA'):.3f} | p75={overall.get('p75','NA'):.3f} "
                 f"| min={overall.get('min','NA'):.3f} | max={overall.get('max','NA'):.3f}" if overall.get('n',0) else "No pairs.")

    # Show a few SRs
    for spmid in srma_pmids[:10]:
        m = meta_all.get(spmid, {})
        t = m.get("title", "")[:140].replace("\n"," ")
        stats = by_sr.get(spmid, {})
        logging.info(f"[SR/MA {spmid}] {t} ... | pairs={stats.get('n',0)} | mean={stats.get('mean','NA'):.3f} | med={stats.get('median','NA'):.3f}")

    # A few example pairs (top-5 by cosine)
    top5 = sorted(pairwise, key=lambda d: d["cosine"], reverse=True)[:5]
    if top5:
        logging.info("=== Top 5 SR–PRIMARY pairs by cosine ===")
        for e in top5:
            s = e["sr_pmid"]; p = e["primary_pmid"]
            st = meta_all.get(s, {}).get("title","")[:110].replace("\n"," ")
            pt = meta_all.get(p, {}).get("title","")[:110].replace("\n"," ")
            logging.info(f"cos={e['cosine']:.3f} | SR={s}: {st} ... || Primary={p}: {pt} ... | in_common_with_seed={e['in_common_with_seed']}")

if __name__ == "__main__":
    main()