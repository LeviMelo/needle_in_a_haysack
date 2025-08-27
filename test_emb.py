#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, time, json, math, hashlib, itertools, logging
from typing import Dict, List, Set, Tuple, Any, Optional
from xml.etree import ElementTree as ET
import requests

# ---------------------------- Config ---------------------------- #

SEED_PMID = 33591115

ICITE_BASE = "https://icite.od.nih.gov/api/pubs"
ICITE_TIMEOUT = 40

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EUTILS_TOOL  = os.environ.get("EUTILS_TOOL", "sr_semantic_eval")
EUTILS_EMAIL = os.environ.get("EUTILS_EMAIL", "your_email@example.com")
EUTILS_KEY   = os.environ.get("NCBI_API_KEY")  # optional
EUTILS_TIMEOUT = 40
EUTILS_SLEEP   = 0.35  # polite

LMSTUDIO_EMB_URL = "http://127.0.0.1:1234/v1/embeddings"
EMB_MODEL = "text-embedding-qwen3-embedding-0.6b@f16"
EMB_TIMEOUT = 120
EMB_MAX_CHARS = 6000

OUT_DIR = "out"
OUT_JSON = os.path.join(OUT_DIR, "sr_semantic_eval.json")
LOG_LEVEL = logging.INFO

PRIMARY_PT = {
    "Randomized Controlled Trial", "Controlled Clinical Trial",
    "Clinical Trial", "Clinical Trial, Phase I", "Clinical Trial, Phase II",
    "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial", "Observational Study", "Case Reports",
    "Comparative Study", "Evaluation Study", "Multicenter Study",
    "Validation Study", "Prospective Studies", "Retrospective Studies",
    "Cohort Studies", "Cross-Sectional Studies", "Case-Control Studies", "Twin Study",
}
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
    if not s: return []
    return [s[i:i+max_chars] for i in range(0, len(s), max_chars)]

def l2_norm(v: List[float]) -> float:
    return math.sqrt(sum(x*x for x in v))

def l2_normalize(v: List[float]) -> List[float]:
    n = l2_norm(v)
    return [x/n for x in v] if n > 0 else v

def cosine_unit(a: List[float], b: List[float]) -> float:
    return sum(x*y for x,y in zip(a,b))

def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

# ---------------------------- iCite (fixed) ---------------------------- #

ICITE_BASE = "https://icite.od.nih.gov/api/pubs"
ICITE_TIMEOUT = 40

def icite_query_pmids(pmids: List[int]) -> List[Dict[str, Any]]:
    """
    Correct, robust call:
      /api/pubs?pmids=...&legacy=false&refs=true&fl=pmid,citedPmids,citedByPmids
    Force-return link fields to avoid records where they're omitted or null.
    """
    if not pmids:
        return []
    params = {
        "pmids": ",".join(str(p) for p in pmids),
        "legacy": "false",         # modern field names
        "refs": "true",            # ensure references are included
        "fl": "pmid,citedPmids,citedByPmids",  # force link fields to be present
        "format": "json",
    }
    r = requests.get(ICITE_BASE, params=params, timeout=ICITE_TIMEOUT)
    r.raise_for_status()
    data = r.json()
    if isinstance(data, dict) and "data" in data:
        return data["data"]
    if isinstance(data, list):
        return data
    if isinstance(data, dict) and "pmid" in data:
        return [data]
    return []

def icite_get_many_modern(pmids: List[int]) -> Dict[int, Dict[str, Any]]:
    """Batch fetch with modern fields + links."""
    out = {}
    for chunk_start in range(0, len(pmids), 1000):
        chunk = pmids[chunk_start:chunk_start+1000]
        for d in icite_query_pmids(chunk):
            if "pmid" in d:
                out[int(d["pmid"])] = d
    return out

def icite_extract_links(doc: Dict[str, Any]) -> Tuple[Set[int], Set[int]]:
    """
    Safely extract refs/citers with null/omitted-field tolerance + legacy fallback.
    """
    def _as_ids(x) -> List[int]:
        if not x:
            return []
        if isinstance(x, list):
            return [int(v) for v in x if str(v).isdigit()]
        # occasionally comes as strings or other – be conservative
        return []

    if not doc:
        return set(), set()

    # Prefer modern names
    refs = set(_as_ids(doc.get("citedPmids")))
    citers = set(_as_ids(doc.get("citedByPmids")))

    # If modern empty, try legacy names (some deployments still use them)
    if not refs and not citers:
        refs = set(_as_ids(doc.get("references")))
        citers = set(_as_ids(doc.get("cited_by")))

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

def parse_pubmed_records(root: ET.Element) -> Dict[int, Dict[str, Any]]:
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

def get_pubmed_metadata(pmids: List[int]) -> Dict[int, Dict[str, Any]]:
    res: Dict[int, Dict[str, Any]] = {}
    BATCH = 150
    for i in range(0, len(pmids), BATCH):
        root = efetch_pubmed_xml(pmids[i:i+BATCH])
        res.update(parse_pubmed_records(root))
        time.sleep(EUTILS_SLEEP)
    return res

# ---------------------------- Embeddings ---------------------------- #

def lmstudio_embed(text: str) -> List[float]:
    def embed_once(chunk: str) -> List[float]:
        payload = {"model": EMB_MODEL, "input": chunk}
        r = requests.post(LMSTUDIO_EMB_URL, json=payload, timeout=EMB_TIMEOUT)
        if r.status_code != 200:
            raise RuntimeError(f"Embeddings HTTP {r.status_code}: {r.text[:200]}")
        data = r.json()
        return [float(x) for x in data["data"][0]["embedding"]]

    chunks = chunk_text(text, EMB_MAX_CHARS)
    if not chunks: return []
    if len(chunks) == 1:
        return embed_once(chunks[0])

    acc = None
    for ch in chunks:
        vec = embed_once(ch)
        if acc is None: acc = [0.0]*len(vec)
        for i,v in enumerate(vec): acc[i] += v
    return [v/len(chunks) for v in acc]

# ---------------------------- Logic ---------------------------- #

def is_srma(pub_types: List[str]) -> bool:
    s = set(pub_types or [])
    return bool(s & SRMA_PT)

def is_primary(pub_types: List[str]) -> bool:
    s = set(pub_types or [])
    if s & PRIMARY_PT: return True
    return any("Trial" in pt for pt in s)

def join_title_abstract(m: Dict[str,Any]) -> str:
    return f"{m.get('title','').strip()}\n\n{m.get('abstract','').strip()}".strip()

def summarize(vals: List[float]) -> Dict[str, Any]:
    if not vals: return {"n": 0}
    v = sorted(vals); n = len(v)
    def q(p):
        if n==1: return v[0]
        k=(n-1)*p; f=math.floor(k); c=math.ceil(k)
        return v[f] if f==c else v[f]+(v[c]-v[f])*(k-f)
    return {"n":n, "mean":sum(v)/n, "median":q(0.5), "p10":q(0.10), "p25":q(0.25),
            "p75":q(0.75), "p90":q(0.90), "min":v[0], "max":v[-1]}
    
# ---------------------------- Logging Helpers ---------------------------- #

def pick_primary_set_for_compact(pairwise, roles):
    """Return the unique set of primary PMIDs that appear in SR–primary pairs."""
    prim_ids = {p["primary_pmid"] for p in pairwise}
    return sorted([pid for pid in prim_ids if roles.get(pid, {}).get("is_primary")])

def build_compact_payload(
    seed_pmid: int,
    seed_refs: set,
    seed_citers: set,
    emb_model: str,
    universe: set,
    sr_neighbors: dict,
    srma_pmids: list,
    meta_all: dict,
    roles: dict,
    pairwise: list,
    overall: dict,
    by_sr: dict,
    topk: int = 5
):
    # SR/MA cards
    sr_cards = []
    for spmid in srma_pmids:
        m = meta_all.get(spmid, {})
        neigh = sr_neighbors.get(spmid, {"refs": [], "citers": [], "common_with_seed": []})
        sr_cards.append({
            "pmid": spmid,
            "title": m.get("title", ""),
            "pub_types": m.get("pub_types", []),
            "refs_count": len(neigh.get("refs", [])),
            "citers_count": len(neigh.get("citers", [])),
            "common_with_seed_count": len(neigh.get("common_with_seed", [])),
        })

    # Primary list limited to those actually used in pairwise eval
    primary_pmids = pick_primary_set_for_compact(pairwise, roles)
    prim_cards = []
    for pid in primary_pmids:
        m = meta_all.get(pid, {})
        prim_cards.append({
            "pmid": pid,
            "title": m.get("title", ""),
            "pub_types": m.get("pub_types", []),
        })

    # Similarity rows, but add primary_title inline for auditability
    # (titles are convenient for reading diffs or browsing)
    sim_rows = []
    for row in sorted(pairwise, key=lambda d: d["cosine"], reverse=True):
        p_title = meta_all.get(row["primary_pmid"], {}).get("title", "")
        sim_rows.append({
            "sr_pmid": row["sr_pmid"],
            "primary_pmid": row["primary_pmid"],
            "primary_title": p_title[:300],  # keep compact
            "cosine": round(row["cosine"], 6),
            "in_common_with_seed": bool(row.get("in_common_with_seed", False))
        })

    # Top/bottom convenience lists (optional)
    top5 = sim_rows[:topk]
    bottom5 = list(reversed(sim_rows))[:topk] if sim_rows else []

    compact = {
        "seed": {
            "pmid": seed_pmid,
            "title": meta_all.get(seed_pmid, {}).get("title", ""),
            "refs_count": len(seed_refs),
            "citers_count": len(seed_citers),
        },
        "model": emb_model,
        "universe_size": len(universe),
        "srmas": sr_cards,
        "primary_articles": prim_cards,
        "similarities": sim_rows,
        "metrics": {
            "overall": {k: (round(v, 3) if isinstance(v, float) else v) for k, v in overall.items()},
            "by_sr": {str(k): {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in val.items()}
                      for k, val in by_sr.items()},
            "top5_pairs": top5,
            "bottom5_pairs": bottom5
        }
    }
    return compact

def write_compact_json(
    path: str,
    seed_pmid: int,
    seed_refs: set,
    seed_citers: set,
    emb_model: str,
    universe: set,
    sr_neighbors: dict,
    srma_pmids: list,
    meta_all: dict,
    roles: dict,
    pairwise: list,
    overall: dict,
    by_sr: dict
):
    compact = build_compact_payload(
        seed_pmid, seed_refs, seed_citers, emb_model, universe,
        sr_neighbors, srma_pmids, meta_all, roles, pairwise, overall, by_sr
    )
    with open(path, "w", encoding="utf-8") as f:
        json.dump(compact, f, ensure_ascii=False, indent=2)


def main():
    setup_logging()
    t0 = time.time()
    logging.info(f"Seed PMID: {SEED_PMID}")

    # --- iCite seed (modern fields + references) ---
    seed_doc = icite_get_many_modern([SEED_PMID]).get(SEED_PMID, {})  # new
    seed_refs, seed_citers = icite_extract_links(seed_doc)

    # Assert visibility (debug)
    logging.info(f"iCite raw keys for seed: {list(seed_doc.keys())[:12] if seed_doc else 'NONE'}")
    logging.info(f"Seed refs: {len(seed_refs)} | Seed citers: {len(seed_citers)}")
    if len(seed_refs)==0 and len(seed_citers)==0:
        logging.warning("iCite returned 0 refs and 0 citers for the seed. "
                        "This is unusual—double-check connectivity and that pmid is in OCC. "
                        "Proceeding anyway with downstream steps.")

    # --- Find SR/MAs among citers ---
    citer_list = sorted(seed_citers)
    meta_citers = get_pubmed_metadata(citer_list) if citer_list else {}
    srma_pmids = [p for p,m in meta_citers.items() if is_srma(m.get("pub_types", []))]
    logging.info(f"SR/MAs citing seed: {len(srma_pmids)}")
    
    if srma_pmids:
        first = srma_pmids[0]
        dbg = icite_get_many_modern([first]).get(first, {})
        logging.info(f"iCite raw keys for first SR {first}: {sorted(list(dbg.keys()))[:20]}")

    # --- Neighborhoods for each SR/MA ---
    sr_docs = icite_get_many_modern(srma_pmids) if srma_pmids else {}
    seed_neighbors = seed_refs | seed_citers
    sr_neighbors: Dict[int, Dict[str, Any]] = {}
    for spmid, doc in sr_docs.items():
        srefs, sciters = icite_extract_links(doc)
        common = (seed_neighbors) & (srefs | sciters)
        sr_neighbors[spmid] = {
            "refs": sorted(srefs),
            "citers": sorted(sciters),
            "common_with_seed": sorted(common),
        }

    # --- Universe ---
    universe: Set[int] = {SEED_PMID} | set(srma_pmids) | set(seed_neighbors)
    for s in sr_neighbors.values():
        universe.update(s["common_with_seed"])
        universe.update(s["refs"])  # ensure SR reference set pulled in

    logging.info(f"Universe size (seed + SR/MAs + commons + SR refs): n={len(universe)}")

    # --- PubMed metadata ---
    meta_all = get_pubmed_metadata(sorted(universe)) if universe else {}
    for spmid in srma_pmids:
        if spmid not in meta_all and spmid in meta_citers:
            meta_all[spmid] = meta_citers[spmid]

    roles = {}
    for pmid, m in meta_all.items():
        pts = m.get("pub_types", [])
        roles[pmid] = {
            "is_seed": (pmid == SEED_PMID),
            "is_srma": is_srma(pts),
            "is_primary": is_primary(pts),
        }

    # --- Embeddings ---
    embeddings: Dict[int, Dict[str, Any]] = {}
    for pmid, m in meta_all.items():
        text = join_title_abstract(m)
        try:
            vec = lmstudio_embed(text) if text else []
        except Exception as e:
            logging.error(f"Embedding failed for PMID {pmid}: {e}")
            vec = []
        unit = l2_normalize(vec) if vec else []
        embeddings[pmid] = {
            "dim": len(vec),
            "norm": l2_norm(vec) if vec else 0.0,
            "hash": hashlib.md5((",".join(f"{x:.6f}" for x in unit)).encode("utf-8")).hexdigest() if unit else None,
            "vec": unit,
        }

    # --- SR vs its referenced PRIMARY evidence (only if SR actually cites it) ---
    pairwise = []
    for spmid in srma_pmids:
        s_vec = embeddings.get(spmid, {}).get("vec") or []
        if not s_vec: continue
        s_ref_set = set(sr_neighbors.get(spmid, {}).get("refs", []))
        prims = [p for p in s_ref_set if roles.get(p, {}).get("is_primary") and p in embeddings]
        for p in prims:
            p_vec = embeddings[p]["vec"]
            if not p_vec: continue
            pairwise.append({
                "sr_pmid": spmid,
                "primary_pmid": p,
                "cosine": float(cosine_unit(s_vec, p_vec)),
                "in_common_with_seed": p in set(sr_neighbors.get(spmid, {}).get("common_with_seed", []))
            })

    # --- Metrics ---
    overall = summarize([x["cosine"] for x in pairwise])
    by_sr: Dict[int, Dict[str, Any]] = {}
    for spmid, group in itertools.groupby(sorted(pairwise, key=lambda d: d["sr_pmid"]), key=lambda d: d["sr_pmid"]):
        vals = [g["cosine"] for g in group]
        by_sr[spmid] = summarize(vals)

    # --- JSON ---
    data = {
        "seed_pmid": SEED_PMID,
        "timestamp_utc": now_iso(),
        "icite_endpoint": f"{ICITE_BASE}?pmids=...&legacy=false&refs=true",
        "icite_fields_note": "Modern fields: citedByPmids (citers), citedPmids (references). Legacy: cited_by, references.",
        "lmstudio": {"endpoint": LMSTUDIO_EMB_URL, "model": EMB_MODEL},
        "seed_neighbors": {"refs": sorted(seed_refs), "citers": sorted(seed_citers)},
        "srma_pmids": sorted(srma_pmids),
        "sr_neighbors": sr_neighbors,
        "universe_pmids": sorted(list(universe)),
        "pubmed_meta": meta_all,
        "roles": roles,
        "embeddings": embeddings,     # normalized vecs
        "sr_primary_pairs": pairwise,
        "metrics": {
            "overall_pairwise_cosine": overall,
            "by_sr": by_sr
        }
    }

    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)

    logging.info(f"Saved JSON → {OUT_JSON} | elapsed {time.time()-t0:.1f}s")
    logging.info("=== Topline ===")
    logging.info(f"SR/MAs citing seed: {len(srma_pmids)} | Universe size: {len(universe)} | SR–primary pairs: {overall.get('n',0)}")
    if overall.get("n",0):
        logging.info(f"Cosine mean={overall['mean']:.3f} | median={overall['median']:.3f} | p25={overall['p25']:.3f} | p75={overall['p75']:.3f} | min={overall['min']:.3f} | max={overall['max']:.3f}")
    else:
        logging.info("No pairs. If seed links are still zero, verify iCite access and that the PMID exists in OCC.")
    if not (len(seed_refs) or len(seed_citers)):
        logging.info("DEBUG: Try curl this by hand:\n"
                     f"curl '{ICITE_BASE}?pmids={SEED_PMID}&legacy=false&refs=true' | jq .")

    COMPACT_JSON = os.path.join(OUT_DIR, "sr_semantic_eval.compact.json")
    write_compact_json(
        COMPACT_JSON,
        SEED_PMID,
        seed_refs,
        seed_citers,
        EMB_MODEL,
        universe,
        sr_neighbors,
        srma_pmids,
        meta_all,
        roles,
        pairwise,
        overall,
        by_sr
    )
    logging.info(f"Wrote compact JSON → {COMPACT_JSON}")

if __name__ == "__main__":
    main()
