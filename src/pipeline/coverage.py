# src/pipeline/coverage.py
from __future__ import annotations
from typing import Dict, Any, List, Tuple, Set
import numpy as np
import pandas as pd
import itertools

from cache.icite import ICiteCache
from clients.icite import get_pubs, extract_refs_and_citers
from pipeline.evidence import split_by_kind
from config import COV_LEVELS, LMSTUDIO_EMB_MODEL
from cache.emb import EmbCache

def sr_included_primaries(sr_pmids: List[str],
                          primary_pool: Set[str]) -> Dict[str, Set[str]]:
    """
    Approximate SR 'included studies' as its referenced PMIDs intersected with known primary pool.
    """
    ic = ICiteCache()
    have = ic.get_many(sr_pmids, legacy=True)
    need = [p for p in sr_pmids if p not in have]
    if need:
        fetched = get_pubs(need, fields=["pmid","references"], legacy=True)
        ic.put_many(fetched, legacy=True)
        for rec in fetched:
            have[str(rec.get("pmid") or rec.get("_id") or "")] = rec
    out: Dict[str,Set[str]] = {}
    for s in sr_pmids:
        refs,_ = extract_refs_and_citers(have.get(s, {}))
        out[s] = set(str(x) for x in refs) & primary_pool
    return out

def _doc_map(docs_df: pd.DataFrame, pmids: List[str], abstract_chars: int = 420) -> Dict[str, Any]:
    sub = docs_df.set_index(docs_df["pmid"].astype(str)).loc[pmids]
    out = {}
    for pid, row in sub.iterrows():
        abst = (row.get("abstract") or "")
        out[str(pid)] = {
            "pmid": str(pid),
            "title": row.get("title") or "",
            "abstract": abst[:abstract_chars] + ("…" if len(abst) > abstract_chars else ""),
            "year": int(row.get("year")) if pd.notna(row.get("year")) else None,
            "doi": (row.get("doi") or None),
            "journal": (row.get("journal") or None),
            "pub_types": list(row.get("pub_types") or []),
        }
    return out

def _semantic_stats(theme: Dict[str,Any], pmids: List[str]) -> Dict[str, Any]:
    """
    Compute cosine to theme centroid for each PMID (if embedding cached).
    Returns {pmid: cos}, plus aggregate stats (mean/median).
    """
    cent = np.array(theme.get("centroid") or [], dtype="float32")
    if cent.size == 0:
        return {"per_doc": {}, "mean": None, "median": None, "n": 0}
    cent = cent / (np.linalg.norm(cent) + 1e-12)

    cache = EmbCache()
    vecs = cache.get_many(LMSTUDIO_EMB_MODEL, pmids)
    cos = {}
    vals = []
    for pid in pmids:
        v = vecs.get(pid)
        if v is None: 
            continue
        v = v / (np.linalg.norm(v) + 1e-12)
        c = float(v.dot(cent))
        cos[pid] = c
        vals.append(c)
    if not vals:
        return {"per_doc": cos, "mean": None, "median": None, "n": 0}
    arr = np.array(vals, dtype="float32")
    return {"per_doc": cos, "mean": float(arr.mean()), "median": float(np.median(arr)), "n": int(arr.size)}

def _bib_coupling_stats(pmids: List[str], max_pairs: int = 20000) -> Dict[str, Any]:
    """
    Pairwise Jaccard over references on the set (sampled if large).
    """
    m = len(pmids)
    if m < 2:
        return {"pairs": 0, "mean": None, "p90": None}
    ic = ICiteCache()
    have = ic.get_many(pmids, legacy=True)
    need = [p for p in pmids if p not in have]
    if need:
        fetched = get_pubs(need, fields=["pmid","references"], legacy=True)
        ic.put_many(fetched, legacy=True)
        for rec in fetched:
            have[str(rec.get("pmid") or rec.get("_id") or "")] = rec

    refsets = {p: set(extract_refs_and_citers(have.get(p, {}))[0]) for p in pmids}
    pairs = list(itertools.combinations(pmids, 2))
    if len(pairs) > max_pairs:
        # uniform sample
        rng = np.random.default_rng(17)
        pairs = list(rng.choice(pairs, size=max_pairs, replace=False))

    jvals = []
    for a,b in pairs:
        Ra, Rb = refsets.get(a, set()), refsets.get(b, set())
        inter = len(Ra & Rb)
        if inter == 0:
            j = 0.0
        else:
            uni = len(Ra | Rb) or 1
            j = inter / uni
        jvals.append(j)
    if not jvals:
        return {"pairs": 0, "mean": None, "p90": None}
    arr = np.array(jvals, dtype="float32")
    return {"pairs": int(arr.size), "mean": float(arr.mean()), "p90": float(np.percentile(arr, 90.0))}

def coverage_for_theme(theme: Dict[str,Any], docs_df: pd.DataFrame, include_docs: bool = True) -> Dict[str,Any]:
    """
    Compute coverage metrics for one theme:
      - E: primaries in theme
      - S: SR/MA in theme
      - SR → included primaries (approx via references)
      - CoverageRatio = |∪covered| / |E|
      - NewPrimaryCount since max SR year
      - Optional: attach doc details + semantic & bibliographic stats for auditing
    """
    members = theme["members_pmids"]
    sub = docs_df[docs_df["pmid"].astype(str).isin(members)].copy()
    prim, sr, _ = split_by_kind(sub.to_dict(orient="records"))
    E = set(prim)
    S = list(sr)

    # empty E quick return
    if not E:
        row = {"theme_id": theme["theme_id"], "E_size": 0, "S_count": len(S),
               "coverage_ratio": 0.0, "covered": [], "sr_map": {}, "new_primary_count": 0,
               "last_sr_year": None, "coverage_level": "NONE", "E": [], "S": []}
        if include_docs:
            row["details"] = {
                "E_docs": {}, "S_docs": {}, "sem": {}, "bib": {},
                "sr_breakdown": []
            }
        return row

    # SR → included primaries (approx by references)
    sr_map = sr_included_primaries(S, E) if S else {}
    covered_union: Set[str] = set()
    for s in S:
        covered_union |= sr_map.get(s, set())
    cov_ratio = (len(covered_union) / len(E)) if E else 0.0

    # recency: SR "last search" proxy = max(SR year)
    years = pd.to_numeric(sub.set_index("pmid").loc[S]["year"], errors="coerce") if S else pd.Series([], dtype=float)
    last_sr_year = int(years.max()) if not years.empty and years.notna().any() else None

    # new primary count
    if last_sr_year is None:
        new_prim = len(E)  # treat as all new (no SR exists)
    else:
        p_years = pd.to_numeric(sub.set_index("pmid").loc[list(E)]["year"], errors="coerce")
        new_prim = int((p_years > last_sr_year).sum()) if p_years.notna().any() else 0

    # coverage level
    level = "NONE"
    for name, thr in COV_LEVELS.items():
        if cov_ratio < thr:
            level = name
            break

    row = {
        "theme_id": theme["theme_id"],
        "E_size": len(E),
        "S_count": len(S),
        "coverage_ratio": cov_ratio,
        "covered": sorted(list(covered_union)),
        "sr_map": {k: sorted(list(v)) for k,v in sr_map.items()},
        "new_primary_count": new_prim,
        "last_sr_year": last_sr_year,
        "coverage_level": level,
        "E": sorted(list(E)),
        "S": sorted(list(S)),
    }

    if include_docs:
        E_docs = _doc_map(docs_df, sorted(list(E)))
        S_docs = _doc_map(docs_df, sorted(list(S)))
        # semantic alignments to theme centroid
        sem_E = _semantic_stats(theme, list(E_docs.keys()))
        sem_S = _semantic_stats(theme, list(S_docs.keys()))
        # biblio coupling among E
        bib_E = _bib_coupling_stats(sorted(list(E)), max_pairs=20000)
        # SR breakdown rows
        sr_breakdown = []
        for s in S:
            sdoc = S_docs.get(s, {})
            inc = sorted(list(sr_map.get(s, set())))
            sr_breakdown.append({
                "pmid": s,
                "year": sdoc.get("year"),
                "title": sdoc.get("title"),
                "included_primaries": inc,
                "coverage_of_E": float(len(inc) / max(1, len(E)))
            })
        row["details"] = {
            "E_docs": E_docs,
            "S_docs": S_docs,
            "sem": {"E": sem_E, "S": sem_S},
            "bib": {"E_jaccard_refs": bib_E},
            "sr_breakdown": sr_breakdown
        }
    return row
