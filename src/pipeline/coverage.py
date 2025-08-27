# src/pipeline/coverage.py
from __future__ import annotations
from typing import Dict, Any, List, Tuple, Set
import numpy as np
import pandas as pd

from cache.icite import ICiteCache
from clients.icite import get_pubs, extract_refs_and_citers
from pipeline.evidence import split_by_kind
from config import COV_LEVELS

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

def coverage_for_theme(theme: Dict[str,Any], docs_df: pd.DataFrame) -> Dict[str,Any]:
    """
    Compute coverage metrics for one theme:
      - E: primaries in theme
      - S: SR/MA in theme
      - For each SR: CoveredPrimaries = included ∩ E
      - CoverageRatio = |∪covered| / |E|
      - NewPrimaryCount since max SR year
    """
    members = theme["members_pmids"]
    sub = docs_df[docs_df["pmid"].astype(str).isin(members)].copy()
    prim, sr, _ = split_by_kind(sub.to_dict(orient="records"))
    E = set(prim)
    S = list(sr)
    if not E:
        return {"theme_id": theme["theme_id"], "E_size": 0, "S_count": len(S),
                "coverage_ratio": 0.0, "covered": [], "sr_map": {}, "new_primary_count": 0,
                "last_sr_year": None, "coverage_level": "NONE"}
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
    return {
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
