# src/pipeline/coverage.py
from __future__ import annotations
from typing import Dict, Any, List, Tuple, Set
import numpy as np
import pandas as pd
import itertools

from cache.icite import ICiteCache
from clients.icite import get_pubs, extract_refs_and_citers
from clients.entrez import efetch_abstracts  # NEW
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
    """
    Map pmid -> {title, abstract(snippet), year, doi, journal, pub_types}.
    If some pmids are not present in docs_df, hydrate them via efetch to avoid KeyError.
    """
    pmids = [str(p) for p in pmids]
    idx = docs_df["pmid"].astype(str)
    present_set = set(idx.tolist())
    present = [p for p in pmids if p in present_set]
    missing = [p for p in pmids if p not in present_set]

    out: Dict[str,Any] = {}
    if present:
        df_idxed = docs_df.set_index(idx)
        for pid in present:
            row = df_idxed.loc[pid]
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

    if missing:
        # Hydrate missing docs minimally for display/audit
        fetched = efetch_abstracts(missing)
        for pid, rec in fetched.items():
            abst = (rec.get("abstract") or "")
            out[str(pid)] = {
                "pmid": str(pid),
                "title": rec.get("title") or "",
                "abstract": abst[:abstract_chars] + ("…" if len(abst) > abstract_chars else ""),
                "year": rec.get("year"),
                "doi": rec.get("doi"),
                "journal": rec.get("journal"),
                "pub_types": list(rec.get("pub_types") or []),
            }
        # For any still-missing, insert minimal placeholders
        for pid in missing:
            if str(pid) not in out:
                out[str(pid)] = {"pmid": str(pid), "title": "", "abstract": "", "year": None,
                                 "doi": None, "journal": None, "pub_types": []}
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

def _harvest_sr_from_E_citers(E: List[str]) -> List[str]:
    """
    Pull SR/Review PMIDs that CITE any primary in E (across PubMed, even if outside the theme).
    """
    if not E: return []
    # 1) get citers for E
    recs = get_pubs(E, fields=["pmid","cited_by"], legacy=True)
    citers = sorted({str(c) for r in recs for c in (r.get("cited_by") or [])})
    if not citers: return []
    # 2) fetch pub types/years to filter SR/MA
    metas = get_pubs(citers, fields=["pmid","pub_types","year","title"], legacy=True)
    out = []
    for m in metas:
        pts = [p.lower() for p in (m.get("pub_types") or [])]
        txt = " ".join(pts)
        if ("systematic review" in txt) or ("meta-analysis" in txt) or ("review" in txt):
            out.append(str(m["pmid"]))
    return sorted(set(out))

def coverage_for_theme(theme: Dict[str,Any], docs_df: pd.DataFrame, include_docs: bool = True) -> Dict[str,Any]:
    """
    Compute coverage metrics for one theme:
      - E: primaries in theme
      - S: SR/MA in theme  ∪  SR/MA that cite E (outside theme allowed)
      - SR → included primaries (approx via references)
      - CoverageRatio = |∪covered| / |E|
      - NewPrimaryCount since max SR year
      - Optional: attach doc details + semantic & bibliographic stats for auditing
    """
    members = theme["members_pmids"]
    sub = docs_df[docs_df["pmid"].astype(str).isin(members)].copy()
    prim, sr, _ = split_by_kind(sub.to_dict(orient="records"))
    E = set(prim)

    # empty E quick return
    if not E:
        row = {"theme_id": theme["theme_id"], "E_size": 0, "S_count": 0,
               "coverage_ratio": 0.0, "covered": [], "sr_map": {}, "new_primary_count": 0,
               "last_sr_year": None, "coverage_level": "NONE", "E": [], "S": []}
        if include_docs:
            row["details"] = {
                "E_docs": {}, "S_docs": {}, "sem": {}, "bib": {},
                "sr_breakdown": []
            }
        return row

    # SR candidates:
    #   (A) SR/Review papers inside the theme
    S_in_theme = list(sr)
    #   (B) SRs that cite E anywhere in PubMed (brings in 40450020, etc.)
    S_from_citers = _harvest_sr_from_E_citers(sorted(list(E)))
    S_all = sorted(set(S_in_theme) | set(S_from_citers))

    # SR → included primaries (approx by references)
    sr_map = sr_included_primaries(S_all, E) if S_all else {}
    covered_union: Set[str] = set()
    for s in S_all:
        covered_union |= sr_map.get(s, set())
    cov_ratio = (len(covered_union) / len(E)) if E else 0.0

    # recency: SR "last search" proxy = max(SR year)
    if S_all:
        df_idx = docs_df.set_index(docs_df["pmid"].astype(str))
        present = [s for s in S_all if s in df_idx.index]
        years_present = pd.to_numeric(
            df_idx.loc[present, "year"], errors="coerce"
        ) if len(present) > 0 else pd.Series([], dtype=float)

        missing = [s for s in S_all if s not in df_idx.index]
        extra = efetch_abstracts(missing) if len(missing) > 0 else {}

        sr_years: List[int] = []
        if len(present) > 0:
            sr_years.extend([int(y) for y in years_present.dropna().astype(int).tolist()])
        for s in missing:
            y = extra.get(s, {}).get("year")
            if y is not None:
                sr_years.append(int(y))

        last_sr_year = max(sr_years) if sr_years else None
    else:
        last_sr_year = None

    # new primary count since last SR year
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
        "S_count": len(S_all),
        "coverage_ratio": cov_ratio,
        "covered": sorted(list(covered_union)),
        "sr_map": {k: sorted(list(v)) for k,v in sr_map.items()},
        "new_primary_count": new_prim,
        "last_sr_year": last_sr_year,
        "coverage_level": level,
        "E": sorted(list(E)),
        "S": sorted(list(S_all)),
    }

    if include_docs:
        E_docs = _doc_map(docs_df, sorted(list(E)))
        S_docs = _doc_map(docs_df, sorted(list(S_all)))
        # semantic alignments to theme centroid
        sem_E = _semantic_stats(theme, list(E_docs.keys()))
        sem_S = _semantic_stats(theme, list(S_docs.keys()))
        # biblio coupling among E
        bib_E = _bib_coupling_stats(sorted(list(E)), max_pairs=20000)
        # SR breakdown rows
        sr_breakdown = []
        for s in S_all:
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
