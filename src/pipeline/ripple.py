# src/pipeline/ripple.py
from __future__ import annotations
from typing import List, Dict, Any, Set, Tuple
import numpy as np, pandas as pd

from cache.icite import ICiteCache
from cache.emb import EmbCache
from clients.icite import get_pubs, extract_refs_and_citers
from clients.entrez import efetch_abstracts
from clients.lmstudio import LMEmbeddings
from themes.hybrid_graph import cosine_sim_matrix
from config import LMSTUDIO_EMB_MODEL

def ripple_expand_from_primaries(seed_pmids: List[str],
                                 allowed_since_year: int | None,
                                 max_expand: int = 300,
                                 prefer: str = "citers") -> Dict[str,Any]:
    """
    Expand around seed primary studies, preferring recent citers (or refs), with a strict cap.
    Returns fetched docs {pmid->meta}, the ordered candidate list, and STATS for auditing.
    """
    ic = ICiteCache()
    have = ic.get_many(seed_pmids, legacy=True)
    need = [p for p in seed_pmids if p not in have]
    if need:
        fetched = get_pubs(need, fields=["pmid","references","cited_by","year"], legacy=True)
        ic.put_many(fetched, legacy=True)
        for rec in fetched:
            have[str(rec.get("pmid") or rec.get("_id") or "")] = rec

    # gather pool per seed for auditing
    per_seed = []
    cand: List[int] = []
    for s in seed_pmids:
        refs,citers = extract_refs_and_citers(have.get(s, {}))
        pool = citers if prefer=="citers" else refs
        per_seed.append({
            "seed": str(s), 
            "refs": len(refs), 
            "citers": len(citers), 
            "used": "citers" if prefer=="citers" else "refs",
            "pool_size": len(pool)
        })
        cand.extend(pool)

    # unique and (optional) filter by year
    cand = list(dict.fromkeys(cand))
    n_before_year = len(cand)

    passed_year = cand
    if allowed_since_year is not None:
        have2 = ic.get_many([str(x) for x in cand], legacy=True)
        need2 = [str(x) for x in cand if str(x) not in have2]
        if need2:
            fetched2 = get_pubs(need2, fields=["pmid","year"], legacy=True)
            ic.put_many(fetched2, legacy=True)
            for rec in fetched2:
                have2[str(rec.get("pmid") or rec.get("_id") or "")] = rec
        passed_year = [c for c in cand if (have2.get(str(c),{}).get("year") or 0) >= allowed_since_year]

    # cap
    capped = passed_year[:max_expand]

    # fetch metadata for candidates
    meta = efetch_abstracts([str(x) for x in capped])

    # embed (cache)
    cache = EmbCache()
    pid = [str(x) for x in meta.keys()]
    texts = [(meta[p]["title"] or "") + "\n" + (meta[p]["abstract"] or "") for p in pid]
    cached = cache.get_many(LMSTUDIO_EMB_MODEL, pid)
    miss_idx, miss_texts, miss_pmids = [], [], []
    for i,p in enumerate(pid):
        if p not in cached:
            miss_idx.append(i); miss_texts.append(texts[i]); miss_pmids.append(p)
    if miss_texts:
        emb = LMEmbeddings()
        new_vecs = emb.encode(miss_texts, batch_size=32)
        cache.put_many(LMSTUDIO_EMB_MODEL, [(miss_pmids[j], new_vecs[j]) for j in range(len(miss_pmids))])

    stats = {
        "seeds": per_seed,
        "prefer": prefer,
        "allowed_since_year": allowed_since_year,
        "pool_unique": n_before_year,
        "after_year_gate": len(passed_year),
        "after_cap": len(capped),
        "reason_zero": None if len(capped)>0 else (
            "no-citers-or-refs" if all(s["pool_size"]==0 for s in per_seed) else
            "year-gate-too-strict" if n_before_year>0 and len(passed_year)==0 else
            "capped-to-zero" if n_before_year>0 and len(capped)==0 else
            "unknown"
        )
    }

    return {"meta": meta, "candidates": [str(x) for x in capped], "stats": stats}
