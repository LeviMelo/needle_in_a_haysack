# src/pipeline/gap.py
from __future__ import annotations
from typing import Dict, Any, List
import numpy as np
import pandas as pd
import re
from sklearn.feature_extraction.text import ENGLISH_STOP_WORDS

STOP = ENGLISH_STOP_WORDS
TOK = re.compile(r"[A-Za-z0-9]+")

def top_terms(titles: List[str], k: int = 8) -> List[str]:
    tf = {}
    for t in titles:
        for w in TOK.findall((t or "").lower()):
            if len(w) < 3 or w in STOP: continue
            tf[w] = tf.get(w, 0) + 1
    return [w for w,_ in sorted(tf.items(), key=lambda x: x[1], reverse=True)[:k]]

def simple_questions(theme_title_terms: List[str]) -> List[str]:
    # produce a couple of templated question sketches from term list
    if not theme_title_terms: return []
    t = theme_title_terms[:4]
    out = []
    if len(t) >= 2:
        out.append(f"In {t[0]} patients, does {t[1]} improve outcomes vs standard care?")
    if len(t) >= 3:
        out.append(f"Does {t[0]} {t[1]} reduce {t[2]} compared with usual practice?")
    if len(t) >= 4:
        out.append(f"What is the effect of {t[0]} {t[1]} on {t[2]} in {t[3]} settings?")
    return out

def gap_score(coverage_ratio: float, new_primary_count: int, E_size: int, last_sr_year: int | None, now_year: int) -> float:
    """
    Deterministic ranking: higher when coverage is low, new primaries exist, and E has mass.
    """
    cov_term = 1.0 - coverage_ratio
    recency = 0.0 if last_sr_year is None else max(0.0, min(1.0, (now_year - last_sr_year) / 6.0))  # 6y horizon
    mass = np.tanh(E_size / 30.0)  # saturate after ~30
    newp = np.tanh(new_primary_count / 10.0)
    return 0.5*cov_term + 0.2*recency + 0.2*newp + 0.1*mass

def rank_gaps(universe: Dict[str,Any], coverage_rows: List[Dict[str,Any]], now_year: int) -> List[Dict[str,Any]]:
    df = pd.DataFrame(universe["docs"])
    theme_by_id = {t["theme_id"]: t for t in universe["themes"]}
    rows = []
    for row in coverage_rows:
        tid = row["theme_id"]; t = theme_by_id[tid]
        members_idx = t["members_idx"]
        titles = [df.iloc[i]["title"] for i in members_idx]
        terms = top_terms(titles, k=8)
        qs = simple_questions(terms)
        score = gap_score(row["coverage_ratio"], row["new_primary_count"], row["E_size"], row["last_sr_year"], now_year)
        rows.append({
            "theme_id": tid,
            "gap_score": float(score),
            "coverage_ratio": row["coverage_ratio"],
            "coverage_level": row["coverage_level"],
            "E_size": row["E_size"],
            "new_primary_count": row["new_primary_count"],
            "last_sr_year": row["last_sr_year"],
            "terms": terms,
            "questions": qs,
            "E": row.get("E",[]),
            "S": row.get("S",[]),
        })
    rows.sort(key=lambda x: x["gap_score"], reverse=True)
    return rows
