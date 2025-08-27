# src/pipeline/evidence.py
from __future__ import annotations
import re
from typing import List, Dict, Any, Tuple

PRIM_HINTS = {"Randomized Controlled Trial","Clinical Trial","Trial","Cohort","Case-Control","Observational Study","Prospective Studies"}
SR_HINTS   = {"Systematic Review","Meta-Analysis","Review"}

_re_rct = re.compile(r"\b(randomi[sz]ed|rct)\b", re.I)
_re_trial = re.compile(r"\b(trial|phase\s*[IiVv]+)\b", re.I)
_re_cohort = re.compile(r"\bcohort\b", re.I)
_re_case_control = re.compile(r"\bcase[- ]control\b", re.I)
_re_meta = re.compile(r"\bmeta-?analysis\b", re.I)
_re_syst = re.compile(r"\bsystematic review\b", re.I)
_re_review = re.compile(r"\breview\b", re.I)

def paper_kind(title: str, pub_types: List[str]) -> str:
    s = set(pub_types or [])
    t = title or ""
    # Strong publication type first
    if s & SR_HINTS: return "sr"
    if s & PRIM_HINTS: return "primary"
    # Fallbacks from title
    if _re_meta.search(t) or _re_syst.search(t): return "sr"
    if _re_rct.search(t) or _re_trial.search(t) or _re_cohort.search(t) or _re_case_control.search(t):
        return "primary"
    # default
    return "other"

def split_by_kind(docs: List[Dict[str,Any]]) -> Tuple[List[str], List[str], List[str]]:
    """
    Return (primary_pmids, sr_pmids, other_pmids)
    """
    prim, sr, oth = [], [], []
    for d in docs:
        k = paper_kind(d.get("title",""), d.get("pub_types",[]))
        pid = str(d.get("pmid"))
        if k=="primary": prim.append(pid)
        elif k=="sr": sr.append(pid)
        else: oth.append(pid)
    return prim, sr, oth
