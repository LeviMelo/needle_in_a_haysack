# src/clients/icite.py
from __future__ import annotations
from typing import List, Dict, Any, Iterable, Tuple
import requests
from config import ICITE_BASE, HTTP_TIMEOUT, USER_AGENT

HEADERS = {"User-Agent": USER_AGENT, "Accept": "application/json"}

def get_pubs(pmids: Iterable[int | str],
             fields: List[str] | None = None,
             legacy: bool = True) -> List[Dict[str, Any]]:
    pmids = [str(p) for p in pmids]
    out: List[Dict[str, Any]] = []
    if not pmids:
        return out
    for i in range(0, len(pmids), 800):
        batch = pmids[i:i+800]
        params = {"pmids": ",".join(batch)}
        if fields:
            params["fl"] = ",".join(fields)
        params["legacy"] = "true" if legacy else "false"
        r = requests.get(f"{ICITE_BASE}/pubs", headers=HEADERS, params=params, timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        data = r.json()
        if isinstance(data, list):
            out.extend(data)
        elif isinstance(data, dict):
            out.append(data)
    return out

def extract_refs_and_citers(rec: Dict[str, Any]) -> Tuple[List[int], List[int]]:
    refs = rec.get("citedPmids") or rec.get("references") or []
    citers = rec.get("citedByPmids") or rec.get("cited_by") or []
    def _norm(x):
        try: return int(x)
        except Exception: return None
    refs = [v for v in (_norm(v) for v in refs) if v]
    citers = [v for v in (_norm(v) for v in citers) if v]
    return refs, citers
