# src/clients/icite.py
from __future__ import annotations
import pathlib, sqlite3, json, re, requests
from typing import List, Dict, Any, Iterable, Tuple
from config import ICITE_BASE, HTTP_TIMEOUT, USER_AGENT, ENTREZ_EMAIL, ENTREZ_API_KEY

# -------------------------
# Cache (unchanged)
# -------------------------
PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[2]
CACHE_DIR = PROJECT_ROOT / "data" / "cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
DB_PATH = CACHE_DIR / "icite.sqlite3"

class ICiteCache:
    """
    Cache for iCite /pubs responses.
    Key: pmid (TEXT). We store the JSON blob and a 'legacy' flag (0/1).
    """
    def __init__(self, db_path: pathlib.Path = DB_PATH):
        self.db_path = db_path
        self._conn = sqlite3.connect(str(db_path))
        self._conn.execute("""
            CREATE TABLE IF NOT EXISTS pubs(
                pmid TEXT PRIMARY KEY,
                legacy INTEGER NOT NULL,
                json TEXT NOT NULL
            )
        """)
        self._conn.commit()

    def get_many(self, pmids: Iterable[str], legacy: bool = True) -> Dict[str, Dict[str,Any]]:
        pmids = [str(p) for p in pmids]
        out: Dict[str, Dict[str,Any]] = {}
        if not pmids: return out
        qmarks = ",".join(["?"]*len(pmids))
        cur = self._conn.execute(
            f"SELECT pmid, json FROM pubs WHERE legacy=? AND pmid IN ({qmarks})",
            [1 if legacy else 0] + pmids
        )
        for pmid, blob in cur.fetchall():
            try:
                out[pmid] = json.loads(blob)
            except Exception:
                pass
        return out

    def put_many(self, rows: Iterable[Dict[str,Any]], legacy: bool = True) -> int:
        data = []
        for rec in rows:
            pmid = str(rec.get("pmid") or rec.get("_id") or "")
            if not pmid:
                continue
            data.append((pmid, 1 if legacy else 0, json.dumps(rec)))
        if not data: return 0
        self._conn.executemany("INSERT OR REPLACE INTO pubs(pmid,legacy,json) VALUES(?,?,?)", data)
        self._conn.commit()
        return len(data)

    def close(self):
        try: self._conn.close()
        except Exception: pass

# -------------------------
# iCite HTTP helpers
# -------------------------
HEADERS = {"User-Agent": USER_AGENT, "Accept": "application/json"}

def _icite_url(path: str) -> str:
    """
    Build a correct iCite URL whether ICITE_BASE ends with '/api' or not.
    """
    base = ICITE_BASE.rstrip('/')
    # If the base already ends with '/api', don't add another '/api'
    if base.endswith('/api'):
        return f"{base}/{path.lstrip('/')}"
    # Otherwise, append '/api'
    return f"{base}/api/{path.lstrip('/')}"

def _icite_pubs_fetch(pmids: List[str]) -> List[Dict[str,Any]]:
    """Fetch raw iCite /pubs JSON (batched). Returns list of records."""
    pmids = [str(p) for p in pmids]
    out: List[Dict[str,Any]] = []
    url = _icite_url("pubs")
    for i in range(0, len(pmids), 200):
        sub = pmids[i:i+200]
        params = {"pmids": ",".join(sub), "format": "json"}  # format=json is important
        r = requests.get(url, headers=HEADERS, params=params, timeout=HTTP_TIMEOUT)
        if r.status_code >= 400:
            # Helpful debug on auth/routing mistakes
            raise requests.HTTPError(f"{r.status_code} {r.reason} for url: {r.url}\nBody: {r.text[:300]}")
        data = r.json()
        rows = data.get("data") if isinstance(data, dict) else data
        if not isinstance(rows, list): rows = []
        out.extend(rows)
    return out

# -------- Entrez ELink fallback (unchanged interface) --------
_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
def _elink(pmids: List[str], linkname: str) -> Dict[str, List[str]]:
    pmids = [str(p) for p in pmids]
    out: Dict[str,List[str]] = {p: [] for p in pmids}
    for i in range(0, len(pmids), 200):
        sub = pmids[i:i+200]
        params = {
            "dbfrom": "pubmed", "db": "pubmed",
            "linkname": linkname,
            "id": ",".join(sub),
            "retmode": "json",
            "email": ENTREZ_EMAIL
        }
        if ENTREZ_API_KEY: params["api_key"] = ENTREZ_API_KEY
        r = requests.get(_EUTILS, headers=HEADERS, params=params, timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        js = r.json()

        for ls in (js.get("linksets", []) or []):
            ids_list = ls.get("ids") or []
            uid = str(ids_list[0]) if ids_list else ""

            linkdbs = ls.get("linksetdbs") or []
            if isinstance(linkdbs, dict):
                linkdbs = [linkdbs]

            for grp in linkdbs:
                if not isinstance(grp, dict):
                    continue
                if grp.get("linkname") != linkname:
                    continue

                ids: List[str] = []

                # links can be a list of dicts ({"id": "123"}), a list of strings, or absent
                links = grp.get("links")
                if isinstance(links, list) and links:
                    for x in links:
                        if isinstance(x, dict):
                            v = x.get("id") or x.get("Id")
                            if v: ids.append(str(v))
                        else:
                            ids.append(str(x))

                # some payloads expose only "ids": [...]
                if not ids:
                    only_ids = grp.get("ids")
                    if isinstance(only_ids, list):
                        ids.extend(str(z) for z in only_ids if z)

                if uid and ids:
                    out.setdefault(uid, []).extend(ids)

    # dedupe+sort for stability
    return {k: sorted(set(v)) for k,v in out.items()}

def _elink_refs(pmids: List[str]) -> Dict[str, List[str]]:
    return _elink(pmids, "pubmed_pubmed_refs")

def _elink_citers(pmids: List[str]) -> Dict[str, List[str]]:
    return _elink(pmids, "pubmed_pubmed_citedin")


# -------- Public API --------
_DIGITS = re.compile(r"\d+")
def _coerce_pmid_list(val) -> List[int]:
    if val is None: return []
    if isinstance(val, list):
        flat = []
        for x in val:
            if isinstance(x, dict) and "pmid" in x: x = x["pmid"]
            if isinstance(x, (int, float)) or (isinstance(x, str) and x.isdigit()):
                flat.append(int(x))
            elif isinstance(x, str):
                flat.extend([int(m.group()) for m in _DIGITS.finditer(x)])
        return sorted(set(flat))
    if isinstance(val, str):
        return sorted(set(int(m.group()) for m in _DIGITS.finditer(val)))
    return []

def get_pubs(pmids: Iterable[str|int],
             fields: List[str] | None = None,
             legacy: bool = True) -> List[Dict[str,Any]]:
    ids = [str(p) for p in pmids]
    if not ids: return []

    raw = _icite_pubs_fetch(ids)

    # normalize (keep explicit empties so we can distinguish "present but empty" from "missing")
    norm: Dict[str, Dict[str,Any]] = {}
    had_refs_key: set[str] = set()
    had_citers_key: set[str] = set()

    for rec in raw:
        pmid = str(rec.get("pmid") or rec.get("_id") or "").strip()
        if not pmid: 
            continue
        norm.setdefault(pmid, {})
        norm[pmid].update(rec)

        if "references" in rec:
            had_refs_key.add(pmid)
            norm[pmid]["references"] = _coerce_pmid_list(rec.get("references"))

        if "cited_by" in rec:
            had_citers_key.add(pmid)
            norm[pmid]["cited_by"] = _coerce_pmid_list(rec.get("cited_by"))

    # -------- Fallbacks only for PMIDs where the field was truly missing in iCite --------
    missing_refs   = [p for p in ids if p not in had_refs_key]
    if missing_refs:
        refs_map = _elink_refs(missing_refs)
        for p in missing_refs:
            if refs_map.get(p) is not None:
                norm.setdefault(p, {})["references"] = [int(x) for x in refs_map.get(p, [])]

    missing_citers = [p for p in ids if p not in had_citers_key]
    if missing_citers:
        cit_map = _elink_citers(missing_citers)
        for p in missing_citers:
            if cit_map.get(p) is not None:
                norm.setdefault(p, {})["cited_by"] = [int(x) for x in cit_map.get(p, [])]

    # -------- Build output (preserve requested fields) --------
    out: List[Dict[str,Any]] = []
    for p in ids:
        rec = norm.get(p, {"pmid": p})
        if fields:
            out.append({k: rec.get(k) for k in fields if k in rec or k == "pmid"})
            out[-1]["pmid"] = p
        else:
            out.append(rec)
    return out


def extract_refs_and_citers(rec: Dict[str,Any]) -> Tuple[List[int], List[int]]:
    refs = _coerce_pmid_list(rec.get("references"))
    cit  = _coerce_pmid_list(rec.get("cited_by"))
    return refs, cit