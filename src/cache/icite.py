# src/cache/icite.py
from __future__ import annotations
import sqlite3, pathlib, json
from typing import Iterable, Dict, Any

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[2]
CACHE_DIR = PROJECT_ROOT / "data" / "cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
DB_PATH = CACHE_DIR / "icite.sqlite3"

class ICiteCache:
    """Cache for iCite /pubs responses keyed by PMID and legacy flag."""
    def __init__(self, db_path: pathlib.Path = DB_PATH):
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
