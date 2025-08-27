# src/cache/emb.py
from __future__ import annotations
import pathlib, numpy as np
from typing import Dict, Iterable, List, Tuple

ROOT = pathlib.Path(__file__).resolve().parents[2]
BASE = ROOT / "data" / "cache" / "emb"
BASE.mkdir(parents=True, exist_ok=True)

class EmbCache:
    """
    Very simple file cache:
      data/cache/emb/{model}/{id}.npy  (float32 vector)
    """
    def __init__(self, base: pathlib.Path = BASE):
        self.base = base

    def _model_dir(self, model: str) -> pathlib.Path:
        d = self.base / model.replace("/", "_")
        d.mkdir(parents=True, exist_ok=True)
        return d

    def get_many(self, model: str, ids: Iterable[str]) -> Dict[str, np.ndarray]:
        d = self._model_dir(model)
        out: Dict[str, np.ndarray] = {}
        for pid in ids:
            p = d / f"{str(pid)}.npy"
            if p.exists():
                try:
                    out[str(pid)] = np.load(p)
                except Exception:
                    pass
        return out

    def put_many(self, model: str, rows: List[Tuple[str, np.ndarray]]) -> int:
        d = self._model_dir(model)
        n = 0
        for pid, vec in rows:
            p = d / f"{str(pid)}.npy"
            try:
                np.save(p, vec.astype("float32"))
                n += 1
            except Exception:
                pass
        return n
