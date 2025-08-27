# src/utils/io.py
from __future__ import annotations
import pathlib, orjson
from typing import Any

try:
    import numpy as np  # optional, only to recognize numpy types
except Exception:
    np = None

def _default(o):
    if np is not None:
        if isinstance(o, (np.integer, np.floating)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
    if isinstance(o, set):
        return list(o)
    if isinstance(o, pathlib.Path):
        return str(o)
    # Let orjson raise for anything unexpected
    raise TypeError(f"Type not serializable: {type(o)}")

def jdump(obj: Any, path: pathlib.Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        f.write(orjson.dumps(
            obj,
            option=orjson.OPT_INDENT_2 | orjson.OPT_SERIALIZE_NUMPY,
            default=_default
        ))
