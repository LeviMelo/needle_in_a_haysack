# src/clients/lmstudio.py
from __future__ import annotations
import os, subprocess, shutil
import requests, numpy as np
from typing import List, Optional

from config import LMSTUDIO_BASE, LMSTUDIO_EMB_MODEL, LMSTUDIO_CHAT_MODEL, HTTP_TIMEOUT, USER_AGENT

HEADERS_JSON = {"Content-Type": "application/json", "User-Agent": USER_AGENT}

# Prefer LM Studio Python SDK if available (pip install lmstudio)
try:
    import lmstudio as lms  # official SDK
    _HAVE_LMSDK = True
except Exception:
    _HAVE_LMSDK = False

_USE_CLI = os.getenv("LMSTUDIO_USE_CLI", "0") == "1"
_LMS = shutil.which("lms")  # LM Studio CLI path, if present

def _cli_unload_all() -> None:
    if _LMS:
        try:
            subprocess.run([_LMS, "unload", "--all"], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except Exception:
            pass

def _cli_load(model_key: str, ttl: int = 900) -> None:
    if not _LMS:
        raise RuntimeError("LM Studio CLI 'lms' not found. Either install SDK (`pip install lmstudio`) "
                           "or make the 'lms' CLI available and set LMSTUDIO_USE_CLI=1.")
    cmd = [_LMS, "load", model_key, "--ttl", str(ttl)]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

class LMEmbeddings:
    """
    Embedding client with on-demand model management:
      - SDK path: load (with TTL), embed, explicit unload()
      - REST path: optional CLI auto load/unload around the call
    Returns L2-normalized float32 ndarray [N, D].
    """
    def __init__(self,
                 base: str = LMSTUDIO_BASE,
                 model: str = LMSTUDIO_EMB_MODEL,
                 ttl_sec: int = 900):
        self.base = base.rstrip("/")
        self.model = model
        self.ttl_sec = ttl_sec

    def _encode_sdk(self, texts: List[str], batch_size: int = 64) -> np.ndarray:
        mdl = lms.embedding_model(self.model, ttl=self.ttl_sec)  # auto-load w/ TTL
        vecs: list[list[float]] = []
        # SDK exposes per-text .embed(); keep batches small to tame VRAM spikes
        for i in range(0, len(texts), batch_size):
            for t in texts[i:i+batch_size]:
                vecs.append(mdl.embed(t))
        try:
            mdl.unload()  # free VRAM immediately
        except Exception:
            pass
        arr = np.array(vecs, dtype="float32")
        arr /= (np.linalg.norm(arr, axis=1, keepdims=True) + 1e-12)
        return arr

    def _encode_rest(self, texts: List[str], batch_size: int = 128) -> np.ndarray:
        # Optionally ensure only the embedding model is loaded
        if _USE_CLI:
            _cli_unload_all()
            _cli_load(self.model, ttl=self.ttl_sec)

        vecs: list[list[float]] = []
        url = f"{self.base}/v1/embeddings"
        for i in range(0, len(texts), batch_size):
            body = {"model": self.model, "input": texts[i:i+batch_size]}
            r = requests.post(url, headers=HEADERS_JSON, json=body, timeout=HTTP_TIMEOUT)
            if r.status_code != 200:
                msg = r.text
                if "model_not_found" in msg or "Failed to load model" in msg:
                    hint = ("Embeddings model isn't loaded. Install the lmstudio SDK (preferred) "
                            "or set LMSTUDIO_USE_CLI=1 so we auto load/unload via CLI.")
                    raise RuntimeError(f"LM Studio embeddings error: {msg}\n{hint}")
                r.raise_for_status()
            data = r.json()
            vecs.extend(d["embedding"] for d in data["data"])

        if _USE_CLI:
            _cli_unload_all()  # free VRAM

        arr = np.array(vecs, dtype="float32")
        arr /= (np.linalg.norm(arr, axis=1, keepdims=True) + 1e-12)
        return arr

    def encode(self, texts: List[str], batch_size: int = 64) -> np.ndarray:
        if _HAVE_LMSDK:
            return self._encode_sdk(texts, batch_size=batch_size)
        return self._encode_rest(texts, batch_size=batch_size)

class LMChat:
    def __init__(self, base: str = LMSTUDIO_BASE, model: str = LMSTUDIO_CHAT_MODEL):
        self.base = base.rstrip("/")
        self.model = model
    def chat(self, system: str, user: str, temperature: float = 0.2, max_tokens: int = 384) -> str:
        body = {
            "model": self.model,
            "messages": [{"role":"system","content":system},{"role":"user","content":user}],
            "temperature": temperature,
            "max_tokens": max_tokens,
            "stream": False
        }
        r = requests.post(f"{self.base}/v1/chat/completions", headers=HEADERS_JSON, json=body, timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        return r.json()["choices"][0]["message"]["content"]