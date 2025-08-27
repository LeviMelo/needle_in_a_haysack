from __future__ import annotations
import os

LMSTUDIO_BASE = os.getenv("LMSTUDIO_BASE", "http://127.0.0.1:1234")
LMSTUDIO_EMB_MODEL = os.getenv("LMSTUDIO_EMB_MODEL", "text-embedding-qwen3-embedding-0.6b")
LMSTUDIO_CHAT_MODEL = os.getenv("LMSTUDIO_CHAT_MODEL", "gemma-3n-e2b-it")

ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL", "you@example.com")
ENTREZ_API_KEY = os.getenv("ENTREZ_API_KEY", "")

ICITE_BASE = os.getenv("ICITE_BASE", "https://icite.od.nih.gov/api")
HTTP_TIMEOUT = int(os.getenv("HTTP_TIMEOUT", "30"))
USER_AGENT = os.getenv("USER_AGENT", "litgap-poc/0.1 (+https://example.org)")

KNN_K = int(os.getenv("KNN_K", "20"))
HYBRID_ALPHA = float(os.getenv("HYBRID_ALPHA", "0.6"))
HYBRID_BETA  = float(os.getenv("HYBRID_BETA", "0.4"))


# Coverage thresholds
COV_LEVELS = {
    "NONE": 0.2,
    "LOW": 0.5,
    "SUBSTANTIAL": 0.8,
    "NEAR_FULL": 0.95,
    "FULL": 1.01,   # anything >=0.95 treated as full
}
