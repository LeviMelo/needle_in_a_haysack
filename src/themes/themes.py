from __future__ import annotations
import numpy as np

def soft_membership(
    doc_vecs: np.ndarray,
    labels: np.ndarray,
    kNN_idx: np.ndarray,
    kNN_sims: np.ndarray,
    topm: int = 2,
    lam: float = 0.5
):
    # Centroids
    uniq = sorted(set(int(x) for x in labels if x >= 0))
    centroids = {}
    for t in uniq:
        M = doc_vecs[labels==t]
        c = M.mean(axis=0)
        centroids[t] = c / (np.linalg.norm(c)+1e-12)

    # cosine to centroids
    cos_to_t = np.zeros((doc_vecs.shape[0], len(uniq)), dtype="float32")
    for j,t in enumerate(uniq):
        c = centroids[t]
        cos_to_t[:, j] = (doc_vecs @ c) / (np.linalg.norm(doc_vecs,axis=1)+1e-12)

    # neighbor similarity per theme (avg sim to neighbors with that label)
    avg_nei = np.zeros_like(cos_to_t)
    for i in range(doc_vecs.shape[0]):
        nei = kNN_idx[i]; sims = kNN_sims[i]
        for j,t in enumerate(uniq):
            mask = (labels[nei] == t)
            avg_nei[i,j] = sims[mask].mean() if mask.any() else 0.0

    raw = lam * cos_to_t + (1-lam) * avg_nei
    raw = raw - raw.max(axis=1, keepdims=True)
    e = np.exp(raw)
    W = e / (e.sum(axis=1, keepdims=True) + 1e-12)

    # keep only top-m
    top_idx = np.argsort(-W, axis=1)[:, :topm]
    W_sparse = np.zeros_like(W)
    rows = np.arange(W.shape[0])[:,None]
    W_sparse[rows, top_idx] = W[rows, top_idx]
    return uniq, W_sparse
