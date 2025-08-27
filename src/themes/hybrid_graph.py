from __future__ import annotations
import numpy as np
from typing import Tuple

def cosine_sim_matrix(X: np.ndarray, Y: np.ndarray | None = None) -> np.ndarray:
    if Y is None: Y = X
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
    Yn = Y / (np.linalg.norm(Y, axis=1, keepdims=True) + 1e-12)
    return Xn @ Yn.T

def build_knn(sims: np.ndarray, k: int = 20) -> Tuple[np.ndarray, np.ndarray]:
    n = sims.shape[0]
    k = min(k, n-1) if n > 1 else 1
    idx = np.argpartition(-sims, kth=k-1, axis=1)[:, :k]
    vals = np.take_along_axis(sims, idx, axis=1)
    return idx, vals

def hybrid_weights(cosine_knn_sims: np.ndarray, coupling_knn_sims: np.ndarray | None, alpha: float, beta: float) -> np.ndarray:
    if coupling_knn_sims is None:
        return cosine_knn_sims
    return alpha * cosine_knn_sims + beta * coupling_knn_sims
