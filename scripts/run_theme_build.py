# scripts/run_theme_build.py
from __future__ import annotations

import argparse, pathlib, sys, re, time
from typing import List, Dict, Any
import numpy as np, pandas as pd, networkx as nx
from tqdm import tqdm
from sklearn.feature_extraction.text import ENGLISH_STOP_WORDS


# --- Make both project root and src/ importable ---
ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

# --- Local imports ---
from config import KNN_K, HYBRID_ALPHA, HYBRID_BETA, LMSTUDIO_EMB_MODEL
from cache.emb import EmbCache
from utils.io import jdump
from clients.entrez import esearch, efetch_abstracts
from clients.icite import get_pubs, extract_refs_and_citers
from clients.lmstudio import LMEmbeddings
from themes.hybrid_graph import cosine_sim_matrix, build_knn, hybrid_weights
from themes.themes import soft_membership

TOKEN_RE = re.compile(r"[A-Za-zÀ-ÿ0-9_]+")

PRIM_HINT = {"Randomized Controlled Trial","Clinical Trial","Cohort","Case-Control"}
SR_HINT   = {"Systematic Review","Meta-Analysis","Review"}

def summarize_pubtypes(pt_lists):
    prim=sr=other=0
    for pts in pt_lists:
        s=set(pts or [])
        if s & PRIM_HINT: prim+=1
        elif s & SR_HINT: sr+=1
        else: other+=1
    return prim, sr, other

STOP = ENGLISH_STOP_WORDS

def top_keywords(texts, topn: int = 8):
    tf = {}
    for t in texts:
        for tok in TOKEN_RE.findall((t or "").lower()):
            if len(tok) < 4 or tok in STOP:
                continue
            tf[tok] = tf.get(tok, 0) + 1
    return [w for w,_ in sorted(tf.items(), key=lambda x: x[1], reverse=True)[:topn]]

def nearest_to_centroid(vecs: np.ndarray, idxs: List[int], centroid: np.ndarray, k: int = 3) -> List[int]:
    # return indices (within idxs) of the k closest docs to centroid
    sub = vecs[idxs]
    sims = (sub @ centroid) / (np.linalg.norm(sub, axis=1) + 1e-12)
    order = np.argsort(-sims)[:k]
    return [idxs[i] for i in order]

def build(args):
    outdir = pathlib.Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    # 1) PubMed search → PMIDs
    pmids = esearch(args.query, retmax=args.retmax, mindate=args.year_min, maxdate=args.year_max, sort="date")
    print(f"PubMed IDs: {len(pmids)}")
    if not pmids:
        print("Nothing found for this query/time window.")
        return

    # 2) Fetch metadata (title, abstract, year, pub types, doi)
    meta = efetch_abstracts(pmids)
    if not meta:
        print("efetch returned no metadata.")
        return
    df = pd.DataFrame.from_records(list(meta.values()))
    df = df.dropna(subset=["title"]).reset_index(drop=True)
    n = len(df)
    print(f"Fetched metadata for {n} items")

    # 3) Embeddings (cached)
    cache = EmbCache()
    pmid_list = [str(x) for x in df["pmid"].tolist()]
    texts = (df["title"].fillna("") + "\n" + df["abstract"].fillna("")).tolist()

    tE = time.perf_counter()
    cached = cache.get_many(LMSTUDIO_EMB_MODEL, pmid_list)
    hits = len(cached); miss = len(pmid_list) - hits
    print(f"Embedding cache: hits={hits} miss={miss}")

    vecs = np.zeros((len(pmid_list), 0), dtype="float32")  # placeholder for shape
    miss_idx, miss_texts, miss_pmids = [], [], []
    for i, pmid in enumerate(pmid_list):
        if pmid in cached:
            continue
        miss_idx.append(i); miss_texts.append(texts[i]); miss_pmids.append(pmid)

    if miss_texts:
        emb = LMEmbeddings()
        new_vecs = emb.encode(miss_texts, batch_size=args.emb_batch)
        # persist
        cache.put_many(LMSTUDIO_EMB_MODEL, [(pmid, new_vecs[j]) for j, pmid in enumerate(miss_pmids)])
    else:
        new_vecs = np.empty((0, 0), dtype="float32")

    # build full matrix in original order
    # first, obtain dim from either cached or new_vecs
    if hits > 0:
        any_vec = next(iter(cached.values()))
        dim = any_vec.size
    elif miss > 0:
        dim = new_vecs.shape[1]
    else:
        raise RuntimeError("No documents to embed.")
    vecs = np.zeros((len(pmid_list), dim), dtype="float32")

    # fill from cache
    for i, pmid in enumerate(pmid_list):
        if pmid in cached:
            vecs[i] = cached[pmid]
    # fill newly computed
    for pos, j in enumerate(miss_idx):
        vecs[j] = new_vecs[pos]

    # L2-normalize (safety)
    vecs /= (np.linalg.norm(vecs, axis=1, keepdims=True) + 1e-12)
    tE = time.perf_counter() - tE
    n_dim = int(vecs.shape[1]) if vecs.ndim == 2 else -1
    print(f"Embeddings: shape={vecs.shape}, batch={args.emb_batch}, time={tE:.2f}s")


    # 4) iCite references/citers for coupling
    pubs = get_pubs(df["pmid"].tolist(), fields=["pmid","cited_by","references","doi","year"], legacy=True)
    ref_sets: dict[str, set[int]] = {}
    for rec in pubs:
        refs, citers = extract_refs_and_citers(rec)
        pmid = str(rec.get("pmid") or rec.get("_id") or "")
        ref_sets[pmid] = set(refs)

    # 5) Cosine similarity (semantic)
    cos = cosine_sim_matrix(vecs)

    # 6) Bibliographic coupling (Jaccard) on kNN pairs
    knn_idx_cos, knn_cos = build_knn(cos, k=args.knn_k or KNN_K)
    bc_knn = np.zeros_like(knn_cos, dtype="float32")
    tC = time.perf_counter()
    for i in range(n):
        Ri = ref_sets.get(str(df.iloc[i]["pmid"]), set())
        if not Ri:
            continue
        for t, j in enumerate(knn_idx_cos[i]):
            Rj = ref_sets.get(str(df.iloc[j]["pmid"]), set())
            if not Rj:
                continue
            inter = len(Ri & Rj)
            if inter:
                uni = len(Ri | Rj) or 1
                bc_knn[i, t] = inter / uni
    tC = time.perf_counter() - tC
    print(f"Coupling (kNN pairs): n={n}, k={knn_idx_cos.shape[1]}, time={tC:.2f}s")

    # 7) Hybrid weights for the kNN graph
    hyb = hybrid_weights(knn_cos, bc_knn, alpha=args.alpha or HYBRID_ALPHA, beta=args.beta or HYBRID_BETA)

    # 8) Build weighted kNN graph
    G = nx.Graph()
    for i in range(n):
        for t, j in enumerate(knn_idx_cos[i]):
            w = float(hyb[i, t])
            if w <= 0:
                continue
            if G.has_edge(i, j):
                if G[i][j]["weight"] < w:
                    G[i][j]["weight"] = w
            else:
                G.add_edge(i, j, weight=w)

    # 9) Cluster: Leiden → HDBSCAN → thresholded components
    labels = None
    method = None
    try:
        import igraph as ig  # type: ignore
        import leidenalg as la  # type: ignore
        mapping = {n:i for i,n in enumerate(G.nodes())}
        edges = [(mapping[u], mapping[v]) for u,v in G.edges()]
        weights = [G[u][v].get("weight", 1.0) for u,v in G.edges()]
        g = ig.Graph(n=len(mapping), edges=edges)
        g.es["weight"] = weights
        part = la.find_partition(g, la.RBConfigurationVertexPartition, weights="weight", resolution_parameter=args.resolution)
        labels_ig = np.zeros(len(mapping), dtype=int)
        for comm_id, members in enumerate(part):
            for v in members:
                labels_ig[v] = comm_id
        labels = labels_ig
        method = "leiden"
    except Exception:
        try:
            import importlib.util
            if importlib.util.find_spec("hdbscan") is None:
                raise ImportError("hdbscan not installed")
            import hdbscan  # type: ignore
            labels = hdbscan.HDBSCAN(min_cluster_size=max(10, n // 50), metric="euclidean").fit_predict(vecs)
            method = "hdbscan"
        except Exception:
            TH = float(args.threshold or 0.4)
            H = nx.Graph((u,v,d) for u,v,d in G.edges(data=True) if d.get("weight",0.0) >= TH)
            labels = -1 * np.ones(n, dtype=int)
            cid = 0
            for comp in nx.connected_components(H):
                for i in comp:
                    labels[i] = cid
                cid += 1
            method = f"components@{TH}"

    # 10) Soft membership (top-2 themes per doc)
    unique, W = soft_membership(vecs, labels, knn_idx_cos, hyb, topm=2, lam=0.5)

    # 11) Persist artifacts
    themes = []
    for t in unique:
        members = np.where(labels == t)[0].tolist()
        if not members:
            continue
        cent = vecs[members].mean(axis=0)
        cent = cent / (np.linalg.norm(cent) + 1e-12)
        yrs = pd.to_numeric(df.iloc[members]["year"], errors="coerce")
        theme = {
            "theme_id": int(t),
            "size": len(members),
            "members_idx": members,
            "members_pmids": [str(df.iloc[i]["pmid"]) for i in members],
            "centroid": cent.tolist(),
            "year_stats": {
                "min": int(yrs.min()) if yrs.notna().any() else None,
                "max": int(yrs.max()) if yrs.notna().any() else None,
                "median": float(yrs.median()) if yrs.notna().any() else None,
            },
        }
        themes.append(theme)

    out = {
        "query": args.query,
        "count": n,
        "cluster_method": method,
        "themes": themes,
        "docs": df.to_dict(orient="records"),
    }
    jdump(out, outdir / "themes.json")

    # 12) HUMAN-FRIENDLY TERMINAL SUMMARY
    print("\n=== RUN SUMMARY ===")
    print(f"Query: {args.query}")
    print(f"Window: {args.year_min}–{args.year_max}  | retmax={args.retmax}")
    print(f"Docs: {n} | Embedding dim: {n_dim} | kNN-k={args.knn_k} | alpha={args.alpha} beta={args.beta}")
    print(f"Cluster method: {method} | Themes: {len(themes)}")
    # quick global stats
    yrs_all = pd.to_numeric(df["year"], errors="coerce")
    if yrs_all.notna().any():
        print(f"Year range (all docs): {int(yrs_all.min())}–{int(yrs_all.max())}  median={float(yrs_all.median()):.1f}")

    # theme-wise preview with pubtypes & journals
    for theme in themes:
        tid = theme["theme_id"]; members = theme["members_idx"]
        y = theme["year_stats"]

        # keywords
        t_texts = df.iloc[members]["title"].fillna("").tolist()
        kws = top_keywords(t_texts, topn=8)

        # representatives
        cent = np.array(theme["centroid"], dtype="float32")
        reps = nearest_to_centroid(vecs, members, cent, k=3)
        rep_titles = [f"- {df.iloc[i]['title'][:140].strip()}" for i in reps]

        # pub types summary
        prim, sr, other = summarize_pubtypes(df.iloc[members]["pub_types"].tolist())

        # top journals
        top_j = (df.iloc[members]["journal"].fillna("")
                .replace("", np.nan).dropna()
                .value_counts().head(3).to_dict())

        print(f"\nTheme {tid}  | size={theme['size']} | years {y['min']}–{y['max']} (med {y['median']})")
        print("  keywords:", ", ".join(kws) if kws else "(n/a)")
        print(f"  pub types: primary={prim}  SR/Review={sr}  other={other}")
        if top_j:
            tj = ", ".join([f"{k}×{v}" for k,v in top_j.items()])
            print(f"  top journals: {tj}")
        print("  reps:")
        for rt in rep_titles:
            print("   ", rt)


    dt = time.perf_counter() - t0
    print(f"\n✔ wrote {outdir/'themes.json'}  | total time {dt:.1f}s")
    print("===============")

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True, help="PubMed query")
    ap.add_argument("--year-min", dest="year_min", type=int, default=None)
    ap.add_argument("--year-max", dest="year_max", type=int, default=None)
    ap.add_argument("--retmax", type=int, default=600)
    ap.add_argument("--knn-k", type=int, default=20)
    ap.add_argument("--alpha", type=float, default=0.6)
    ap.add_argument("--beta", type=float, default=0.4)
    ap.add_argument("--resolution", type=float, default=0.6)
    ap.add_argument("--threshold", type=float, default=0.4, help="edge weight threshold if components fallback is used")
    ap.add_argument("--emb-batch", dest="emb_batch", type=int, default=48, help="embedding batch size (lower to reduce VRAM)")
    ap.add_argument("--outdir", default="runs/theme_build")
    args = ap.parse_args()
    build(args)
