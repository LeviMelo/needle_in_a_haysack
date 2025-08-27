# src/pipeline/universe.py
from __future__ import annotations
import time, pathlib, sys
from typing import List, Dict, Any, Literal, Set
import numpy as np, pandas as pd, networkx as nx

# import roots
ROOT = pathlib.Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
if str(SRC)  not in sys.path: sys.path.insert(0, str(SRC))

from clients.entrez import esearch, efetch_abstracts
from clients.icite import get_pubs, extract_refs_and_citers
from cache.icite import ICiteCache
from clients.lmstudio import LMEmbeddings
from themes.hybrid_graph import cosine_sim_matrix, build_knn, hybrid_weights
from themes.themes import soft_membership
from config import KNN_K, HYBRID_ALPHA, HYBRID_BETA, LMSTUDIO_EMB_MODEL
from cache.emb import EmbCache

def hydrate_pmids(seed_pmids: List[str],
                  mode: Literal["none","refs","citers","both"]="none",
                  hops: int = 1,
                  per_seed_budget: int = 200) -> List[str]:
    """
    Optional node expansion using iCite refs/citers with hard budgets.
    Returns the union list (seed + expansion), de-duplicated, as strings.
    """
    if mode == "none" or hops <= 0:
        return [str(x) for x in seed_pmids]
    cache = ICiteCache()
    frontier: Set[str] = set(str(x) for x in seed_pmids)
    universe: Set[str] = set(frontier)
    for _ in range(hops):
        need = list(frontier)
        have = cache.get_many(need, legacy=True)
        missing = [p for p in need if p not in have]
        if missing:
            fetched = get_pubs(missing, fields=["pmid","references","cited_by"], legacy=True)
            cache.put_many(fetched, legacy=True)
            for rec in fetched:
                have[str(rec.get("pmid") or rec.get("_id") or "")] = rec
        new: Set[str] = set()
        for p in need:
            rec = have.get(p, {})
            refs, citers = extract_refs_and_citers(rec)
            pool: List[int] = []
            if mode in ("refs","both"):
                pool.extend(refs[:per_seed_budget//(2 if mode=="both" else 1)])
            if mode in ("citers","both"):
                pool.extend(citers[:per_seed_budget//(2 if mode=="both" else 1)])
            for q in pool:
                q = str(q)
                if q not in universe:
                    new.add(q)
        frontier = new
        universe |= new
        if not frontier: break
    return list(universe)

def build_universe(queries: List[str],
                   year_min: int | None,
                   year_max: int | None,
                   retmax: int = 500,
                   hydrate: Literal["none","refs","citers","both"]="none",
                   hops: int = 1,
                   per_seed_budget: int = 150,
                   knn_k: int = KNN_K,
                   alpha: float = HYBRID_ALPHA,
                   beta: float  = HYBRID_BETA,
                   resolution: float = 0.8,
                   threshold: float = 0.4,
                   emb_batch: int = 48) -> Dict[str,Any]:
    """
    Multi-query -> (optional) hydration -> efetch -> cached embeddings -> hybrid kNN -> clustering.
    Returns full universe dict {docs, themes, cluster_method, params...}.
    """
    t0 = time.perf_counter()
    # 1) union PMIDs from all queries
    pmid_set: Set[str] = set()
    for q in queries:
        ids = esearch(q, retmax=retmax, mindate=year_min, maxdate=year_max, sort="date")
        pmid_set.update(str(x) for x in ids)
    pmids = list(pmid_set)
    # 2) hydrate optionally
    pmids_h = hydrate_pmids(pmids, mode=hydrate, hops=hops, per_seed_budget=per_seed_budget)
    # 3) efetch metadata
    meta = efetch_abstracts(pmids_h)
    df = pd.DataFrame.from_records(list(meta.values()))
    df = df.dropna(subset=["title"]).reset_index(drop=True)
    # 4) embeddings (cached)
    cache = EmbCache()
    pid = [str(x) for x in df["pmid"].tolist()]
    texts = (df["title"].fillna("") + "\n" + df["abstract"].fillna("")).tolist()
    cached = cache.get_many(LMSTUDIO_EMB_MODEL, pid)
    miss_idx, miss_texts, miss_pmids = [], [], []
    for i, p in enumerate(pid):
        if p not in cached:
            miss_idx.append(i); miss_texts.append(texts[i]); miss_pmids.append(p)
    if miss_texts:
        emb = LMEmbeddings()
        new_vecs = emb.encode(miss_texts, batch_size=emb_batch)
        cache.put_many(LMSTUDIO_EMB_MODEL, [(miss_pmids[j], new_vecs[j]) for j in range(len(miss_pmids))])
    else:
        new_vecs = np.empty((0,0), dtype="float32")
    # assemble vecs
    if cached:
        dim = next(iter(cached.values())).size
    else:
        dim = new_vecs.shape[1]
    vecs = np.zeros((len(pid), dim), dtype="float32")
    for i,p in enumerate(pid):
        if p in cached:
            vecs[i] = cached[p]
    for pos, j in enumerate(miss_idx):
        vecs[j] = new_vecs[pos]
    vecs /= (np.linalg.norm(vecs, axis=1, keepdims=True) + 1e-12)
    # 5) hybrid kNN graph
    cos = cosine_sim_matrix(vecs)
    knn_idx, knn_cos = build_knn(cos, k=knn_k)
    # bibliographic coupling on kNN pairs
    icache = ICiteCache()
    have = icache.get_many(pid, legacy=True)
    need = [p for p in pid if p not in have]
    if need:
        fetched = get_pubs(need, fields=["pmid","references","cited_by"], legacy=True)
        icache.put_many(fetched, legacy=True)
        for rec in fetched:
            have[str(rec.get("pmid") or rec.get("_id") or "")] = rec
    ref_sets: Dict[str,set[int]] = {}
    for p in pid:
        refs,_ = extract_refs_and_citers(have.get(p, {}))
        ref_sets[p] = set(refs)
    bc_knn = np.zeros_like(knn_cos, dtype="float32")
    for i in range(len(pid)):
        Ri = ref_sets.get(pid[i], set())
        if not Ri: continue
        for t, j in enumerate(knn_idx[i]):
            Rj = ref_sets.get(pid[j], set())
            if not Rj: continue
            inter = len(Ri & Rj)
            if inter:
                uni = len(Ri | Rj) or 1
                bc_knn[i,t] = inter/uni
    hyb = hybrid_weights(knn_cos, bc_knn, alpha=alpha, beta=beta)
    # 6) cluster
    G = nx.Graph()
    for i in range(len(pid)):
        for t, j in enumerate(knn_idx[i]):
            w = float(hyb[i,t])
            if w<=0: continue
            if G.has_edge(i,j):
                if G[i][j]["weight"] < w:
                    G[i][j]["weight"] = w
            else:
                G.add_edge(i,j,weight=w)
    labels=None; method=None
    try:
        import igraph as ig, leidenalg as la
        mapping = {n:i for i,n in enumerate(G.nodes())}
        edges = [(mapping[u],mapping[v]) for u,v in G.edges()]
        weights = [G[u][v].get("weight",1.0) for u,v in G.edges()]
        g = ig.Graph(n=len(mapping), edges=edges)
        g.es["weight"] = weights
        part = la.find_partition(g, la.RBConfigurationVertexPartition, weights="weight", resolution_parameter=resolution)
        labels = np.zeros(len(mapping), dtype=int)
        for cid, members in enumerate(part):
            for v in members:
                labels[v] = cid
        method="leiden"
    except Exception:
        try:
            import importlib.util, hdbscan  # type: ignore
            labels = hdbscan.HDBSCAN(min_cluster_size=max(10, len(pid)//50), metric="euclidean").fit_predict(vecs)
            method="hdbscan"
        except Exception:
            TH = float(threshold)
            H = nx.Graph((u,v,d) for u,v,d in G.edges(data=True) if d.get("weight",0.0)>=TH)
            labels = -1*np.ones(len(pid), dtype=int); cid=0
            for comp in nx.connected_components(H):
                for i in comp: labels[i]=cid
                cid+=1
            method=f"components@{TH}"
    uniq, W = soft_membership(vecs, labels, knn_idx, hyb, topm=2, lam=0.5)
    # 7) package
    themes=[]
    for t in uniq:
        members = np.where(labels==t)[0].tolist()
        if not members: continue
        cent = vecs[members].mean(axis=0)
        cent = cent/(np.linalg.norm(cent)+1e-12)
        yrs = pd.to_numeric(pd.Series([meta[str(df.iloc[i]['pmid'])]['year'] if str(df.iloc[i]['pmid']) in meta else df.iloc[i]['year'] for i in members]), errors="coerce")
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
        "queries": queries,
        "params": {"year_min":year_min,"year_max":year_max,"retmax":retmax,"hydrate":hydrate,"hops":hops,
                   "per_seed_budget":per_seed_budget,"knn_k":knn_k,"alpha":alpha,"beta":beta,"resolution":resolution},
        "count": len(pid),
        "cluster_method": method,
        "themes": themes,
        "docs": df.to_dict(orient="records"),
    }
    return out
