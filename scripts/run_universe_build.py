# scripts/run_universe_build.py
from __future__ import annotations
import argparse, pathlib, sys, json
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
if str(SRC)  not in sys.path: sys.path.insert(0, str(SRC))

from pipeline.universe import build_universe
from utils.io import jdump

def main(args):
    queries = [q.strip() for q in args.queries.split("||") if q.strip()]
    uni = build_universe(
        queries=queries,
        year_min=args.year_min,
        year_max=args.year_max,
        retmax=args.retmax,
        hydrate=args.hydrate,
        hops=args.hops,
        per_seed_budget=args.per_seed_budget,
        knn_k=args.knn_k,
        alpha=args.alpha,
        beta=args.beta,
        resolution=args.resolution,
        threshold=args.threshold,
        emb_batch=args.emb_batch
    )
    outdir = pathlib.Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "universe.json"
    jdump(uni, path)

    df = pd.DataFrame(uni["docs"])
    print(f"✔ wrote {path}  | docs={uni['count']} themes={len(uni['themes'])} method={uni['cluster_method']}")
    yrs = pd.to_numeric(df["year"], errors="coerce")
    if yrs.notna().any():
        print(f"Years: {int(yrs.min())}–{int(yrs.max())} (med {float(yrs.median()):.1f})")
    sizes = [(t["theme_id"], t["size"]) for t in uni["themes"]]
    sizes.sort(key=lambda x: x[1], reverse=True)
    print("Theme sizes:", ", ".join([f"{tid}:{sz}" for tid,sz in sizes[:8]]))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--queries", required=True, help='Pipe multiple queries with "||"')
    ap.add_argument("--year-min", dest="year_min", type=int, default=None)
    ap.add_argument("--year-max", dest="year_max", type=int, default=None)
    ap.add_argument("--retmax", type=int, default=500)
    ap.add_argument("--hydrate", choices=["none","refs","citers","both"], default="none")
    ap.add_argument("--hops", type=int, default=1)
    ap.add_argument("--per-seed-budget", dest="per_seed_budget", type=int, default=150)
    ap.add_argument("--knn-k", type=int, default=20)
    ap.add_argument("--alpha", type=float, default=0.6)
    ap.add_argument("--beta", type=float, default=0.4)
    ap.add_argument("--resolution", type=float, default=0.6)
    ap.add_argument("--threshold", type=float, default=0.4)
    ap.add_argument("--emb-batch", dest="emb_batch", type=int, default=48)
    ap.add_argument("--outdir", default="runs/universe")
    args = ap.parse_args()
    main(args)
