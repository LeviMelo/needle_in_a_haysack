# scripts/run_ripple_boost.py
from __future__ import annotations
import argparse, pathlib, sys, json
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
if str(SRC)  not in sys.path: sys.path.insert(0, str(SRC))

from pipeline.ripple import ripple_expand_from_primaries
from pipeline.evidence import split_by_kind
from utils.io import jdump

def main(args):
    uni = json.loads(pathlib.Path(args.universe).read_text(encoding="utf-8"))
    docs_df = pd.DataFrame(uni["docs"]).copy()
    t = next((x for x in uni["themes"] if int(x["theme_id"])==int(args.theme_id)), None)
    if not t:
        print(f"Theme {args.theme_id} not found.")
        return
    members = set(t["members_pmids"])
    sub = docs_df[docs_df["pmid"].astype(str).isin(members)].copy()
    prim, sr, _ = split_by_kind(sub.to_dict(orient="records"))
    if not prim:
        print("No primaries in theme; nothing to ripple from.")
        return

    res = ripple_expand_from_primaries(prim, allowed_since_year=args.since_year,
                                       max_expand=args.max_expand, prefer=args.prefer)

    outdir = pathlib.Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"ripple_theme{args.theme_id}.json"
    jdump(res, outpath)

    # terminal stats
    st = res.get("stats", {})
    print(f"Seeds: {len(prim)} | prefer={st.get('prefer')} | since>={st.get('allowed_since_year')}")
    if st:
        print(f"Pool unique: {st.get('pool_unique')} | after year gate: {st.get('after_year_gate')} | after cap: {st.get('after_cap')}")
        if st.get("reason_zero"):
            print("Reason zero candidates:", st["reason_zero"])
        # per-seed preview (first 5)
        for row in st.get("seeds", [])[:5]:
            print(f"  seed {row['seed']}: refs={row['refs']} citers={row['citers']} used={row['used']} pool={row['pool_size']}")

    print(f"✔ wrote {outpath} | candidates={len(res['candidates'])}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", required=True)
    ap.add_argument("--theme-id", type=int, required=True)
    ap.add_argument("--since-year", type=int, default=None, help="only include items >= year")
    ap.add_argument("--max-expand", type=int, default=300)
    ap.add_argument("--prefer", choices=["citers","refs","both"], default="citers")  # keep CLI superset (pipeline ignores 'both')
    ap.add_argument("--outdir", default="runs/ripple")
    args = ap.parse_args()
    main(args)
