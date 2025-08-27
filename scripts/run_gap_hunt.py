# scripts/run_gap_hunt.py
from __future__ import annotations
import argparse, pathlib, sys, json
from typing import Dict, Any, List
import pandas as pd
from datetime import datetime, timezone

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
if str(SRC)  not in sys.path: sys.path.insert(0, str(SRC))

from utils.io import jdump
from pipeline.coverage import coverage_for_theme
from pipeline.gap import rank_gaps, top_terms
from clients.lmstudio import LMChat

LABEL_SYS = "You summarize biomedical literature themes. Be concise and precise."
LABEL_TMPL = """Given these paper titles (newline separated), return:
1) A short theme label (≤7 words, no punctuation noise)
2) 2–3 candidate research questions (1 sentence each), neutral and specific.

Titles:
{titles}

Return YAML with keys: label, questions"""

def maybe_llm_label(chat: LMChat, titles: List[str]) -> Dict[str, Any]:
    try:
        resp = chat.chat(LABEL_SYS, LABEL_TMPL.format(titles="\n".join(titles)), temperature=0.2, max_tokens=400)
        return {"llm_yaml": resp}
    except Exception as e:
        return {"llm_error": str(e)}

def main(args):
    uni = json.loads(pathlib.Path(args.universe).read_text(encoding="utf-8"))
    docs_df = pd.DataFrame(uni["docs"])

    # coverage with details for auditing
    cover_rows = [coverage_for_theme(t, docs_df, include_docs=True) for t in uni["themes"]]
    now_year = datetime.now(timezone.utc).year
    ranked = rank_gaps(uni, cover_rows, now_year)

    # Optional LLM labeling for the top themes
    llm = LMChat() if args.llm_label else None
    id2members = {t["theme_id"]: t["members_idx"] for t in uni["themes"]}
    if llm:
        for r in ranked[:args.topk]:
            idxs = id2members[r["theme_id"]]
            titles = [docs_df.iloc[i]["title"] for i in idxs][:30]
            r.update(maybe_llm_label(llm, titles))

    outdir = pathlib.Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    payload = {"universe_file": str(args.universe), "coverage": cover_rows, "ranked": ranked}
    jdump(payload, outdir / "gaps.json")

    # --- Terminal report (auditable) ---
    print(f"✔ wrote {outdir/'gaps.json'}")

    print("\n=== TOP CANDIDATES ===")
    for r in ranked[:args.topk]:
        tid = r["theme_id"]
        cov = r["coverage_ratio"]; lvl = r["coverage_level"]
        E  = r["E"]; S = r["S"]
        last = r["last_sr_year"]
        print(f"- Theme {tid}: GAP={r['gap_score']:.3f} | cov={cov:.2f} ({lvl}) | E={len(E)} | new={r['new_primary_count']} | lastSR={last}")

        if args.llm_label and "llm_yaml" in r:
            print("  LLM label/questions →")
            print("  " + r["llm_yaml"].replace("\n", "\n  "))
        else:
            if r["questions"]:
                print("  Q:", r["questions"][0])

        # semantic + biblio snippets (if present in coverage)
        cov_row = next(cr for cr in cover_rows if cr["theme_id"] == tid)
        det = cov_row.get("details", {})
        semE = det.get("sem", {}).get("E", {})
        semS = det.get("sem", {}).get("S", {})
        bibE = det.get("bib", {}).get("E_jaccard_refs", {})
        if semE.get("n"):
            print(f"  sem(E): mean={semE['mean']:.3f} median={semE['median']:.3f} n={semE['n']}")
        if semS.get("n"):
            print(f"  sem(S): mean={semS['mean']:.3f} median={semS['median']:.3f} n={semS['n']}")
        if bibE.get("pairs"):
            print(f"  bib(E Jaccard refs): mean={bibE['mean']:.3f} p90={bibE['p90']:.3f} pairs={bibE['pairs']}")

        # show a few primaries and SRs with titles
        E_docs = det.get("E_docs", {})
        S_docs = det.get("S_docs", {})
        if E_docs:
            print("  Primaries (top 5):")
            for pid in list(E_docs.keys())[:5]:
                d = E_docs[pid]
                print(f"    - {pid} [{d.get('year')}] {d.get('title')[:120]}")
        if S_docs:
            print("  SR/Reviews (top 5):")
            for pid in list(S_docs.keys())[:5]:
                d = S_docs[pid]
                print(f"    - {pid} [{d.get('year')}] {d.get('title')[:120]}")

        # per-SR coverage breakdown (first 3)
        br = det.get("sr_breakdown", [])[:3]
        if br:
            print("  SR breakdown:")
            for row in br:
                print(f"    - {row['pmid']} [{row.get('year')}] inc|E={len(row['included_primaries'])}/{len(E)} ({row['coverage_of_E']:.2f})")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", required=True)
    ap.add_argument("--outdir", default="runs/gap_hunt")
    ap.add_argument("--topk", type=int, default=6)
    ap.add_argument("--llm-label", action="store_true", help="use local LLM (LM Studio) to label themes and draft questions")
    args = ap.parse_args()
    main(args)
