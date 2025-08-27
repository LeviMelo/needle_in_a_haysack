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
    cover_rows = [coverage_for_theme(t, docs_df) for t in uni["themes"]]
    now_year = datetime.now(timezone.utc).year
    ranked = rank_gaps(uni, cover_rows, now_year)

    # Optional LLM labeling for the top themes
    llm = LMChat() if args.llm_label else None
    if llm:
        id2members = {t["theme_id"]: t["members_idx"] for t in uni["themes"]}
        for r in ranked[:args.topk]:
            idxs = id2members[r["theme_id"]]
            titles = [docs_df.iloc[i]["title"] for i in idxs][:30]
            r.update(maybe_llm_label(llm, titles))

    outdir = pathlib.Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    jdump({"universe_file": str(args.universe), "coverage": cover_rows, "ranked": ranked}, outdir / "gaps.json")

    print(f"✔ wrote {outdir/'gaps.json'}")
    print("\n=== TOP CANDIDATES ===")
    for r in ranked[:args.topk]:
        print(f"- Theme {r['theme_id']}: GAP={r['gap_score']:.3f} | cov={r['coverage_ratio']:.2f} ({r['coverage_level']}) | E={r['E_size']} | new={r['new_primary_count']} | lastSR={r['last_sr_year']}")
        if args.llm_label and "llm_yaml" in r:
            print("  LLM label/questions →")
            print("  " + r["llm_yaml"].replace("\n", "\n  "))
        else:
            if r["questions"]:
                print("  Q:", r["questions"][0])
            if r["terms"]:
                print("  terms:", ", ".join(r["terms"][:8]))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", required=True)
    ap.add_argument("--outdir", default="runs/gap_hunt")
    ap.add_argument("--topk", type=int, default=6)
    ap.add_argument("--llm-label", action="store_true", help="use local LLM (LM Studio) to label themes and draft questions")
    args = ap.parse_args()
    main(args)
