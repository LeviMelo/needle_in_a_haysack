#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pmid_embed_similarity.py

Fetch PubMed records (title + full abstract), embed them with a local LM Studio embeddings server,
and compute extensive similarity metrics.

Default: compares PMIDs 24548571 and 40016490 using:
  model='text-embedding-qwen3-embedding-0.6b@f16'
  lm-url='http://127.0.0.1:1234/v1/embeddings'

Usage examples
--------------
# 1) Quick run with your two PMIDs and defaults
python pmid_embed_similarity.py

# 2) Any set of PMIDs (space-separated)
python pmid_embed_similarity.py --pmids 24548571 40016490 33591115

# 3) Explicit model / URL (LM Studio)
python pmid_embed_similarity.py \
  --pmids 24548571 40016490 \
  --model text-embedding-qwen3-embedding-0.6b@f16 \
  --lm-url http://127.0.0.1:1234/v1/embeddings

# 4) Save artifacts (JSON and CSV) into an output folder
python pmid_embed_similarity.py --pmids 24548571 40016490 --outdir runs/sim_checks

# 5) Control chunking if texts are very long (mean-pooled embedding across chunks)
python pmid_embed_similarity.py --pmids 24548571 40016490 --max-chars 8000

# 6) Use only title or only abstract (or both; default is both)
python pmid_embed_similarity.py --pmids 24548571 40016490 --include-title --include-abstract
python pmid_embed_similarity.py --pmids 24548571 40016490 --include-title --no-include-abstract
python pmid_embed_similarity.py --pmids 24548571 40016490 --no-include-title --include-abstract

Notes
-----
- Cosine similarity is the main metric for embedding proximity.
- Euclidean/L1 distances are provided for additional texture but are scale-sensitive.
- If an abstract is structured into sections, we join all sections (label-prefixed).
- If there is no abstract, the script gracefully falls back to title-only (unless disabled).
- NCBI etiquette: you may add --ncbi-email to identify your tool.
"""
import argparse
import json
import math
import os
import re
import sys
import time
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Tuple, Optional

import numpy as np
import requests
import xml.etree.ElementTree as ET

EFETCH_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def clamp(x: float, lo: float = -1.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def angle_from_cosine(c: float) -> Tuple[float, float]:
    """Return angle in radians and degrees from cosine value (clamped)."""
    c = clamp(c, -1.0, 1.0)
    rad = math.acos(c)
    deg = rad * 180.0 / math.pi
    return rad, deg


def l2_norm(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


def l1_norm(v: np.ndarray) -> float:
    return float(np.sum(np.abs(v)))


def cosine_similarity(u: np.ndarray, v: np.ndarray) -> float:
    nu = np.linalg.norm(u)
    nv = np.linalg.norm(v)
    if nu == 0.0 or nv == 0.0:
        return float('nan')
    return float(np.dot(u, v) / (nu * nv))


def euclidean_distance(u: np.ndarray, v: np.ndarray) -> float:
    return float(np.linalg.norm(u - v))


def manhattan_distance(u: np.ndarray, v: np.ndarray) -> float:
    return float(np.sum(np.abs(u - v)))


def describe_vector(v: np.ndarray) -> Dict[str, float]:
    return {
        "dim": int(v.shape[0]),
        "min": float(np.min(v)),
        "max": float(np.max(v)),
        "mean": float(np.mean(v)),
        "std": float(np.std(v)),
        "l2_norm": l2_norm(v),
        "l1_norm": l1_norm(v),
    }


def fetch_pubmed_xml(pmid: str, email: Optional[str] = None, tool: str = "pmid_embedder") -> str:
    params = {
        "db": "pubmed",
        "id": pmid,
        "rettype": "xml",
        "retmode": "xml",
        "tool": tool,
    }
    if email:
        params["email"] = email
    r = requests.get(EFETCH_BASE, params=params, timeout=30)
    r.raise_for_status()
    return r.text


def get_text(elem: Optional[ET.Element]) -> str:
    return elem.text if elem is not None and elem.text else ""


def clean_whitespace(t: str) -> str:
    t = re.sub(r"\s+", " ", t or "").strip()
    return t


def parse_pubmed_record(xml_text: str) -> Dict[str, str]:
    """
    Return a dict with keys:
      pmid, title, abstract, year, journal, article_type (pipe-joined pub types)
    """
    root = ET.fromstring(xml_text)
    article = root.find(".//MedlineCitation/Article")
    pmid_el = root.find(".//MedlineCitation/PMID")
    pmid = get_text(pmid_el)

    title = clean_whitespace(get_text(article.find("./ArticleTitle")) if article is not None else "")

    # Abstract (may be structured with multiple AbstractText elements)
    abstract_full_parts = []
    if article is not None:
        for ab in article.findall("./Abstract/AbstractText"):
            label = ab.attrib.get("Label") or ab.attrib.get("NlmCategory")
            text = "".join(ab.itertext())  # AbstractText may contain nested tags
            text = clean_whitespace(text)
            if text:
                if label:
                    abstract_full_parts.append(f"{label}: {text}")
                else:
                    abstract_full_parts.append(text)
    abstract = clean_whitespace("\n\n".join(abstract_full_parts))

    # Journal & Year (best-effort)
    journal = ""
    year = ""
    if article is not None:
        journal = clean_whitespace(get_text(article.find("./Journal/Title")))
        year_el = article.find("./Journal/JournalIssue/PubDate/Year")
        medline_date = get_text(article.find("./Journal/JournalIssue/PubDate/MedlineDate"))
        if year_el is not None and year_el.text:
            year = clean_whitespace(year_el.text)
        elif medline_date:
            m = re.search(r"(\d{4})", medline_date)
            if m:
                year = m.group(1)

    # Publication types
    pub_types = []
    if article is not None:
        for pt in article.findall("./PublicationTypeList/PublicationType"):
            pt_text = clean_whitespace("".join(pt.itertext()))
            if pt_text:
                pub_types.append(pt_text)
    article_type = " | ".join(pub_types) if pub_types else ""

    return {
        "pmid": pmid,
        "title": title,
        "abstract": abstract,
        "year": year,
        "journal": journal,
        "article_type": article_type,
    }


def build_text_for_embedding(rec: Dict[str, str],
                             include_title: bool = True,
                             include_abstract: bool = True) -> str:
    parts = []
    if include_title and rec.get("title"):
        parts.append(rec["title"])
    if include_abstract and rec.get("abstract"):
        parts.append(rec["abstract"])
    return "\n\n".join(parts).strip()


def chunk_by_max_chars(text: str, max_chars: int) -> List[str]:
    if max_chars is None or max_chars <= 0 or len(text) <= max_chars:
        return [text]
    words = text.split()
    chunks = []
    cur = []
    cur_len = 0
    for w in words:
        wlen = len(w) + (1 if cur_len > 0 else 0)  # space if needed
        if cur_len + wlen <= max_chars:
            cur.append(w)
            cur_len += wlen
        else:
            if cur:
                chunks.append(" ".join(cur))
            cur = [w]
            cur_len = len(w)
    if cur:
        chunks.append(" ".join(cur))
    return chunks


def embed_texts_lmstudio(texts: List[str], model: str, lm_url: str) -> List[np.ndarray]:
    """
    Sends all texts in one request if supported by the server.
    Falls back to per-text requests if needed.
    Returns a list of np.ndarray embeddings in the same order.
    """
    headers = {"Content-Type": "application/json"}
    payload = {"model": model, "input": texts}
    t0 = time.time()
    r = requests.post(lm_url, headers=headers, data=json.dumps(payload), timeout=120)
    t1 = time.time()
    if r.status_code != 200:
        # Fallback: try per-text (some servers might not accept batch input)
        embs = []
        for t in texts:
            payload_one = {"model": model, "input": t}
            rr = requests.post(lm_url, headers=headers, data=json.dumps(payload_one), timeout=120)
            rr.raise_for_status()
            data = rr.json()
            emb = np.array(data["data"][0]["embedding"], dtype=np.float32)
            embs.append(emb)
        return embs

    data = r.json()
    # Expect OpenAI-like: {"data": [{"embedding": [...], "index": 0}, ...]}
    out = []
    for item in sorted(data.get("data", []), key=lambda x: x.get("index", 0)):
        out.append(np.array(item["embedding"], dtype=np.float32))
    # Timing can be useful
    payload["_elapsed_sec"] = round(t1 - t0, 3)
    return out


def embed_article_text(text: str, model: str, lm_url: str, max_chars: int) -> Tuple[np.ndarray, Dict]:
    """
    Chunk long text, embed each chunk, then mean-pool to a single vector.
    Returns (embedding, stats_dict).
    """
    chunks = chunk_by_max_chars(text, max_chars)
    chunk_embs = embed_texts_lmstudio(chunks, model=model, lm_url=lm_url)
    # Mean-pool (simple, robust; avoids overweighting)
    M = np.vstack(chunk_embs)
    pooled = np.mean(M, axis=0)
    stats = {
        "n_chunks": len(chunks),
        "chars_total": len(text),
        "chars_per_chunk": [len(c) for c in chunks],
        "dim": int(pooled.shape[0]),
    }
    return pooled, stats


def pairwise_metrics(vecs: Dict[str, np.ndarray]) -> Dict[Tuple[str, str], Dict[str, float]]:
    """
    Compute symmetric pairwise metrics for all keys in vecs.
    """
    keys = list(vecs.keys())
    results = {}
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = keys[i], keys[j]
            u, v = vecs[a], vecs[b]
            cos = cosine_similarity(u, v)
            rad, deg = angle_from_cosine(cos)
            eu = euclidean_distance(u, v)
            l1 = manhattan_distance(u, v)
            results[(a, b)] = {
                "cosine_similarity": cos,
                "cosine_distance": 1.0 - cos if not math.isnan(cos) else float('nan'),
                "angle_rad": rad,
                "angle_deg": deg,
                "euclidean": eu,
                "manhattan": l1,
                "dot_raw": float(np.dot(u, v)),
                "norm_u": l2_norm(u),
                "norm_v": l2_norm(v),
            }
    return results


def save_artifacts(outdir: str,
                   records: Dict[str, Dict],
                   vecs: Dict[str, np.ndarray],
                   per_article_stats: Dict[str, Dict],
                   pw: Dict[Tuple[str, str], Dict[str, float]]) -> None:
    os.makedirs(outdir, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    # JSON blob
    blob = {
        "timestamp": ts,
        "n_articles": len(records),
        "articles": records,
        "per_article_stats": per_article_stats,
        "pairwise": {f"{a}__{b}": m for (a, b), m in pw.items()},
        "vector_dim": int(next(iter(vecs.values())).shape[0]) if vecs else None,
    }
    json_path = os.path.join(outdir, f"pmid_similarity_{ts}.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(blob, f, ensure_ascii=False, indent=2)

    # CSV for pairwise
    csv_path = os.path.join(outdir, f"pairwise_{ts}.csv")
    header = ["pmid_a", "pmid_b", "cosine_similarity", "cosine_distance",
              "angle_deg", "euclidean", "manhattan", "dot_raw", "norm_u", "norm_v"]
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(",".join(header) + "\n")
        for (a, b), m in pw.items():
            row = [
                a, b,
                f"{m['cosine_similarity']:.6f}",
                f"{m['cosine_distance']:.6f}" if not math.isnan(m['cosine_distance']) else "nan",
                f"{m['angle_deg']:.3f}",
                f"{m['euclidean']:.6f}",
                f"{m['manhattan']:.6f}",
                f"{m['dot_raw']:.6f}",
                f"{m['norm_u']:.6f}",
                f"{m['norm_v']:.6f}",
            ]
            f.write(",".join(row) + "\n")

    print(f"\nArtifacts saved:\n- {json_path}\n- {csv_path}\n")


def main():
    ap = argparse.ArgumentParser(description="Compute cosine similarity between PubMed articles via LM Studio embeddings.")
    ap.add_argument("--pmids", nargs="+",
                    default=["24548571", "40016490"],
                    help="List of PubMed IDs.")
    ap.add_argument("--pmids-file", type=str, default=None,
                    help="Optional file with one PMID per line. Combined with --pmids if both are given.")
    ap.add_argument("--model", type=str,
                    default=os.getenv("LM_EMBED_MODEL", "text-embedding-qwen3-embedding-0.6b@f16"),
                    help="Embedding model id for LM Studio.")
    ap.add_argument("--lm-url", type=str,
                    default=os.getenv("LM_EMBED_URL", "http://127.0.0.1:1234/v1/embeddings"),
                    help="LM Studio embeddings endpoint URL.")
    ap.add_argument("--max-chars", type=int, default=8000,
                    help="Max characters per chunk before mean-pooling. Set <=0 to disable chunking.")
    ap.add_argument("--include-title", action=argparse.BooleanOptionalAction, default=True,
                    help="Include article title in the embedded text.")
    ap.add_argument("--include-abstract", action=argparse.BooleanOptionalAction, default=True,
                    help="Include article abstract in the embedded text.")
    ap.add_argument("--ncbi-email", type=str, default=None,
                    help="Optional contact email for NCBI E-utilities etiquette.")
    ap.add_argument("--outdir", type=str, default=None,
                    help="If set, write JSON and CSV artifacts here.")

    args = ap.parse_args()

    pmids = list(args.pmids)
    if args.pmids_file:
        with open(args.pmids_file, "r", encoding="utf-8") as f:
            file_pmids = [ln.strip() for ln in f if ln.strip()]
        pmids.extend(file_pmids)
    # Deduplicate while preserving order
    seen = set()
    pmids = [p for p in pmids if (p not in seen and not seen.add(p))]

    print("=== Configuration ===")
    print(f"PMIDs           : {pmids}")
    print(f"Model           : {args.model}")
    print(f"LM URL          : {args.lm_url}")
    print(f"Include         : title={args.include_title} | abstract={args.include_abstract}")
    print(f"Max chars/chunk : {args.max_chars}")
    print(f"NCBI email      : {args.ncbi_email or '(none)'}")
    print()

    # 1) Fetch and parse records
    records: Dict[str, Dict] = {}
    for pmid in pmids:
        try:
            t0 = time.time()
            xml = fetch_pubmed_xml(pmid, email=args.ncbi_email)
            rec = parse_pubmed_record(xml)
            t1 = time.time()
            records[pmid] = rec
            print(f"[Fetch] PMID={pmid} | title_len={len(rec['title'])} | abstract_len={len(rec['abstract'])} | {t1 - t0:.2f}s")
        except Exception as e:
            print(f"[ERROR] Failed to fetch/parse PMID {pmid}: {e}", file=sys.stderr)

    # 2) Build embedding texts
    texts: Dict[str, str] = {}
    for pmid, rec in records.items():
        text = build_text_for_embedding(rec,
                                        include_title=args.include_title,
                                        include_abstract=args.include_abstract)
        if not text:
            print(f"[WARN] PMID {pmid} has empty text under current include flags; skipping.")
            continue
        texts[pmid] = text

    if len(texts) < 2:
        print("[ABORT] Need at least two valid articles to compare.", file=sys.stderr)
        sys.exit(2)

    # 3) Embed with chunking + mean pooling
    vecs: Dict[str, np.ndarray] = {}
    per_article_stats: Dict[str, Dict] = {}
    for pmid, text in texts.items():
        try:
            emb, stats = embed_article_text(text,
                                            model=args.model,
                                            lm_url=args.lm_url,
                                            max_chars=args.max_chars)
            vecs[pmid] = emb
            per_article_stats[pmid] = {
                **stats,
                "vector_stats": describe_vector(emb),
                "char_len": len(text),
                "word_count": len(text.split()),
            }
            print(f"[Embed] PMID={pmid} | chunks={stats['n_chunks']} | dim={stats['dim']} | char_len={len(text)}")
        except Exception as e:
            print(f"[ERROR] Embedding failed for PMID {pmid}: {e}", file=sys.stderr)

    # 4) Pairwise metrics
    print("\n=== Pairwise Metrics ===")
    pw = pairwise_metrics(vecs)
    for (a, b), m in pw.items():
        print(f"{a}  vs  {b}")
        print(f"  cosine_similarity : {m['cosine_similarity']:.6f}")
        print(f"  cosine_distance   : {m['cosine_distance']:.6f}" if not math.isnan(m['cosine_distance']) else "  cosine_distance   : nan")
        print(f"  angle_deg         : {m['angle_deg']:.3f}")
        print(f"  euclidean         : {m['euclidean']:.6f}")
        print(f"  manhattan         : {m['manhattan']:.6f}")
        print(f"  dot_raw           : {m['dot_raw']:.6f}")
        print(f"  ||u||, ||v||      : {m['norm_u']:.6f}, {m['norm_v']:.6f}")
        print()

    # 5) Per-article summary (useful for audits)
    print("=== Per-Article Summary ===")
    for pmid, rec in records.items():
        if pmid not in vecs:
            continue
        vs = per_article_stats[pmid]["vector_stats"]
        print(f"PMID {pmid}")
        print(f"  Title   : {rec['title'][:160]}{'...' if len(rec['title'])>160 else ''}")
        print(f"  Journal : {rec['journal']} ({rec['year']})")
        if rec.get("article_type"):
            print(f"  Types   : {rec['article_type']}")
        print(f"  Text    : chars={per_article_stats[pmid]['char_len']}, words={per_article_stats[pmid]['word_count']}, chunks={per_article_stats[pmid]['n_chunks']}")
        print(f"  Vector  : dim={vs['dim']}, norm(L2)={vs['l2_norm']:.6f}, norm(L1)={vs['l1_norm']:.6f}, mean={vs['mean']:.6f}, std={vs['std']:.6f}, min={vs['min']:.6f}, max={vs['max']:.6f}")
        print()

    # 6) Optional artifacts
    if args.outdir:
        save_artifacts(args.outdir, records, vecs, per_article_stats, pw)


if __name__ == "__main__":
    main()
