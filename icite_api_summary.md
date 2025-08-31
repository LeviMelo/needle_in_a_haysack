Here’s a practical, end-to-end briefing on the NIH iCite API—what it is, exactly which data you get (and don’t), how to query it efficiently, and how to integrate it with PubMed/E-utilities, Crossref, and NIH RePORTER in real research pipelines.

# What iCite is (in one sentence)

iCite is the NIH Office of Portfolio Analysis’ public web service for article-level bibliometrics and open citation links for PubMed-indexed journal articles, built on the NIH Open Citation Collection (NIH-OCC). ([icite.od.nih.gov][1], [PLOS][2], [PMC][3])

---

# What data iCite contains (modules → fields → where they come from)

## 1) Influence (RCR and related metrics)

* **Core idea:** field- and time-normalized influence of an article, via the **Relative Citation Ratio (RCR)**. Also exposes underlying terms used to compute RCR. ([support.icite.nih.gov][4])
* **Key fields (API “legacy=false” names in parentheses):**

  * `pmid` (`pmid`), `title` (`title`), `authors` (`authors`), `journal` ISO abbrev (`journalNameIso`), `pubYear` (`pubYear`), `doi` (`doi`).
  * **RCR** (`rcr`), **Article Citation Rate**—citations/year (`acr`), **Field Citation Rate** (`fcr`), **Expected Citation Rate** (`ecr`), NIH RCR percentile (`nihRcrPercentile`), **provisional** flag for recent papers (`rcrIsProvisional`). ([support.icite.nih.gov][4])
* **What RCR means:** RCR≈1.0 ≈ median NIH-funded paper in same field/year; 2.0 ≈ twice as influential, etc. Calculation details in Hutchins 2016. ([support.icite.nih.gov][4], [PMC][5])

## 2) Translation (bench-to-bedside orientation)

* **Core idea:** where the paper sits on the **Triangle of Biomedicine** (Human / Animal / Molecular-Cellular) and **APT** (Approximate Potential to Translate—probability the paper will later be cited by clinical trials/guidelines). ([support.icite.nih.gov][6])
* **Key fields:**

  * Human/Animal/Mol-Cell fractions (`human`, `animal`, `molCell`) derived from MeSH categories; triangle coordinates (`xCoord`, `yCoord`).
  * **APT** (`apt`) ∈ {0.05, 0.25, 0.50, 0.75, 0.95}. ([support.icite.nih.gov][6])
  * Flags for **isClinicalArticle** and list of **citing clinical PMIDs** (trials, guidelines). ([support.icite.nih.gov][6])

## 3) Citations (link-level “who cites whom”)

* **Core idea:** Open, public-domain citation links for PubMed-indexed articles; source aggregation from **Crossref, Medline, PubMed Central, and Entrez**, plus ML extraction from open-access full texts. ([support.icite.nih.gov][7])
* **Key fields:**

  * Total citations (`citedByPmidCount`), **cited-by PMIDs** (`citedByPmids`), **references / cited PMIDs** (`citedPmids`), and sometimes **by-year breakdown** (`citedByPmidsByYear`). ([support.icite.nih.gov][8])

> **Where it all comes from:** The underlying **NIH-OCC** is described in Hutchins et al., PLoS Biology (2019); it consolidates open citations from Medline, PMC, Crossref, and ML-parsed open-access full texts, and iCite disseminates the dataset via the web service and periodic bulk snapshots. ([PMC][3], [support.icite.nih.gov][7])

---

# What iCite does **not** contain (and common misconceptions)

* **No abstracts or full text.** iCite provides concise metadata and metrics—not abstracts or article bodies. (Use **PubMed E-utilities** for titles/abstracts/MeSH terms.) ([Biblioteca Nacional de Medicina][9])
* **No full MeSH heading lists per article.** iCite exposes **derived** Human/Animal/Mol-Cell fractions but does **not** return the article’s MeSH term list. ([support.icite.nih.gov][6])
* **No author disambiguation profiles** (author IDs), affiliations, funding links, or grant data. For funding/grants, use **NIH RePORTER API**; to enrich authors, rely on other sources (e.g., PubMed XML, ORCID, OpenAlex, Scopus). ([api.reporter.nih.gov][10])
* **Coverage is PubMed-indexed journal articles;** nodes are PubMed PMIDs. Citations are gathered from public sources but reported as links **among PubMed articles**; journals outside PubMed aren’t represented as nodes in iCite. ([itools.od.nih.gov][11], [support.icite.nih.gov][7])
* **Temporal coverage:** the iCite database is focused on **1980–present**; older literature is generally out of scope, which matters for longitudinal fields. ([icite.od.nih.gov][12], [ophthalmologyretina.org][13])
* **Web UI query caps vs. API caps:** the **web interface** search/analysis is capped (e.g., typically 10,000 PMIDs; file uploads currently mention up to 50,000), but the **API** has **per-request** limits (e.g., `pmids` ≤ 1000 per request, `limit` ≤ 1000). Don’t conflate UI limits with API limits. ([icite.od.nih.gov][12], [support.icite.nih.gov][14])
* **Rate limits:** iCite support pages don’t publish hard rate-limit numbers; heavy users should **batch**, **cache**, and prefer **bulk snapshots** when feasible. (Snapshots are released on NIH Figshare.) ([support.icite.nih.gov][15], [nih.figshare.com][16])

---

# The API surface (what you can actually call)

## Base concept

A single REST collection: **`/api/pubs`** returning publication-level records that include Influence, Translation, and Citations fields. Response can be **JSON** (default) or **CSV**. You can request specific PMIDs (up to 1000 per call), filter by **year**, paginate using **offset/limit**, select fields via **`fl`**, and choose **legacy vs. database** field names. ([support.icite.nih.gov][15])

### Endpoints & parameters (current support page)

* `GET /api/pubs/{pmid}` → one publication.
* `GET /api/pubs?pmids=PMID1,PMID2,…` (≤1000) → multiple publications.
* `GET /api/pubs?year=YYYY&limit=…&offset=…` → browse by year.
* `&fl=pmid,year,title,rcr,apt,citedByPmidCount,…` → **field selector**.
* `&format=csv` → CSV output (also honored if sent via `Accept: text/csv`).
* `&legacy=false` → switch to database field names (e.g., `rcr` instead of `relative_citation_ratio`).
  Examples are documented in the official support page. ([support.icite.nih.gov][15])

> **Historical note:** Early release notes mentioned a `&refs=true` toggle for exposing citation links. The current API examples already include `citedByPmids` and `citedPmids` in default responses; rely on the present **Bulk Data & API** documentation for the latest behavior. ([itools.od.nih.gov][17], [support.icite.nih.gov][15])

### Response shape (abridged real example)

* JSON includes: `pmid`, `title`, `authors`, `journal`/`journalNameIso`, `pubYear`, `doi`, `rcr`, `nihRcrPercentile`, `acr`, `fcr`, `ecr`, `rcrIsProvisional`, `human`, `animal`, `molCell`, `xCoord`, `yCoord`, `apt`, `isClinicalArticle`, `citingClinicalPmids`, `citedByPmidCount`, `citedByPmidsByYear`, `citedByPmids`, `citedPmids`. ([support.icite.nih.gov][15])

---

# Bulk access (snapshots)

For large-scale analytics or to avoid request caps, NIH publishes **iCite database snapshots** (and a dedicated **citation-links CSV**) on **NIH Figshare**. These snapshots include the NIH-OCC link set and iCite metadata in CSV and JSON tarballs. ([nih.figshare.com][18])

---

# Coverage, sources, and known constraints

* **Nodes (articles):** PubMed-indexed journal articles; iCite integrates with **PubMed** for queries and identifiers. ([icite.od.nih.gov][12])
* **Edges (citations):** Aggregated from **Crossref, Medline, PMC, Entrez** and **ML extraction** of open-access full texts; links are public-domain. Keep in mind: coverage is very broad but still imperfect (e.g., non-OA reference lists missing, publisher metadata lag, pre-1980 gaps). ([support.icite.nih.gov][7])
* **Time window:** Focused on **1980–present** (RCR available across that window). This influences field baselines and longitudinal comparisons. ([icite.od.nih.gov][12])
* **Provisional RCR:** Recent-year RCRs are flagged **provisional** as citation counts stabilize. Treat early RCRs with caution in evaluations. ([support.icite.nih.gov][4])

---

# How to use iCite effectively (researcher’s playbook)

## A) Fetch metrics and open citations for known PMIDs

* **Python (requests):**

  ```python
  import requests
  pmids = ["23456789", "27599104"]
  url = "https://icite.od.nih.gov/api/pubs"
  r = requests.get(url, params={
      "pmids": ",".join(pmids),
      "legacy": "false",                 # use modern field names (rcr, fcr, ecr…)
      "fl": "pmid,title,doi,pubYear,journalNameIso,rcr,acr,fcr,ecr,nihRcrPercentile,"
            "human,animal,molCell,xCoord,yCoord,apt,isClinicalArticle,citingClinicalPmids,"
            "citedByPmidCount,citedByPmids,citedPmids"
  })
  r.raise_for_status()
  data = r.json()
  ```

  Returns a list of per-PMID dicts with metrics and citation links. ([support.icite.nih.gov][15])

* **CSV one-liner (for spreadsheets/SQL ingest):**
  `https://icite.od.nih.gov/api/pubs?pmids=23456789,27599104&format=csv` ([icite.od.nih.gov][19])

## B) Build neighborhoods / citation graphs

Use `citedPmids` (references) and `citedByPmids` (citers) to assemble local networks; `citedByPmidsByYear` helps time-slice growth. For very large ego-nets, prefer **bulk snapshots**. ([support.icite.nih.gov][15], [nih.figshare.com][18])

## C) Benchmark influence & translation

* Compare **RCR** (field-normalized) instead of raw citation counts; consider **NIH percentiles** for interpretability. ([support.icite.nih.gov][4])
* Use **APT** and **citingClinicalPmids** to flag translation toward clinical trials/guidelines. ([support.icite.nih.gov][6])
* Map **Human/Animal/Mol-Cell** fractions and triangle coordinates for portfolio positioning. ([support.icite.nih.gov][6])

## D) Integrate with complementary APIs

* **PubMed E-utilities**: search/query by topic, pull titles/abstracts/MeSH, author affiliations, and other bibliographic details—then resolve PMIDs into iCite. ([Biblioteca Nacional de Medicina][9])
* **NIH RePORTER API**: link PMIDs to NIH grants and funding portfolios. ([api.reporter.nih.gov][10])
* **Crossref / publisher data**: fill DOI gaps, fetch issue/volume, publisher metadata (when needed). (iCite already brings DOIs when present.) ([support.icite.nih.gov][15])

---

# Practical limits, batching, and scaling

* **Per-request constraints (API):** `pmids` ≤ 1000, `limit` ≤ 1000; paginate large pulls with `offset` (by PMID), or iterate by `year`. ([support.icite.nih.gov][15])
* **Web UI constraints:** searching/analysis caps (e.g., 10k PMIDs) and bulk **uploads** up to 50k PMIDs; for bigger jobs or repeated runs, switch to snapshots. ([icite.od.nih.gov][12], [support.icite.nih.gov][14])
* **Rate limiting:** iCite doesn’t document explicit numeric limits; to be robust, **cache responses**, **exponential-backoff** on HTTP 429/5xx, and **prefer snapshots** for heavy network construction. ([support.icite.nih.gov][15], [nih.figshare.com][16])

---

# Data quality cautions (what to watch for)

1. **Edge completeness varies by openness**
   References present in closed full texts may be missing; Crossref/PMC coverage improves it but doesn’t guarantee exhaustiveness. Use OCC snapshots to check coverage. ([support.icite.nih.gov][7])

2. **Field normalization context (RCR)**
   RCR is article-level and field/time-normalized via **co-citation networks**; this is stronger than journal-based metrics but still sensitive to field idiosyncrasies and early-year volatility (hence “provisional”). ([support.icite.nih.gov][4])

3. **Scope bias (1980–present)**
   For mature literatures with pre-1980 classics, RCR and citation counts will miss early foundational edges. Interpret trends with the window in mind. ([icite.od.nih.gov][12])

---

# Side-by-side: iCite vs. adjacent sources (what to use when)

* **iCite API** → *influence* (RCR, ACR/FCR/ECR), *translation* (APT; clinical citations), and *open citation links* for PMIDs; minimal descriptive metadata. ([support.icite.nih.gov][15])
* **PubMed E-utilities** → rich bibliographic data (titles, abstracts, MeSH, affiliations); search construction; PMIDs discovery. Use first to **find** and then feed PMIDs to iCite. ([Biblioteca Nacional de Medicina][9])
* **NIH RePORTER API** → grant-publication linkages, funding info; join on PMID or grant IDs. ([api.reporter.nih.gov][10])
* **Crossref** → DOI registries, reference lists (varying openness), publisher metadata—useful to complement gaps when you need non-PubMed nodes.

---

# Example workflows

### 1) Topic → PMIDs → iCite metrics + links

1. Use PubMed E-utilities to query by PICO/MeSH → collect PMIDs. ([Biblioteca Nacional de Medicina][9])
2. Call `GET /api/pubs?pmids=…&legacy=false&fl=pmid,title,pubYear,rcr,apt,citedByPmidCount,citedByPmids,citedPmids`. ([support.icite.nih.gov][15])
3. Build an ego-net from `citedPmids`/`citedByPmids`, then compute local indicators (new trial-cited neighbors, high-RCR gaps, etc.).

### 2) Portfolio audit (department, grant set)

1. Start from an internal PMID list or NIH RePORTER grant set, de-duplicate; ([api.reporter.nih.gov][10])
2. Pull iCite: RCR distributions, APT spectrum, % cited by clinical, and a clinical-citation list (`citingClinicalPmids`) for translational impact. ([support.icite.nih.gov][6])

### 3) Heavy graph analysis

Download the **NIH-OCC snapshot** (CSV/JSON) and the **citations-only CSV** from Figshare; load into your graph DB / Python network stack for scalable traversal. ([nih.figshare.com][18])

---

# Field reference (selected, with modern names)

* **Identity & bibliographic:** `pmid`, `doi`, `title`, `authors`, `journalNameIso`, `pubYear`. ([support.icite.nih.gov][15])
* **Influence:** `rcr`, `acr` (citations/year), `fcr`, `ecr`, `nihRcrPercentile`, `rcrIsProvisional`. ([support.icite.nih.gov][4])
* **Translation:** `human`, `animal`, `molCell`, `xCoord`, `yCoord`, `apt`, `isClinicalArticle`, `citingClinicalPmids`. ([support.icite.nih.gov][6])
* **Citations:** `citedByPmidCount`, `citedByPmids`, `citedPmids`, `citedByPmidsByYear`. ([support.icite.nih.gov][8])

---

# Notes on the web interface (if you also use it)

* iCite’s **Analysis** UI runs your PubMed query (through PubMed’s API) and then analyzes up to **10,000 PMIDs** (recently, uploads note up to **50,000** PMIDs). Useful for quick checks, but not a substitute for API/snapshots at scale. ([icite.od.nih.gov][12], [support.icite.nih.gov][14])
* The UI shows the three module tabs (Influence/Translation/Citations) with similar definitions as the API. ([info.hsls.pitt.edu][20])

---

# Bottom line

* Use **iCite** when you need **RCR**, **translation signals (APT / clinical citations)**, and **open citation links** among **PubMed** articles. It’s lean, fast, and public-domain for the link data. ([support.icite.nih.gov][7])
* Pair it with **PubMed E-utilities** for abstracts/MeSH and **NIH RePORTER** for grants—this three-piece stack covers most biomedical evidence-mapping pipelines. ([Biblioteca Nacional de Medicina][9], [api.reporter.nih.gov][10])
* For **large graphs** or repeated pulls, consume **Figshare snapshots** rather than hammering the API. ([nih.figshare.com][16])

# Example

## https://icite.od.nih.gov/api/pubs/32162448
```json
{
  "data": [
    {
      "_id": "32162448",
      "animal": 0,
      "authors": [
        {
          "lastName": "Zhang",
          "fullName": "Zhang, Yan-Ling",
          "firstName": "Yan-Ling"
        },
        {
          "lastName": "Li",
          "fullName": "Li, Hui",
          "firstName": "Hui"
        },
        {
          "lastName": "Zeng",
          "fullName": "Zeng, Hua",
          "firstName": "Hua"
        },
        {
          "lastName": "Li",
          "fullName": "Li, Qiao",
          "firstName": "Qiao"
        },
        {
          "lastName": "Qiu",
          "fullName": "Qiu, Li-Ping",
          "firstName": "Li-Ping"
        },
        {
          "lastName": "Dai",
          "fullName": "Dai, Ru-Ping",
          "firstName": "Ru-Ping"
        }
      ],
      "doi": "10.1111/pan.13853",
      "human": 1,
      "pmid": 32162448,
      "title": "Ultrasonographic evaluation of gastric emptying after ingesting carbohydrate-rich drink in young children: A randomized crossover study.",
      "apt": 0.75,
      "citedByClinicalArticle": true,
      "citedByPmidsByYear": [
        {
          "32997821": 2020
        },
        {
          "35874585": 2022
        },
        {
          "40535648": 2025
        },
        {
          "36891367": 2023
        },
        {
          "39773771": 2025
        },
        {
          "40104681": 2025
        },
        {
          "37601499": 2023
        },
        {
          "38993713": 2024
        },
        {
          "38615661": 2024
        },
        {
          "34857683": 2022
        },
        {
          "36629465": 2023
        },
        {
          "39949966": 2024
        },
        {
          "34403386": 2021
        },
        {
          "37477770": 2023
        },
        {
          "39651668": 2025
        },
        {
          "35166373": 2022
        },
        {
          "38350987": 2024
        }
      ],
      "year": 2020,
      "journal": "Paediatr Anaesth",
      "is_research_article": true,
      "citation_count": 17,
      "field_citation_rate": 2.70128652829193,
      "expected_citations_per_year": 1.55624071151373,
      "citations_per_year": 3.4,
      "relative_citation_ratio": 2.18475199552701,
      "nih_percentile": 77.1,
      "molecular_cellular": 0,
      "x_coord": 0,
      "y_coord": 1,
      "is_clinical": true,
      "cited_by_clin": [32997821, 38615661, 35166373],
      "cited_by": [32997821, 35874585, 40535648, 36891367, 39773771, 40104681, 37601499, 38993713, 38615661, 34857683, 36629465, 39949966, 34403386, 37477770, 39651668, 35166373, 38350987],
      "references": [30500064, 27595380, 25940831, 28045707, 17304411, 23302981, 29265174, 20108579, 22277664, 29700894, 21676893, 21596885, 22729296, 15900642, 24497851, 22235107, 20177856, 25501720, 27106968, 24189352, 29077816, 17902094, 21999211, 26514824, 15739210, 30115264, 25495405, 11682427, 28675597, 25260696, 29077812],
      "provisional": false
    }
  ]
}
```

# Closing Notes on iCite (Nuances & Caveats)

## 1. Publication Type Logic (Primary Research vs Clinical)

* **Primary Research (`is_research_article`):**
  iCite applies **whitelist/blacklist rules** on PubMed Publication Types:

  * Qualifying: *Journal Article, Clinical Trial (all phases), RCT, Observational Study, Meta-Analysis, Case Report, Dataset, Preprint, etc.*
  * Disqualifying: *Review, Editorial, Guideline, Practice Guideline, Comment, Erratum, Retraction, etc.*
    Result: A **Systematic Review with “Meta-Analysis” tag** → `is_research_article=true`. A **Narrative Review** → `false`.

* **Clinical Article (`isClinicalArticle`):**
  Separate rule set: flagged *true* if PubMed types include *Clinical Trial (any phase), Clinical Study, Observational Study, Guideline, Practice Guideline*.
  Thus: a *Guideline* is **clinical** but not **primary research**.

* **Cited-by clinical vs. being clinical:**
  Distinction is crucial: `isClinicalArticle=false` + `citedByClinicalArticle=true` is common (e.g., a meta-analysis not itself a trial but cited in many trials/guidelines). This is exactly the scenario that **APT** predicts.

---

## 2. Relative Citation Ratio (RCR) and Benchmarks

* **Field definition:** RCR fields are *not* journal categories but **co-citation neighborhoods**; this better adapts to multidisciplinary articles but requires citers to exist. With **no citers**, `fcr/ecr/rcr=null`.
* **NIH-funded benchmark:**

  * Originally (2016): NIH **R01-funded** corpus.
  * Current production: **all NIH-funded publications** (per OPA docs).
  * Implication: RCR values are tied to NIH’s funding portfolio, not global science.
* **Regression model:** Each year, NIH fits ACR \~ FCR across NIH-funded papers → `ecr = α + β·fcr`. By design, **median NIH-funded article has RCR≈1.0**. This calibration is invisible to the user but underpins comparability.
* **Provisional rules:**

  * No RCR in the current fiscal year unless ≥5 citations.
  * Previous two years flagged `rcrIsProvisional=true`.
  * Fiscal year cutoff: October.

---

## 3. Approximate Potential to Translate (APT)

* **Definition:** ML-predicted probability the paper will be cited by a *clinical trial or guideline*.
* **Reporting:** Discrete bins {0.05, 0.25, 0.50, 0.75, 0.95}.
* **Features:** Paper’s **HAMC profile**, early citing network properties.
* **Interpretation:** \~x% of papers with APT=x will eventually be cited by clinical articles. It is a **population-level probability**, not an individual guarantee.
* **Common confusion:** APT is **not suppressed by RCR’s provisional rules**. A brand-new article with no clinical signals → APT=0.05.
* **Key nuance:** APT can be high (0.95) while `isClinicalArticle=false`. Why? Because the model predicts it will be *used by* clinical trials/guidelines, not that it *is* one.

---

## 4. Human/Animal/Molecular Fractions (HAMC)

* Derived from PubMed **MeSH** indexing.
* Fractions sum to 1.0 across `human`, `animal`, `molCell`.
* Mapped into Cartesian coordinates (`xCoord`, `yCoord`) for visualization in the **Triangle of Biomedicine**.
* Useful for portfolio analyses (bench → bedside positioning) but beware of MeSH indexing lags/errors.

---

## 5. Citations Module — Edges and Completeness

* Sources: Medline, PMC, Crossref, Entrez, + ML parsing of open-access full texts → merged into **NIH-OCC**.
* **Biases:**

  * Closed-access journals with no open references = gaps.
  * Citations reported are **PubMed-to-PubMed links**; non-PubMed nodes (books, arXiv-only preprints, some conference papers) aren’t represented.
* **Fields:** `citedByPmids` (citers), `citedPmids` (references), `citedByPmidsByYear` (timeline).
* **Bulk data:** NIH publishes full OCC snapshots on Figshare (recommended for graph-scale work).

---

## 6. Coverage Constraints

* **Time window:** 1980–present. Pre-1980 classics won’t appear as nodes or citers.
* **RCR availability:** Same 1980+ window; expect nulls for earlier PMIDs.
* **Granularity:** NIH-OCC citations are article-level, not section- or sentence-level. No citation context or sentiment.

---

## 7. Practical API Nuances

* **API limits:** `pmids` ≤1000, `limit` ≤1000 per call. Use `offset` for paging.
* **Web UI caps:** up to 10k PMIDs in search, uploads up to 50k.
* **Rate limits:** undocumented; assume conservative policies (batch, cache, prefer bulk).
* **Legacy fields:** by default, API returns legacy names (`relative_citation_ratio`). Use `legacy=false` for compact modern names (`rcr`, `acr`, `fcr`).
* **CSV option:** `format=csv` is often easier for pipelines.

---

## 8. Data Quality Considerations

* **Indexing lag:** PubMed may assign Publication Types months later; iCite inherits this. Flags like `isClinicalArticle` may be unstable for very recent items.
* **Provisional RCR volatility:** Especially in fast-moving fields (COVID, cancer), early citations distort co-citation neighborhoods.
* **APT floor effects:** Many new papers sit at APT=0.05 until enough citing structure emerges.
* **Missing edges:** Because non-PubMed nodes are excluded, edges out of PubMed (to books, conference proceedings, some datasets) vanish. This slightly shrinks “true” degree.

---

## 9. Integration Patterns (best practices)

* **Discovery**: Always use PubMed E-utilities (search/abstracts/MeSH) → feed PMIDs into iCite.
* **Funding**: NIH RePORTER for grant ↔ PMID joins.
* **DOIs/metadata**: Crossref to backfill non-PubMed data.
* **Graphs**: Build networks with OCC bulk CSVs rather than API-hitting at scale.

---

## 10. Conceptual Warnings

* **RCR ≠ impact factor.** RCR is article-level, co-citation normalized, NIH-benchmarked; Impact Factor is journal-level, 2-year window.
* **APT ≠ clinical benefit.** It only signals *likelihood of clinical citation*, not trial outcomes or practice change.
* **Clinical vs. translational:** `isClinicalArticle` (what it *is*) vs. `citedByClinicalArticle` (how it’s *used* clinically). These must be read distinctly.
* **1980 cutoff bias:** Watch for truncated reference chains in long-lived literatures (e.g., cardiology, genetics).


[1]: https://icite.od.nih.gov/?utm_source=chatgpt.com "iCite | Home"
[2]: https://journals.plos.org/plosbiology/article?id=10.1371%2Fjournal.pbio.3000385&utm_source=chatgpt.com "The NIH Open Citation Collection: A public access, broad ..."
[3]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6786512/?utm_source=chatgpt.com "The NIH Open Citation Collection: A public access, broad ..."
[4]: https://support.icite.nih.gov/hc/en-us/related/click?data=BAh7CjobZGVzdGluYXRpb25fYXJ0aWNsZV9pZGwrCBv1fyt%2FCDoYcmVmZXJyZXJfYXJ0aWNsZV9pZGwrCJssKpCHCDoLbG9jYWxlSSIKZW4tdXMGOgZFVDoIdXJsSSJSL2hjL2VuLXVzL2FydGljbGVzLzkzNDIyODM2NzQ5MDctRGVzY3JpcHRpb25zLW9mLUluZmx1ZW5jZS1Nb2R1bGUtZGF0YS1maWVsZHMGOwhUOglyYW5raQo%3D--8f592cb1a65056847cdd3b1103a67513b12e8965 "Descriptions of Influence Module data fields – iCite Support"
[5]: https://pmc.ncbi.nlm.nih.gov/articles/PMC5012559/?utm_source=chatgpt.com "Relative Citation Ratio (RCR): A New Metric That Uses ..."
[6]: https://support.icite.nih.gov/hc/en-us/related/click?data=BAh7CjobZGVzdGluYXRpb25fYXJ0aWNsZV9pZGwrCBuv9j9%2FCDoYcmVmZXJyZXJfYXJ0aWNsZV9pZGwrCBv1fyt%2FCDoLbG9jYWxlSSIKZW4tdXMGOgZFVDoIdXJsSSJUL2hjL2VuLXVzL2FydGljbGVzLzkzNDI2MjcwMDAwOTEtRGVzY3JpcHRpb25zLW9mLVRyYW5zbGF0aW9uLU1vZHVsZS1kYXRhLWZpZWxkcwY7CFQ6CXJhbmtpCQ%3D%3D--4ee5cf64206a49d95f4708e9af89e6fc5af9f8b3 "Descriptions of Translation Module data fields – iCite Support"
[7]: https://support.icite.nih.gov/hc/en-us/articles/9066501755291 "Citations Module Overview – iCite Support"
[8]: https://support.icite.nih.gov/hc/en-us/related/click?data=BAh7CjobZGVzdGluYXRpb25fYXJ0aWNsZV9pZGwrCBuR%2FoBGCDoYcmVmZXJyZXJfYXJ0aWNsZV9pZGwrCJvRnfU%2BCDoLbG9jYWxlSSIKZW4tdXMGOgZFVDoIdXJsSSJSL2hjL2VuLXVzL2FydGljbGVzLzkwOTg5MDQ4OTk4NjctRGVzY3JpcHRpb25zLW9mLUNpdGF0aW9ucy1Nb2R1bGUtZGF0YS1maWVsZHMGOwhUOglyYW5raQY%3D--a48857119cbeffc451b1b24061965e6c67f42c04 "Descriptions of Citations Module data fields – iCite Support"
[9]: https://www.nlm.nih.gov/dataguide/eutilities/utilities.html?utm_source=chatgpt.com "The 9 E-utilities and Associated Parameters"
[10]: https://api.reporter.nih.gov/?utm_source=chatgpt.com "NIH RePORTER API"
[11]: https://itools.od.nih.gov/icite/user_guide?page_id=ug_infl&utm_source=chatgpt.com "iCite | User Guide | NIH Office of Portfolio Analysis"
[12]: https://icite.od.nih.gov/analysis?utm_source=chatgpt.com "iCite | Analysis"
[13]: https://www.ophthalmologyretina.org/article/S2468-6530%2823%2900002-7/fulltext?utm_source=chatgpt.com "Evaluation of Research Productivity among Academic ..."
[14]: https://support.icite.nih.gov/hc/en-us/articles/9002812384027-Upload-a-spreadsheet-or-text-file-of-PMIDs?utm_source=chatgpt.com "Upload a spreadsheet or text file of PMIDs"
[15]: https://support.icite.nih.gov/hc/en-us/articles/9513563045787-Bulk-Data-and-API "Bulk Data and API – iCite Support"
[16]: https://nih.figshare.com/collections/iCite_Database_Snapshots_NIH_Open_Citation_Collection_/4586573?utm_source=chatgpt.com "iCite Database Snapshots (NIH Open Citation Collection)"
[17]: https://itools.od.nih.gov/icite/user_guide?page_id=ug_release_notes&utm_source=chatgpt.com "iCite | User Guide | NIH Office of Portfolio Analysis"
[18]: https://nih.figshare.com/collections/iCite_Database_Snapshots_NIH_Open_Citation_Collection_/4586573/53?utm_source=chatgpt.com "iCite Database Snapshots (NIH Open Citation Collection)"
[19]: https://icite.od.nih.gov/api/pubs?format=csv&pmids=23456789%2C27599104&utm_source=chatgpt.com "https://icite.od.nih.gov/api/pubs?pmids=23456789,2..."
[20]: https://info.hsls.pitt.edu/updatereport/2020/january-2020/assess-article-influence-and-citation-data-with-nih-icite/?utm_source=chatgpt.com "Assess Article Influence and Citation Data with NIH iCite - HSLS - University of Pittsburgh"
