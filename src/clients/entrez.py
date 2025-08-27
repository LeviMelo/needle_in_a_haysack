from __future__ import annotations
import requests
from typing import Dict, List, Any, Iterable, Optional
from urllib.parse import urlencode

# ⬇⬇⬇ change to absolute import (because 'src/' is on sys.path)
from config import ENTREZ_EMAIL, ENTREZ_API_KEY, HTTP_TIMEOUT, USER_AGENT

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEADERS = {"User-Agent": USER_AGENT, "Accept": "application/json"}

def esearch(query: str, db: str = "pubmed", retmax: int = 10000, mindate: Optional[int]=None, maxdate: Optional[int]=None, sort: str="date") -> List[str]:
    params = {
        "db": db,
        "term": query,
        "retmode": "json",
        "retmax": retmax,
        "sort": sort,
        "email": ENTREZ_EMAIL
    }
    if ENTREZ_API_KEY:
        params["api_key"] = ENTREZ_API_KEY
    if mindate:
        params["mindate"] = str(mindate)
    if maxdate:
        params["maxdate"] = str(maxdate)
    r = requests.get(f"{EUTILS}/esearch.fcgi", headers=HEADERS, params=params, timeout=HTTP_TIMEOUT)
    r.raise_for_status()
    return r.json().get("esearchresult", {}).get("idlist", [])

def esummary(pmids: Iterable[str]) -> Dict[str, Dict[str,Any]]:
    pmids = list(pmids)
    out: Dict[str,Dict[str,Any]] = {}
    for i in range(0, len(pmids), 500):
        chunk = pmids[i:i+500]
        params = {
            "db":"pubmed", "retmode":"json", "id": ",".join(chunk),
            "email": ENTREZ_EMAIL
        }
        if ENTREZ_API_KEY: params["api_key"] = ENTREZ_API_KEY
        r = requests.get(f"{EUTILS}/esummary.fcgi", headers=HEADERS, params=params, timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        data = r.json().get("result", {})
        for k,v in data.items():
            if k == "uids": continue
            out[k] = v
    return out

def efetch_abstracts(pmids: Iterable[str]) -> Dict[str, Dict[str,Any]]:
    pmids = [str(p) for p in pmids]
    out: Dict[str,Dict[str,Any]] = {}

    import re
    import xml.etree.ElementTree as ET

    def _join_itertext(node) -> str:
        """Safely join text from an Element (including nested tags)."""
        if node is None:
            return ""
        try:
            return "".join(node.itertext())
        except AttributeError:
            # Defensive: if someone accidentally passes a list
            if isinstance(node, list):
                parts = []
                for n in node:
                    try:
                        parts.append("".join(n.itertext()))
                    except Exception:
                        parts.append((getattr(n, "text", None) or ""))
                return "".join(parts)
            return (getattr(node, "text", None) or "")

    for i in range(0, len(pmids), 200):
        chunk = pmids[i:i+200]
        params = {
            "db": "pubmed",
            "retmode": "xml",
            "rettype": "abstract",
            "id": ",".join(chunk),
            "email": ENTREZ_EMAIL
        }
        if ENTREZ_API_KEY:
            params["api_key"] = ENTREZ_API_KEY

        r = requests.get(f"{EUTILS}/efetch.fcgi",
                         headers={"User-Agent": USER_AGENT},
                         params=params,
                         timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        root = ET.fromstring(r.text)

        for art in root.findall(".//PubmedArticle"):
            pmid = art.findtext(".//PMID") or ""

            # Title (handles inline formatting tags)
            title_el = art.find(".//ArticleTitle")
            title = _join_itertext(title_el).strip()

            # Abstract (structured abstracts have multiple <AbstractText> nodes)
            abs_nodes = art.findall(".//Abstract/AbstractText")
            abstract = " ".join(_join_itertext(n).strip() for n in abs_nodes) if abs_nodes else ""

            # Year: try several places, grab first 4-digit year
            year = None
            for path in (".//ArticleDate/Year",
                         ".//PubDate/Year",
                         ".//DateCreated/Year",
                         ".//PubDate/MedlineDate"):
                s = art.findtext(path)
                if s:
                    m = re.search(r"\d{4}", s)
                    if m:
                        year = int(m.group(0))
                        break

            journal = art.findtext(".//Journal/Title") or ""
            pubtypes = [pt.text for pt in art.findall(".//PublicationTypeList/PublicationType") if pt.text]

            # DOI
            doi = None
            for idn in art.findall(".//ArticleIdList/ArticleId"):
                if (idn.attrib.get("IdType","").lower() == "doi") and idn.text:
                    doi = idn.text.strip().lower()

            out[pmid] = {
                "pmid": pmid,
                "title": title,
                "abstract": abstract,
                "year": year,
                "pub_types": pubtypes,
                "doi": doi,
                "journal": journal
            }

    return out

