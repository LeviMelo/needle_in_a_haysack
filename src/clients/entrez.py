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
    pmids = list(pmids)
    out: Dict[str,Dict[str,Any]] = {}
    for i in range(0, len(pmids), 200):
        chunk = pmids[i:i+200]
        params = {
            "db":"pubmed", "retmode":"xml", "rettype":"abstract", "id": ",".join(chunk),
            "email": ENTREZ_EMAIL
        }
        if ENTREZ_API_KEY: params["api_key"] = ENTREZ_API_KEY
        r = requests.get(f"{EUTILS}/efetch.fcgi", headers={"User-Agent": USER_AGENT}, params=params, timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        # Minimal XML parsing to extract Title, Abstract, PublicationTypes, Year, DOI (if present)
        import xml.etree.ElementTree as ET
        root = ET.fromstring(r.text)
        for art in root.findall(".//PubmedArticle"):
            pmid = art.findtext(".//PMID")
            title = art.findtext(".//ArticleTitle") or ""
            abst = " ".join([n.text or "" for n in art.findall(".//AbstractText")]) or ""
            year = None
            dp = art.findtext(".//PubDate/Year") or art.findtext(".//PubDate/MedlineDate")
            journal = art.findtext(".//Journal/Title") or ""
            try:
                year = int(dp[:4]) if dp else None
            except Exception:
                year = None
            pubtypes = [pt.text for pt in art.findall(".//PublicationTypeList/PublicationType") if pt.text]
            doi = None
            for idn in art.findall(".//ArticleIdList/ArticleId"):
                if idn.attrib.get("IdType","").lower() == "doi":
                    doi = (idn.text or "").lower()
            out[pmid] = {"pmid": pmid, "title": title, "abstract": abst, "year": year, "pub_types": pubtypes, "doi": doi, "journal": journal}
    return out
