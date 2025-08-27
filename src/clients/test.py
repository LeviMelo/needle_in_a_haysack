# src/clients/test.py
from clients.icite import get_pubs, extract_refs_and_citers

if __name__ == "__main__":
    pmid = "33591115"
    recs = get_pubs([pmid])
    r = recs[0] if recs else {}
    refs, citers = extract_refs_and_citers(r)
    print("pmid:", pmid)
    print("refs:", len(refs))
    print("citers:", len(citers))
    print("first few citers:", citers[:10])
