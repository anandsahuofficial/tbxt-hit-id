"""
Onepot library membership filter — to be run on May 9 once organizers provide the lookup interface.

This is a SHIM. The actual interface is unknown until 1:30 pm announcements. Most likely:
  Option A: Web UI — copy SMILES list, paste, get back in-library subset
  Option B: REST API — POST SMILES, GET in-library boolean
  Option C: File-based — provided as bulk catalog SMILES file we filter against locally

This script implements Option C (the most general) and stubs A/B for adaptation on the day.
"""
import argparse
import csv
import sys
from pathlib import Path

from rdkit import Chem


def canonicalize(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None: return None
    return Chem.MolToSmiles(m)


def filter_against_catalog_file(candidates_smiles, catalog_path):
    """Option C: filter against a catalog file (one SMILES per line, optionally with id)."""
    print(f"  Loading onepot catalog from {catalog_path}...")
    catalog = set()
    with open(catalog_path) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line: continue
            smi = line.split()[0]  # in case there's id
            canon = canonicalize(smi)
            if canon: catalog.add(canon)
            if (i+1) % 100000 == 0:
                print(f"    ...{i+1} catalog SMILES loaded")
    print(f"  Catalog size: {len(catalog):,}")

    in_library = []
    not_in = []
    for cid, smi in candidates_smiles:
        c = canonicalize(smi)
        if c in catalog:
            in_library.append((cid, smi))
        else:
            not_in.append((cid, smi))
    return in_library, not_in


def filter_against_api(candidates_smiles, api_url):
    """Option B: query API. Stub — fill in actual API on the day."""
    import requests  # noqa
    in_library = []
    not_in = []
    for cid, smi in candidates_smiles:
        try:
            r = requests.post(api_url, json={"smiles": smi}, timeout=10)
            r.raise_for_status()
            data = r.json()
            if data.get("in_library"):
                in_library.append((cid, smi))
            else:
                not_in.append((cid, smi))
        except Exception as e:
            print(f"  Error for {cid}: {e}")
            not_in.append((cid, smi))
    return in_library, not_in


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True, help="CSV with id,smiles columns")
    p.add_argument("--mode", choices=["catalog", "api", "manual"], default="catalog")
    p.add_argument("--catalog-file", help="Path to onepot catalog SMILES file (mode=catalog)")
    p.add_argument("--api-url", help="Onepot membership API endpoint (mode=api)")
    p.add_argument("--out", required=True, help="Output CSV with in-library candidates")
    args = p.parse_args()

    candidates = []
    with open(args.candidates) as f:
        for r in csv.DictReader(f):
            if r.get("id") and r.get("smiles"):
                candidates.append((r["id"], r["smiles"]))
    print(f"Candidates to check: {len(candidates)}")

    if args.mode == "catalog":
        if not args.catalog_file:
            print("ERROR: --catalog-file required for mode=catalog"); sys.exit(2)
        in_lib, not_in = filter_against_catalog_file(candidates, args.catalog_file)
    elif args.mode == "api":
        if not args.api_url:
            print("ERROR: --api-url required for mode=api"); sys.exit(2)
        in_lib, not_in = filter_against_api(candidates, args.api_url)
    elif args.mode == "manual":
        # Print SMILES list for manual paste-into-onepot-UI workflow
        print("\nPaste this SMILES list into the Onepot search interface:")
        for cid, smi in candidates:
            print(f"  {cid}\t{smi}")
        print("\nThen create the input CSV manually and re-run with mode=catalog using a hand-built catalog file.")
        sys.exit(0)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "smiles"])
        w.writerows(in_lib)

    print(f"\nResults:")
    print(f"  In onepot library: {len(in_lib)} / {len(candidates)} ({100*len(in_lib)/len(candidates):.1f}%)")
    print(f"  Wrote: {out_path}")
    if in_lib:
        print(f"\n  In-library compounds:")
        for cid, smi in in_lib[:20]:
            print(f"    {cid}: {smi[:60]}")


if __name__ == "__main__":
    main()
