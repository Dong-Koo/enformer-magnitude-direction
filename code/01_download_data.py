"""
01_download_data.py
-------------------
Download all raw data files needed for the eQTL / Enformer U-shape analysis.

Downloads:
  1. Enformer per-gene Pearson R correlations (Sasse et al. 2023)
  2. Supplementary Table 1 (Sasse et al. 2023)
  3. GTEx v8 Brain Cortex eGenes via GTEx Portal API v2

Outputs (relative to project root eqtl_enformer_ushape/):
  data/raw/enformer_correlations.txt
  data/raw/supp_table1.tsv
  data/raw/gtex_brain_cortex_egenes.txt
"""

import json
import os
import urllib.request

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
os.makedirs(RAW_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Source URLs
# ---------------------------------------------------------------------------
ENFORMER_URL = (
    "https://raw.githubusercontent.com/mostafavilabuw/EnformerAssessment/"
    "main/Data/Prediction_correlationsCageAdultBrain_Allstats.txt"
)
SUPP1_URL = (
    "https://raw.githubusercontent.com/mostafavilabuw/EnformerAssessment/"
    "main/Data/SupplementaryTable1.tsv"
)
GTEX_API_BASE = (
    "https://gtexportal.org/api/v2/association/egene"
    "?tissueSiteDetailId=Brain_Cortex&datasetId=gtex_v8"
    "&itemsPerPage=2000&page={page}"
)


def download_file(url: str, dest_path: str) -> None:
    """Download *url* to *dest_path*, printing progress."""
    filename = os.path.basename(dest_path)
    print(f"  Downloading {filename} ...", end=" ", flush=True)
    headers = {"User-Agent": "Mozilla/5.0 (eqtl-ushape-analysis)"}
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=120) as response:
        data = response.read()
    with open(dest_path, "wb") as fh:
        fh.write(data)
    kb = len(data) / 1024
    print(f"done ({kb:.1f} KB -> {dest_path})")


def fetch_gtex_egenes(output_path: str) -> int:
    """
    Paginate through GTEx Portal API v2 to collect all Brain Cortex eGene IDs.
    Strips version suffix (.X) from GENCODE IDs.
    Returns the total number of unique eGenes written.
    """
    print("  Fetching GTEx Brain Cortex eGenes via API ...")
    egene_ids: set[str] = set()
    page = 0
    while True:
        url = GTEX_API_BASE.format(page=page)
        print(f"    page {page} -> {url}", flush=True)
        req = urllib.request.Request(
            url, headers={"User-Agent": "Mozilla/5.0 (eqtl-ushape-analysis)"}
        )
        with urllib.request.urlopen(req, timeout=120) as response:
            payload = json.loads(response.read().decode("utf-8"))

        # GTEx v2 API wraps results in {"data": [...], ...}
        data = payload.get("data", [])
        print(f"    received {len(data)} records")
        if not data:
            break

        for record in data:
            raw_id = record.get("gencodeId", "")
            # Strip version suffix: ENSG00000XXXXXX.Y -> ENSG00000XXXXXX
            gene_id = raw_id.split(".")[0] if raw_id else ""
            if gene_id:
                egene_ids.add(gene_id)

        if len(data) < 2000:
            break
        page += 1

    sorted_ids = sorted(egene_ids)
    with open(output_path, "w") as fh:
        fh.write("\n".join(sorted_ids) + "\n")

    print(f"    Written {len(sorted_ids)} unique eGenes -> {output_path}")
    return len(sorted_ids)


def main() -> None:
    print("=" * 60)
    print("01_download_data.py  —  Downloading raw data")
    print("=" * 60)

    # 1. Enformer correlations
    enformer_path = os.path.join(RAW_DIR, "enformer_correlations.txt")
    download_file(ENFORMER_URL, enformer_path)

    # Count non-comment lines
    with open(enformer_path) as fh:
        n_genes = sum(1 for line in fh if line.strip() and not line.startswith("#"))
    print(f"  -> {n_genes} gene records in enformer_correlations.txt")

    # 2. Supplementary Table 1
    supp_path = os.path.join(RAW_DIR, "supp_table1.tsv")
    download_file(SUPP1_URL, supp_path)

    with open(supp_path) as fh:
        n_supp = sum(1 for line in fh if line.strip() and not line.startswith("#"))
    # subtract header row
    n_supp = max(0, n_supp - 1)
    print(f"  -> {n_supp} gene records in supp_table1.tsv")

    # 3. GTEx eGenes
    gtex_path = os.path.join(RAW_DIR, "gtex_brain_cortex_egenes.txt")
    n_egenes = fetch_gtex_egenes(gtex_path)

    print()
    print("Summary")
    print("-------")
    print(f"  Enformer gene records : {n_genes}")
    print(f"  Supp Table 1 records  : {n_supp}")
    print(f"  GTEx Brain Cortex eGenes: {n_egenes}")
    print()
    print("All raw data files are in:", RAW_DIR)
    print("Done.")


if __name__ == "__main__":
    main()
