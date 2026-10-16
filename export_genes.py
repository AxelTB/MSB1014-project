"""Step 1: Extract unique genes from ClinVar variant associations and export to TSV."""

import sys
from pathlib import Path

from config import UNIQUE_GENES_FILE, VARIANT_GENE_NETWORK_FILE
from utils.data_utils import extract_and_export_unique_genes


def main() -> None:
    """Extract unique genes from variant data and write to data/input/unique_genes.tsv."""
    print("=== Step 1: Extracting Unique Genes ===")
    print(f"Reading variants from: {VARIANT_GENE_NETWORK_FILE}")

    if not VARIANT_GENE_NETWORK_FILE.is_file():
        sys.exit(f"ERROR: Input variant file not found: {VARIANT_GENE_NETWORK_FILE}")

    unique_genes = extract_and_export_unique_genes(
        input_file=VARIANT_GENE_NETWORK_FILE,
        output_file=UNIQUE_GENES_FILE,
        gene_col="Gene(s)",
    )

    print(f"Extracted {len(unique_genes)} unique genes.")
    print(f"Exported gene list to: {UNIQUE_GENES_FILE}")
    print(f"Sample genes: {unique_genes[:5]}")


if __name__ == "__main__":
    main()
