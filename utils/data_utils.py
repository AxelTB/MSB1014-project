"""Data input/output, parsing, and export helper functions."""

import csv
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Union

import pandas as pd


def load_unique_genes(path: Union[str, Path]) -> List[str]:
    """Load list of unique gene identifiers from a TSV or CSV file.

    Expects a single column (or column named 'id' or 'Gene'), strips whitespace,
    removes empty lines, and returns a deduplicated list of gene symbols.
    """
    path_str = str(path)
    if not os.path.exists(path_str):
        raise FileNotFoundError(f"Gene file not found: {path_str}")

    df = pd.read_csv(path_str, sep="\t", header=0)
    # Pick the first column regardless of its column name
    col_name = df.columns[0]
    genes = df[col_name].dropna().astype(str).str.strip()
    unique_genes = genes[genes != ""].drop_duplicates().tolist()
    return unique_genes


def extract_and_export_unique_genes(
    input_file: Union[str, Path],
    output_file: Union[str, Path],
    gene_col: str = "Gene(s)",
) -> List[str]:
    """Extract unique gene symbols from variant association table and save to TSV.

    Splits multiple genes joined with '|' into individual symbols, deduplicates,
    and exports a single-column table with header 'Gene'.
    """
    df = pd.read_csv(str(input_file), sep="\t", dtype=str)
    if gene_col not in df.columns:
        # Fallback if already named 'Gene'
        if "Gene" in df.columns:
            gene_col = "Gene"
        else:
            raise KeyError(f"Column '{gene_col}' not found in {input_file}. Columns: {df.columns.tolist()}")

    unique_genes = set()
    for value in df[gene_col].dropna():
        for gene in str(value).split("|"):
            cleaned = gene.strip()
            if cleaned:
                unique_genes.add(cleaned)

    sorted_genes = sorted(unique_genes)
    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    out_df = pd.DataFrame(sorted_genes, columns=["Gene"])
    out_df.to_csv(out_path, sep="\t", index=False)
    return sorted_genes


def clean_variant_data(
    input_file: Union[str, Path],
    output_file: Union[str, Path],
) -> pd.DataFrame:
    """Clean ClinVar variant table, explode multi-gene entries, and save to TSV."""
    variant_df_raw = pd.read_csv(str(input_file), sep="\t", index_col=False, low_memory=False)

    # Drop columns that are completely empty
    variant_df_cleaned = variant_df_raw.dropna(axis=1, how="all")

    # Explode multi-gene rows separated by '|'
    if "Gene(s)" in variant_df_cleaned.columns:
        variant_df_cleaned = (
            variant_df_cleaned.assign(
                **{"Gene(s)": variant_df_cleaned["Gene(s)"].str.split("|")}
            )
            .explode("Gene(s)")
            .reset_index(drop=True)
        )
        variant_df_cleaned["Gene"] = variant_df_cleaned["Gene(s)"].astype(str)
        variant_df_cleaned = variant_df_cleaned.drop(columns=["Gene(s)"])

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    variant_df_cleaned.to_csv(out_path, sep="\t", index=False)
    return variant_df_cleaned


def write_csv_rows(
    path: Union[str, Path],
    fieldnames: List[str],
    rows: Iterable[Dict[str, Any]],
) -> Path:
    """Safely write a list of dictionaries to a CSV file."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return out_path
