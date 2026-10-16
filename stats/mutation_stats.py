"""Mutation stats: total and nonsense mutations per gene from clinvar_variant_filtered."""

import pandas as pd

from config import CLINVAR_FILTERED_FILE, DATA_OUTPUT_DIR

OUT_PATH = DATA_OUTPUT_DIR / "mutation_stats.csv"

df = pd.read_csv(CLINVAR_FILTERED_FILE, sep="\t", dtype=str, index_col=False)
df["is_nonsense"] = df["Molecular consequence"].fillna("").str.split("|").apply(
    lambda c: "nonsense" in c
)
df["gene"] = df["Gene(s)"].fillna("").str.split("|")

table = (
    df.explode("gene")
    .query("gene != ''")
    .groupby("gene")
    .agg(
        total_mutations=("is_nonsense", "size"),
        nonsense_mutations=("is_nonsense", "sum"),
    )
    .astype(int)
    .reset_index()
    .rename(columns={"gene": "affected_gene"})
)
table["nonsense_percentage"] = (
    100 * table["nonsense_mutations"] / table["total_mutations"]
).round(2)
table = table.sort_values(
    ["nonsense_mutations", "total_mutations"], ascending=False
).reset_index(drop=True)

print(f"Total mutations: {len(df)}")
print(f"Nonsense mutations: {int(df['is_nonsense'].sum())}\n")
print(table.to_string(index=False))

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
table.to_csv(OUT_PATH, index=False)
print(f"\nExported to {OUT_PATH}")
