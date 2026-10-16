"""STRING functional enrichment per MCODE cluster of the created network.

For each cluster in the exported network (GraphML) this script:
  1. runs STRING functional enrichment (REST API),
  2. keeps only GO, KEGG and WikiPathways terms and writes a table with the enriched terms (and the proteins behind each term),
  3. draws a graph with the enriched functions on y, the number of cluster genes
     annotated to each function on x, and every protein identified by name.
"""

import argparse
import textwrap
from pathlib import Path
from typing import List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from config import ENRICHMENT_DIR, OUTPUT_GRAPHML_FILE, SPECIES, cluster_color, cluster_name
from pathway_analysis import get_clusters_from_graphml, run_string_api_enrichment

FDR_THRESHOLD = 0.05
FIG_WIDTH = 16  # inches; identical for all graphs
OUT_DIR = ENRICHMENT_DIR / "string_functional"
ALLOWED_CATEGORIES = {"Process", "Function", "KEGG", "WikiPathways"}  # GO BP/MF, KEGG, WikiPathways
TABLE_COLUMNS = ["category", "term_id", "term_description", "number_of_genes",
                 "number_of_genes_in_background", "p_value", "fdr", "genes"]


def _split_genes(value) -> List[str]:
    if isinstance(value, (list, tuple)):
        return [str(g) for g in value]
    return [g.strip() for g in str(value).split(",") if g.strip()]


def plot_cluster(cluster: str, genes: List[str], df: pd.DataFrame, path: Path, top_n: int) -> None:
    """Bar plot: term (y) vs number of module genes annotated to it (x), labelled with FDR."""
    # most genes first (top of the plot); ties broken by lower FDR
    top = df.sort_values(["number_of_genes", "fdr"], ascending=[False, True]).head(top_n)
    top = top.iloc[::-1]  # barh draws the first row at the bottom
    plt.rcParams["font.size"] = 11  # +10% over the matplotlib default of 10
    fig, ax = plt.subplots(figsize=(FIG_WIDTH, 0.4 * len(top) + 1.5))
    color = cluster_color(cluster)
    labels = [textwrap.shorten(f"{d} [{c}]", 70) for d, c in zip(top["term_description"], top["category"])]
    counts = top["number_of_genes"]
    ax.barh(range(len(top)), counts, height=0.7, color=color, edgecolor="black", zorder=3)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=8.8)
    ax.grid(axis="x", linestyle=":", alpha=0.5, zorder=0)
    for y, (n, fdr) in enumerate(zip(counts, top["fdr"])):
        ax.annotate(f"FDR {fdr:.1e}", (n, y), xytext=(5, 0), textcoords="offset points",
                    va="center", fontsize=7.7)
    xmin = counts.min() - 1
    ax.set_xlim(xmin, len(genes) * 1.15)
    # axis-break marker: the x axis does not start at 0
    ax.plot([0.012, 0.028], [0, 0], transform=ax.transAxes, color="white", lw=3, clip_on=False, zorder=4)
    for x0 in (0.008, 0.024):
        ax.plot([x0, x0 + 0.008], [-0.012, 0.012], transform=ax.transAxes, color="black", lw=1.2,
                clip_on=False, zorder=5)
    ax.set_ylim(-0.7, len(top) - 0.3)
    ax.xaxis.get_major_locator().set_params(integer=True)
    ax.set_xlabel(f"Number of genes in module ({len(genes)} total)")
    ax.set_ylabel("STRING enriched function")
    ax.set_title(f"{cluster_name(cluster)}: STRING functional enrichment")
    height = fig.get_figheight()
    # fixed margins in inches: same plot width for every cluster, only the height follows the row count
    fig.subplots_adjust(left=0.3, right=0.98, top=1 - 0.5 / height, bottom=0.9 / height)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graphml", type=Path, default=OUTPUT_GRAPHML_FILE)
    parser.add_argument("--fdr", type=float, default=FDR_THRESHOLD, help="FDR threshold")
    parser.add_argument("--top", type=int, default=10, help="terms plotted per cluster")
    parser.add_argument("--out", type=Path, default=OUT_DIR)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    clusters = get_clusters_from_graphml(args.graphml)
    taxon = 9606 if SPECIES == "Homo sapiens" else 9606
    summary = []

    for cluster, genes in clusters.items():
        slug = cluster.lower().replace(" ", "_")
        if len(genes) < 2:
            print(f"{cluster}: skipped (fewer than 2 genes)")
            continue
        df = run_string_api_enrichment(genes, species_taxon_id=taxon)
        if df.empty or "fdr" not in df.columns:
            print(f"{cluster}: no enrichment returned")
            continue
        df = df[df["category"].isin(ALLOWED_CATEGORIES) & (df["fdr"] <= args.fdr)].sort_values("fdr")
        plot_df = df
        cols = [c for c in TABLE_COLUMNS if c in df.columns]
        df[cols].to_csv(args.out / f"{slug}_string_enrichment.csv", index=False)
        print(f"{cluster_name(cluster)}: {len(genes)} genes, {len(df)} significant terms")
        summary.append({"cluster": cluster_name(cluster), "n_genes": len(genes), "n_terms": len(df)})
        if not plot_df.empty:
            plot_cluster(cluster, genes, plot_df, args.out / f"{slug}_string_enrichment.png", args.top)

    pd.DataFrame(summary).to_csv(args.out / "summary.csv", index=False)
    print(f"Results written to {args.out}")


if __name__ == "__main__":
    main()
