"""Pathway and Gene Ontology (GO) Over-Representation Analysis for MCODE Clusters.

Retrieves the protein gene clusters identified in Cytoscape (or exported GraphML),
performs Gene Ontology (Biological Process, Cellular Component, Molecular Function)
over-representation analysis for each cluster, prints the results to characterize
every cluster, and exports enrichment tables to CSV.
"""

import argparse
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import pandas as pd
import requests

from config import (
    CLUSTER_CHARACTERIZATION_FILE,
    ENRICHMENT_DIR,
    ENRICHMENT_PLOTS_DIR,
    EXTENDED_NETWORK_TITLE,
    NETWORK_TITLE,
    OUTPUT_GRAPHML_FILE,
    SPECIES,
    cluster_color,
)
from utils.cytoscape_utils import check_cytoscape
from utils.graph_utils import nodes_by_cluster


STRING_API_ENRICHMENT_URL = "https://version-12-0.string-db.org/api/json/enrichment"
GO_CATEGORIES = {
    "Process": "GO Biological Process",
    "Component": "GO Cellular Component",
    "Function": "GO Molecular Function",
    "GO Biological Process": "GO Biological Process",
    "GO Cellular Component": "GO Cellular Component",
    "GO Molecular Function": "GO Molecular Function",
}


def get_clusters_from_cytoscape(network_title: Optional[str] = None) -> Dict[str, List[str]]:
    """Retrieve protein gene clusters directly from the active Cytoscape network table."""
    import py4cytoscape as p4c

    existing = p4c.get_network_list()
    if not existing:
        raise RuntimeError("No networks found in Cytoscape.")

    target_network = None
    if network_title and network_title in existing:
        target_network = network_title
    else:
        for preferred in [EXTENDED_NETWORK_TITLE, NETWORK_TITLE]:
            if preferred in existing:
                target_network = preferred
                break
        if not target_network:
            target_network = existing[0]

    network_suid = p4c.get_network_suid(target_network)
    rows = p4c.commands.cyrest_get(f"networks/{network_suid}/tables/defaultnode/rows")

    clusters: Dict[str, List[str]] = {}
    for r in rows:
        ctl_type = str(r.get("CTL.Type", "")).strip().lower()
        if ctl_type == "drug":
            continue

        c_name = r.get("MCODE::Cluster Name (1)")
        if not c_name:
            c_val = r.get("MCODE::Clusters (1)")
            if isinstance(c_val, list) and c_val:
                c_name = str(c_val[0])
            elif isinstance(c_val, str) and c_val.strip() and c_val.strip() not in ("nan", "none", "[]"):
                c_name = c_val.strip().strip("[]").strip("'\"")

        if not c_name or str(c_name).strip().lower() in ("", "nan", "none", "null"):
            continue

        cluster_label = str(c_name).strip()
        gene = str(r.get("display name") or r.get("name") or "").strip()
        if gene:
            clusters.setdefault(cluster_label, []).append(gene)

    # Deduplicate and sort
    return {
        c: sorted(list(dict.fromkeys(genes)))
        for c, genes in sorted(clusters.items())
    }


def get_clusters_from_graphml(graphml_path: Path = OUTPUT_GRAPHML_FILE) -> Dict[str, List[str]]:
    """Retrieve protein gene clusters from an exported GraphML file."""
    if not graphml_path.is_file():
        raise FileNotFoundError(f"GraphML file not found: {graphml_path}")

    graph = nx.read_graphml(str(graphml_path))
    clusters_raw = nodes_by_cluster(graph)

    clusters: Dict[str, List[str]] = {}
    for cluster_id, node_ids in clusters_raw.items():
        if cluster_id is None:
            continue
        cluster_name = f"Cluster {cluster_id}" if not str(cluster_id).startswith("Cluster") else str(cluster_id)
        genes = []
        for nid in node_ids:
            node_data = graph.nodes[nid]
            gene = node_data.get("display name") or node_data.get("name") or nid
            genes.append(str(gene))
        clusters[cluster_name] = sorted(list(dict.fromkeys(genes)))

    return dict(sorted(clusters.items()))


def retrieve_cytoscape_string_enrichment(
    network_title: Optional[str] = None,
    min_group_size: int = 3,
) -> Dict[str, pd.DataFrame]:
    """Retrieve cluster enrichment tables using Cytoscape's stringApp group-wise enrichment."""
    import py4cytoscape as p4c

    existing = p4c.get_network_list()
    target_network = None
    for preferred in [EXTENDED_NETWORK_TITLE, NETWORK_TITLE]:
        if preferred in existing:
            target_network = preferred
            break
    if not target_network:
        target_network = existing[0]

    p4c.set_current_network(target_network)
    net_suid = p4c.get_network_suid(target_network)

    cmd = f'string retrieve group-wise enrichment groupColumn="MCODE::Cluster Name (1)" minGroupSize={min_group_size}'
    res = p4c.commands.commands_post(cmd)

    if not isinstance(res, dict) or not res:
        return {}

    tables = p4c.commands.cyrest_get(f"networks/{net_suid}/tables")
    enrichment_dfs: Dict[str, pd.DataFrame] = {}

    for title, table_suid in res.items():
        # title like 'STRING Enrichment: MCODE::Cluster Name (1) Cluster 1'
        cluster_key = title.split()[-2] + " " + title.split()[-1] if "Cluster" in title else title
        t_info = next((t for t in tables if t.get("SUID") == table_suid), None)
        if not t_info:
            continue

        rows = t_info.get("rows", [])
        if not rows:
            enrichment_dfs[cluster_key] = pd.DataFrame()
            continue

        df = pd.DataFrame(rows)
        # Standardize column names
        rename_map = {
            "term name": "term_id",
            "description": "term_description",
            "p-value": "p_value",
            "FDR value": "fdr",
            "# genes": "number_of_genes",
            "# background genes": "number_of_genes_in_background",
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
        enrichment_dfs[cluster_key] = df

    return enrichment_dfs


def run_string_api_enrichment(
    genes: List[str],
    species_taxon_id: int = 9606,
    caller_identity: str = "network_biology_pipeline",
) -> pd.DataFrame:
    """Run functional over-representation analysis via STRING DB REST API."""
    params = {
        "identifiers": "\n".join(genes),
        "species": species_taxon_id,
        "caller_identity": caller_identity,
    }
    response = requests.post(STRING_API_ENRICHMENT_URL, data=params, timeout=30)
    response.raise_for_status()

    data = response.json()
    if not data or not isinstance(data, list):
        return pd.DataFrame()

    df = pd.DataFrame(data)
    rename_map = {
        "term": "term_id",
        "description": "term_description",
        "preferredNames": "genes",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
    return df


def filter_go_enrichment(
    enrichment_df: pd.DataFrame,
    fdr_threshold: float = 0.05,
    top_n: Optional[int] = 10,
) -> Dict[str, pd.DataFrame]:
    """Filter enrichment results for GO terms (BP, CC, MF) below FDR threshold."""
    if enrichment_df.empty or "category" not in enrichment_df.columns:
        return {
            "GO Biological Process": pd.DataFrame(),
            "GO Cellular Component": pd.DataFrame(),
            "GO Molecular Function": pd.DataFrame(),
        }

    go_results: Dict[str, pd.DataFrame] = {}
    for short_cat, full_cat in [
        ("Process", "GO Biological Process"),
        ("Component", "GO Cellular Component"),
        ("Function", "GO Molecular Function"),
    ]:
        mask = (
            enrichment_df["category"].str.contains(short_cat, case=False, na=False)
            | enrichment_df["category"].str.contains(full_cat, case=False, na=False)
        )
        sub_df = enrichment_df[mask].copy()

        if "fdr" in sub_df.columns:
            sub_df = sub_df[sub_df["fdr"] <= fdr_threshold].sort_values("fdr")
        elif "p_value" in sub_df.columns:
            sub_df = sub_df[sub_df["p_value"] <= fdr_threshold].sort_values("p_value")

        if top_n is not None and len(sub_df) > top_n:
            sub_df = sub_df.head(top_n)

        go_results[full_cat] = sub_df

    return go_results


def print_cluster_characterization(
    cluster_name: str,
    genes: List[str],
    go_results: Dict[str, pd.DataFrame],
    other_pathways: Optional[pd.DataFrame] = None,
) -> None:
    """Print a formatted over-representation analysis report characterizing the cluster."""
    print("=" * 78)
    print(f"  {cluster_name.upper()}: OVER-REPRESENTATION ANALYSIS & FUNCTIONAL CHARACTERIZATION")
    print("=" * 78)
    print(f"Proteins/Genes ({len(genes)}): {', '.join(genes)}")
    print("-" * 78)

    has_go = False
    for category, df in go_results.items():
        if df.empty:
            continue
        has_go = True
        print(f"\n[{category}] (Top significant terms):")
        print(f"{'Term ID':<15} {'FDR':<10} {'Genes':<6} {'Description'}")
        print("-" * 75)
        for _, row in df.iterrows():
            term_id = str(row.get("term_id", ""))
            fdr = row.get("fdr", row.get("p_value", 0.0))
            n_genes = row.get("number_of_genes", len(row.get("genes", [])) if isinstance(row.get("genes"), list) else "")
            desc = str(row.get("term_description", row.get("description", "")))
            if len(desc) > 42:
                desc = desc[:39] + "..."
            gene_list = row.get("genes", [])
            gene_str = f" ({', '.join(gene_list)})" if isinstance(gene_list, list) and len(gene_list) <= 5 else ""

            fdr_str = f"{fdr:.2e}" if isinstance(fdr, (int, float)) else str(fdr)
            print(f"{term_id:<15} {fdr_str:<10} {str(n_genes):<6} {desc}{gene_str}")

    if not has_go:
        print("\nNo GO terms met the significance threshold.")

    if other_pathways is not None and not other_pathways.empty:
        print("\n[Significant Pathways / Disease Annotations]:")
        for _, row in other_pathways.head(5).iterrows():
            cat = row.get("category", "")
            desc = row.get("term_description", row.get("description", ""))
            fdr = row.get("fdr", row.get("p_value", 0.0))
            fdr_str = f"{fdr:.2e}" if isinstance(fdr, (int, float)) else str(fdr)
            print(f"  * [{cat}] {desc} (FDR: {fdr_str})")

    print()


def build_characterization_summary(
    cluster_name: str,
    genes: List[str],
    go_results: Dict[str, pd.DataFrame],
    other_pathways: Optional[pd.DataFrame] = None,
) -> Dict[str, Any]:
    """Build a single-row summary dictionary characterizing a cluster."""
    # Find the top overall term by significance (lowest FDR)
    best_term = None
    best_fdr = 1.0

    for cat in ["GO Biological Process", "GO Cellular Component", "GO Molecular Function"]:
        df = go_results.get(cat)
        if df is not None and not df.empty:
            top_row = df.iloc[0]
            fdr_val = top_row.get("fdr", top_row.get("p_value", 1.0))
            if fdr_val < best_fdr:
                best_fdr = fdr_val
                best_term = (
                    f"{top_row.get('term_description', top_row.get('description', ''))} "
                    f"({top_row.get('term_id', '')}, FDR={best_fdr:.2e})"
                )

    if best_term is None and other_pathways is not None and not other_pathways.empty:
        top_row = other_pathways.iloc[0]
        fdr_val = top_row.get("fdr", top_row.get("p_value", 1.0))
        best_term = (
            f"[{top_row.get('category', '')}] {top_row.get('term_description', top_row.get('description', ''))} "
            f"(FDR={fdr_val:.2e})"
        )

    def _get_top_str(df: Optional[pd.DataFrame]) -> str:
        if df is None or df.empty:
            return "N/A"
        row = df.iloc[0]
        fdr = row.get("fdr", row.get("p_value", 1.0))
        fdr_str = f"{fdr:.2e}" if isinstance(fdr, (int, float)) else str(fdr)
        return f"{row.get('term_description', row.get('description', ''))} ({row.get('term_id', '')}, FDR={fdr_str})"

    return {
        "cluster": cluster_name,
        "size": len(genes),
        "genes": ", ".join(genes),
        "primary_characterization": best_term or "No significant terms",
        "top_go_biological_process": _get_top_str(go_results.get("GO Biological Process")),
        "top_go_cellular_component": _get_top_str(go_results.get("GO Cellular Component")),
        "top_go_molecular_function": _get_top_str(go_results.get("GO Molecular Function")),
        "top_pathway_or_disease": _get_top_str(other_pathways),
    }


def plot_cluster_pathways(
    cluster_name: str,
    enrichment_df: pd.DataFrame,
    total_cluster_genes: int,
    output_dir: Path = ENRICHMENT_PLOTS_DIR,
    fdr_threshold: float = 0.05,
    top_n: int = 15,
) -> Optional[Path]:
    """Generate and save a scatter plot of pathways (x-axis) vs number of cluster genes (y-axis)."""
    if enrichment_df.empty or "category" not in enrichment_df.columns:
        return None

    import matplotlib.pyplot as plt

    # Filter for biological pathway and process categories
    pathway_categories = [
        "KEGG Pathways", "Reactome Pathways", "WikiPathways",
        "GO Biological Process", "KEGG", "RCTM",
    ]
    mask = (
        enrichment_df["category"].isin(pathway_categories)
        | enrichment_df["category"].str.contains("Pathway", case=False, na=False)
        | enrichment_df["category"].str.contains("Process", case=False, na=False)
    )
    df_paths = enrichment_df[mask].copy()

    if df_paths.empty:
        df_paths = enrichment_df.copy()

    # Determine significance column
    sig_col = "fdr" if "fdr" in df_paths.columns else ("p_value" if "p_value" in df_paths.columns else None)
    if sig_col:
        df_paths = df_paths[df_paths[sig_col] <= fdr_threshold]

    if df_paths.empty:
        return None

    # Calculate number of genes
    if "number_of_genes" not in df_paths.columns or df_paths["number_of_genes"].isna().all():
        def _count_genes(val: Any) -> int:
            if isinstance(val, list):
                return len(val)
            if isinstance(val, str):
                cleaned = val.strip("[]").replace("'", "").replace('"', "")
                return len([g for g in cleaned.split(",") if g.strip()])
            return 0
        df_paths["number_of_genes"] = df_paths["genes"].apply(_count_genes)

    # Sort primarily by gene count descending, then by significance (ascending FDR)
    sort_cols = ["number_of_genes"]
    ascending_order = [False]
    if sig_col:
        sort_cols.append(sig_col)
        ascending_order.append(True)
    df_paths = df_paths.sort_values(by=sort_cols, ascending=ascending_order)

    # Select top N terms
    plot_data = df_paths.head(top_n).copy()

    # Labels for x-axis
    name_col = "term_description" if "term_description" in plot_data.columns else (
        "description" if "description" in plot_data.columns else "term_id"
    )

    def _truncate_label(label: Any, max_len: int = 40) -> str:
        s = str(label)
        return s if len(s) <= max_len else s[:max_len - 3] + "..."

    labels = [_truncate_label(x) for x in plot_data[name_col]]
    gene_counts = plot_data["number_of_genes"].tolist()

    # Determine cluster color theme
    point_color = cluster_color(cluster_name)
    if point_color in ("#D9D9D9", "#FFFFFF"):
        point_color = "#1f77b4"

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 0.75 + 2), 6.5))

    x_positions = list(range(len(labels)))
    scatter = ax.scatter(
        x_positions,
        gene_counts,
        s=120,
        c=point_color,
        edgecolors="#222222",
        linewidths=1.5,
        zorder=3,
    )

    # Vertical reference stems connecting the x-axis to each point
    for x_pos, y_val in zip(x_positions, gene_counts):
        ax.plot([x_pos, x_pos], [0, y_val], color="#B0BEC5", linestyle="--", linewidth=1.0, zorder=2)
        # Value label above point
        ax.text(x_pos, y_val + 0.25, str(int(y_val)), ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9.5)
    ax.set_xlabel("Enriched Biological Pathway / Process", fontsize=11, fontweight="bold", labelpad=10)
    ax.set_ylabel(f"Number of Cluster Genes (Total = {total_cluster_genes})", fontsize=11, fontweight="bold")
    ax.set_title(
        f"{cluster_name}: Enriched Pathways vs. Number of Associated Genes",
        fontsize=13,
        fontweight="bold",
        pad=15,
    )

    ax.set_ylim(0, max(gene_counts) + 1.5)
    ax.yaxis.get_major_locator().set_params(integer=True)
    ax.grid(axis="y", linestyle=":", alpha=0.6, zorder=1)

    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    safe_name = cluster_name.replace(" ", "_").lower()
    plot_file = output_dir / f"{safe_name}_pathway_gene_count.png"
    fig.savefig(plot_file, dpi=300)
    plt.close(fig)

    return plot_file


def run_pathway_analysis(
    fdr_cutoff: float = 0.05,
    output_dir: Path = ENRICHMENT_DIR,
    plots_dir: Path = ENRICHMENT_PLOTS_DIR,
    summary_file: Path = CLUSTER_CHARACTERIZATION_FILE,
    export_csv: bool = True,
    export_plots: bool = True,
) -> Tuple[Dict[str, Dict[str, pd.DataFrame]], pd.DataFrame]:
    """Retrieve clusters, execute GO over-representation analysis, print summary, and export CSVs."""
    print("=================================================================")
    print("      PATHWAY & GENE ONTOLOGY (GO) OVER-REPRESENTATION ANALYSIS  ")
    print("=================================================================\n")

    # 1. Retrieve clusters from Cytoscape if online, else fallback to GraphML
    clusters: Dict[str, List[str]] = {}
    use_cytoscape = False
    try:
        check_cytoscape(delete_existing_networks=False)
        clusters = get_clusters_from_cytoscape()
        use_cytoscape = True
        print(f"Successfully retrieved {len(clusters)} clusters directly from active Cytoscape network.\n")
    except Exception as exc:
        print(f"Notice: Cytoscape retrieval not available ({exc}). Falling back to exported GraphML...")
        clusters = get_clusters_from_graphml()
        print(f"Retrieved {len(clusters)} clusters from GraphML file.\n")

    if not clusters:
        sys.exit("ERROR: No clusters found to analyze.")

    for c, g in clusters.items():
        print(f"  - {c} ({len(g)} genes): {', '.join(g)}")
    print()

    # 2. Perform enrichment analysis for each cluster
    all_cluster_results: Dict[str, Dict[str, pd.DataFrame]] = {}
    summary_rows: List[Dict[str, Any]] = []
    if export_csv:
        output_dir.mkdir(parents=True, exist_ok=True)

    # Try Cytoscape stringApp group enrichment if available
    cy_enrichment: Dict[str, pd.DataFrame] = {}
    if use_cytoscape:
        try:
            print("Querying stringApp group-wise enrichment in Cytoscape...")
            cy_enrichment = retrieve_cytoscape_string_enrichment(min_group_size=min(len(g) for g in clusters.values()))
            print(f"Retrieved Cytoscape enrichment tables for {len(cy_enrichment)} cluster(s).\n")
        except Exception as exc:
            print(f"Cytoscape group enrichment failed ({exc}). Using STRING online API fallback.\n")

    for cluster_name, genes in clusters.items():
        enrichment_df = cy_enrichment.get(cluster_name)
        if enrichment_df is None or enrichment_df.empty:
            # Run via STRING DB API
            try:
                enrichment_df = run_string_api_enrichment(genes)
            except Exception as exc:
                print(f"Enrichment query failed for {cluster_name}: {exc}")
                enrichment_df = pd.DataFrame()

        go_results = filter_go_enrichment(enrichment_df, fdr_threshold=fdr_cutoff)

        # Other prominent biological pathways (KEGG, Reactome, WikiPathways, DISEASES)
        other_pathways = pd.DataFrame()
        if not enrichment_df.empty and "category" in enrichment_df.columns:
            pathway_mask = (
                enrichment_df["category"].isin(["KEGG", "KEGG Pathways", "RCTM", "Reactome Pathways", "WikiPathways", "DISEASES"])
                & (enrichment_df.get("fdr", 1.0) <= fdr_cutoff)
            )
            other_pathways = enrichment_df[pathway_mask].sort_values("fdr" if "fdr" in enrichment_df.columns else "p_value")

        # Print detailed report for characterization
        print_cluster_characterization(
            cluster_name=cluster_name,
            genes=genes,
            go_results=go_results,
            other_pathways=other_pathways,
        )

        all_cluster_results[cluster_name] = {
            "all_enrichment": enrichment_df,
            **go_results,
        }

        # Build characterization summary row
        summary_row = build_characterization_summary(
            cluster_name=cluster_name,
            genes=genes,
            go_results=go_results,
            other_pathways=other_pathways,
        )
        summary_rows.append(summary_row)

        # Export CSVs
        if export_csv and not enrichment_df.empty:
            safe_cluster_name = cluster_name.replace(" ", "_").lower()
            csv_path = output_dir / f"{safe_cluster_name}_enrichment.csv"
            enrichment_df.to_csv(csv_path, index=False)

            for cat_name, go_df in go_results.items():
                if not go_df.empty:
                    cat_slug = cat_name.lower().replace(" ", "_")
                    go_path = output_dir / f"{safe_cluster_name}_{cat_slug}.csv"
                    go_df.to_csv(go_path, index=False)

        # Plot pathway scatterplot
        if export_plots and not enrichment_df.empty:
            plot_file = plot_cluster_pathways(
                cluster_name=cluster_name,
                enrichment_df=enrichment_df,
                total_cluster_genes=len(genes),
                output_dir=plots_dir,
                fdr_threshold=fdr_cutoff,
            )
            if plot_file:
                print(f"  Generated pathway scatter plot -> {plot_file}")

    summary_df = pd.DataFrame(summary_rows)

    if export_csv:
        summary_file = Path(summary_file)
        summary_file.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(summary_file, index=False)
        print("=" * 78)
        print("         CLUSTER OVER-REPRESENTATION CHARACTERIZATION SUMMARY")
        print("=" * 78)
        print(f"Exported characterization summary table to: {summary_file}")
        print(f"Detailed enrichment tables saved in:        {output_dir}")
        if export_plots:
            print(f"Pathway gene-count scatter plots saved in: {plots_dir}")
        print()

    return all_cluster_results, summary_df


def main() -> None:
    """CLI entry point for pathway and cluster GO over-representation analysis."""
    parser = argparse.ArgumentParser(
        description="Retrieve MCODE protein clusters from Cytoscape and perform GO over-representation analysis."
    )
    parser.add_argument(
        "--fdr",
        type=float,
        default=0.05,
        help="False Discovery Rate (FDR) significance threshold (default: 0.05).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ENRICHMENT_DIR,
        help="Output directory to save enrichment CSV files.",
    )
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=ENRICHMENT_PLOTS_DIR,
        help="Output directory to save pathway scatter plots.",
    )
    parser.add_argument(
        "--summary-file",
        type=Path,
        default=CLUSTER_CHARACTERIZATION_FILE,
        help="Output CSV path for the high-level cluster characterization summary table.",
    )
    parser.add_argument(
        "--no-export",
        action="store_true",
        help="Do not export CSV tables, only print the report to console.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Do not generate scatter plots.",
    )
    args = parser.parse_args()

    run_pathway_analysis(
        fdr_cutoff=args.fdr,
        output_dir=args.output_dir,
        plots_dir=args.plots_dir,
        summary_file=args.summary_file,
        export_csv=not args.no_export,
        export_plots=not args.no_plots and not args.no_export,
    )


if __name__ == "__main__":
    main()
