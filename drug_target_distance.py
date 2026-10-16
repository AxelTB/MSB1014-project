"""Step 3: Calculate drug-cluster shortest paths, export distance tables, and assign drug borders.

MUST be executed before `set_graph_style.py` so that `border_color` and `border_width`
columns are populated in Cytoscape before visual style mappings are applied.
"""

import statistics
from pathlib import Path
from typing import Any, Dict, List

import networkx as nx
import pandas as pd
import py4cytoscape as p4c

from config import (
    BORDER_COLOR_COLUMN,
    BORDER_WIDTH_COLUMN,
    CLUSTER_TABLES_DIR,
    CLUSTER_TOP_DRUGS_FILE,
    DEFAULT_BORDER_WIDTH,
    HIGHLIGHTED_BORDER_WIDTH,
    OUTPUT_GRAPHML_FILE,
    cluster_color,
    cluster_name,
)
from utils.cytoscape_utils import export_network_graphml, remove_legacy_path_edges
from utils.data_utils import write_csv_rows
from utils.graph_utils import (
    assign_drug_border_clusters,
    build_cluster_distance_tables,
    classify_nodes,
    compute_drug_cluster_distances,
    get_cluster_top_drugs,
    nodes_by_cluster,
    proteins_by_cluster,
)


def print_graph_statistics(graph: nx.Graph) -> None:
    """Print high-level topological statistics of the imported graph."""
    undirected = graph.to_undirected()
    node_count = graph.number_of_nodes()
    edge_count = graph.number_of_edges()

    print("\n=== Graph Statistics ===")
    print(f"Directed:           {graph.is_directed()}")
    print(f"Nodes:              {node_count}")
    print(f"Edges:              {edge_count}")
    print(f"Self-loops:         {nx.number_of_selfloops(graph)}")
    print(f"Isolated nodes:     {nx.number_of_isolates(graph)}")
    print(f"Density:            {nx.density(graph):.6f}")

    if node_count > 0:
        degrees = [d for _, d in undirected.degree()]
        print(f"Degree min/max:     {min(degrees)} / {max(degrees)}")
        print(f"Degree mean/med:    {statistics.mean(degrees):.2f} / {statistics.median(degrees)}")
        comps = sorted(nx.connected_components(undirected), key=len, reverse=True)
        print(f"Connected comps:    {len(comps)} (largest has {len(comps[0])} nodes)")
        print(f"Avg clustering:     {nx.average_clustering(nx.Graph(undirected)):.4f}")


def verify_network_against_cytoscape(graph: nx.Graph) -> None:
    """Compare imported NetworkX graph counts with the active Cytoscape network."""
    cy_nodes = p4c.get_node_count()
    cy_edges = p4c.get_edge_count()
    nx_nodes = graph.number_of_nodes()
    nx_edges = graph.number_of_edges()

    print("\n=== Verification vs. Cytoscape ===")
    print(f"Nodes: GraphML={nx_nodes} | Cytoscape={cy_nodes} -> {'OK' if nx_nodes == cy_nodes else 'MISMATCH'}")
    print(f"Edges: GraphML={nx_edges} | Cytoscape={cy_edges} -> {'OK' if nx_edges == cy_edges else 'MISMATCH'}")


def get_drug_categories(graph: nx.Graph, drugs: List[str]) -> Dict[str, str]:
    """Map drug name -> DrugBank categories; fall back to the Cytoscape node table when missing."""
    categories: Dict[str, str] = {}
    for _, data in graph.nodes(data=True):
        name = str(data.get("name") or data.get("shared name") or "").strip()
        cats = str(data.get("CTL.categories") or "").strip()
        if name and cats:
            categories[name] = cats

    missing = [d for d in drugs if d not in categories]
    if missing:
        try:
            table = p4c.get_table_columns("node", columns=["name", "CTL.categories"])
            for name, cats in zip(table["name"], table["CTL.categories"]):
                if isinstance(cats, str) and cats.strip():
                    categories.setdefault(str(name), cats.strip())
        except Exception as exc:
            print(f"  Could not retrieve categories from Cytoscape: {exc}")
    return categories


def export_cluster_tables(
    tables: Dict[str, List[Dict[str, Any]]],
    graph: nx.Graph,
    output_file: Path = CLUSTER_TOP_DRUGS_FILE,
    top_n: int = 5,
) -> Path:
    """Write a single CSV with the top N closest drugs of every cluster, with categories and distances."""
    top = {
        cluster: [r for r in rows if r["min_distance"] != "NA"][:top_n]
        for cluster, rows in tables.items()
    }
    categories = get_drug_categories(graph, [r["drug"] for rows in top.values() for r in rows])

    fieldnames = ["cluster", "drug", "category", "min_distance", "max_distance", "mean_distance"]
    out_rows = [
        {
            "cluster": cluster_name(cluster),
            "drug": row["drug"],
            "category": categories.get(row["drug"], "Uncategorized"),
            "min_distance": row["min_distance"],
            "max_distance": row["max_distance"],
            "mean_distance": row["mean_distance"],
        }
        for cluster, rows in top.items()
        for row in rows
    ]
    output_file.parent.mkdir(parents=True, exist_ok=True)
    write_csv_rows(output_file, fieldnames, out_rows)

    print(f"\n=== Top {top_n} Drugs per Cluster ===")
    print(pd.DataFrame(out_rows, columns=fieldnames).to_string(index=False))
    print(f"Table exported to: {output_file}")
    return output_file


def add_drug_border_colors_in_cytoscape(
    graph: nx.Graph, tables: Dict[str, List[Dict[str, Any]]]
) -> int:
    """Set border widths and assigned cluster border colors on drug nodes in Cytoscape."""
    assignments = assign_drug_border_clusters(tables)
    node_table = p4c.get_table_columns("node", columns=["name"])
    suids_by_name = {name: suid for suid, name in node_table["name"].items()}

    cluster_color_map = {c: cluster_color(c) for c in nodes_by_cluster(graph) if c is not None}

    table_rows = [
        {
            "SUID": suid,
            BORDER_COLOR_COLUMN: None,
            BORDER_WIDTH_COLUMN: DEFAULT_BORDER_WIDTH,
        }
        for suid in node_table.index
    ]
    rows_by_suid = {row["SUID"]: row for row in table_rows}
    colored = 0

    for drug, assignment in assignments.items():
        suid = suids_by_name.get(drug)
        if suid is None:
            continue
        assigned_cluster = assignment["cluster"]
        color = cluster_color_map.get(assigned_cluster)
        if color is None:
            continue
        rows_by_suid[suid][BORDER_COLOR_COLUMN] = color
        rows_by_suid[suid][BORDER_WIDTH_COLUMN] = HIGHLIGHTED_BORDER_WIDTH
        colored += 1

    counts = {str(c): 0 for c in tables}
    for assignment in assignments.values():
        counts[str(assignment["cluster"])] = counts.get(str(assignment["cluster"]), 0) + 1
    summary = pd.DataFrame(
        [{"cluster": c, "assigned_drugs": n} for c, n in counts.items()]
    )
    summary_path = CLUSTER_TABLES_DIR.parent / "cluster_assigned_drug_counts.csv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)
    print("\n=== Drugs Assigned per Cluster ===")
    print(summary.to_string(index=False))
    print(f"Total drugs assigned to a cluster: {len(assignments)}")
    print(f"Table exported to: {summary_path}")

    p4c.load_table_data(
        pd.DataFrame(table_rows),
        data_key_column="SUID",
        table="node",
        table_key_column="SUID",
    )
    print(
        f"\n=== Drug Border Styling Loaded in Cytoscape ===\n"
        f"{colored} drugs assigned border width {HIGHLIGHTED_BORDER_WIDTH} & cluster colors; "
        f"all other nodes assigned default width {DEFAULT_BORDER_WIDTH}."
    )
    return colored


def main() -> None:
    """Execute Step 3: Export GraphML, calculate distances, write tables, and assign borders."""
    print("=== Step 3: Drug Target Topological Distance Calculation ===")
    net_name = p4c.get_network_name()
    print(f"Active Cytoscape network: {net_name!r} (SUID {p4c.get_network_suid()})")

    # Clean legacy path edges from earlier versions
    remove_legacy_path_edges()

    # Export live network to GraphML
    export_network_graphml(OUTPUT_GRAPHML_FILE)
    print(f"Exported live network to: {OUTPUT_GRAPHML_FILE}")

    # Load into NetworkX
    graph = nx.read_graphml(str(OUTPUT_GRAPHML_FILE))
    print_graph_statistics(graph)
    verify_network_against_cytoscape(graph)

    # Classify nodes
    snp_genes, interactors, drugs = classify_nodes(graph)
    print(f"\nNode classification: {len(snp_genes)} SNPs/genes, {len(interactors)} interactors, {len(drugs)} drugs.")

    # Calculate shortest path distances
    distances = compute_drug_cluster_distances(graph)
    tables = build_cluster_distance_tables(distances)

    # Export cluster tables to CSV
    export_cluster_tables(tables, graph)

    # Assign drug borders in Cytoscape (MUST run before set_graph_style.py)
    add_drug_border_colors_in_cytoscape(graph, tables)
    print("Drug distance calculation and border attribute assignment complete.")


if __name__ == "__main__":
    main()
