"""SNP closest-drug analysis for genes with nonsense variants.

For each affected gene the gene node (and its edges) is deleted from the full graph, then drugs are
ranked by minimum distance to the gene's MCODE cluster. Reads the GraphML exported by
`drug_target_distance.py` and does not touch Cytoscape.
"""

import networkx as nx
import pandas as pd

from config import (
    CLINVAR_FILTERED_FILE,
    OUTPUT_GRAPHML_FILE,
    SNP_SPECIFIC_TABLE_FILE,
    cluster_name,
)
from utils.data_utils import write_csv_rows
from utils.graph_utils import (
    build_cluster_distance_tables,
    compute_drug_cluster_distances,
    find_mcode_column,
    parse_cluster_list,
)

TOP_N = 10
CONSEQUENCE = "nonsense"
HEADER = [
    "gene",
    "cluster",
    "rank",
    "drug",
    "min_distance",
    "mean_distance",
    "max_distance",
]


def load_affected_genes() -> list:
    """Sorted unique genes carrying at least one nonsense variant."""
    df = pd.read_csv(CLINVAR_FILTERED_FILE, sep="\t", dtype=str, index_col=False)
    is_nonsense = df["Molecular consequence"].fillna("").str.split("|").apply(
        lambda consequences: CONSEQUENCE in consequences
    )
    genes = df.loc[is_nonsense, "Gene(s)"].str.split("|").explode().str.strip()
    return sorted(genes.unique())


def find_gene_nodes(graph: nx.Graph) -> dict:
    """Map gene symbol (display name) to graph node ids."""
    nodes: dict = {}
    for node_id, data in graph.nodes(data=True):
        symbol = data.get("display name")
        if symbol:
            nodes.setdefault(str(symbol), []).append(node_id)
    return nodes


def closest_drugs_after_removal(graph: nx.Graph, node_id: str, mcode_col: str) -> list:
    """For each cluster of the node: delete the node, then rank drugs by distance to the cluster.

    Drugs with identical min/mean/max distances share a rank (1, 1, 3, ...); at most TOP_N drugs are shown.
    """
    clusters = parse_cluster_list(graph.nodes[node_id].get(mcode_col))
    if not clusters:
        return [{"cluster": "NA"}]

    reduced = graph.copy()
    reduced.remove_node(node_id)  # also removes every edge touching it
    tables = build_cluster_distance_tables(compute_drug_cluster_distances(reduced))

    rows = []
    for cluster in clusters:
        ranked = [r for r in tables.get(cluster, []) if r["min_distance"] != "NA"][:TOP_N]
        if not ranked:
            rows.append({"cluster": cluster_name(cluster)})
        rank, previous_key = 0, None
        for position, r in enumerate(ranked, start=1):
            key = (r["min_distance"], r["mean_distance"], r["max_distance"])
            if key != previous_key:  # drugs with identical distances share a rank
                rank, previous_key = position, key
            rows.append({
                "cluster": cluster_name(cluster),
                "rank": rank,
                "drug": r["drug"],
                "min_distance": r["min_distance"],
                "mean_distance": r["mean_distance"],
                "max_distance": r["max_distance"],
            })
    return rows


def main() -> None:
    """Top-10 closest drugs per gene with nonsense variants after deleting the gene node."""
    print("=== SNP Drug Distance Analysis ===")
    graph = nx.read_graphml(str(OUTPUT_GRAPHML_FILE))
    mcode_col = find_mcode_column(graph)
    gene_nodes = find_gene_nodes(graph)

    genes = load_affected_genes()
    print(f"Affected genes (nonsense variants): {len(genes)}")

    rows = []
    missing = []
    for gene in genes:
        node_ids = gene_nodes.get(gene, [])
        if not node_ids:
            missing.append(gene)
        results = [
            r for n in node_ids for r in closest_drugs_after_removal(graph, n, mcode_col)
        ] or [{"cluster": "NA"}]
        for result in results:
            rows.append({"gene": gene, **{col: "NA" for col in HEADER[1:]}, **result})

    if missing:
        print(f"Genes not in the network (reported as NA): {', '.join(missing)}")

    write_csv_rows(SNP_SPECIFIC_TABLE_FILE, HEADER, rows)
    print(f"Exported closest drugs ({len(rows)} rows) -> {SNP_SPECIFIC_TABLE_FILE}")


if __name__ == "__main__":
    main()
