"""Cluster stats: protein counts and per-cluster MCODE topology metrics from the exported network."""

import networkx as nx
import pandas as pd

from config import DATA_OUTPUT_DIR, OUTPUT_GRAPHML_FILE, cluster_name
from utils.graph_utils import find_mcode_column, nodes_by_cluster

OUT_PATH = DATA_OUTPUT_DIR / "cluster_stats.csv"


def cluster_metrics(graph: nx.Graph, node_ids: list) -> dict:
    """Topology metrics of the subgraph induced by one cluster's proteins."""
    sub = graph.subgraph(node_ids)
    size = sub.number_of_nodes()
    density = nx.density(sub) if size > 1 else 0.0
    # Diameter is computed on the largest connected component if the cluster is split
    core = sub.subgraph(max(nx.connected_components(sub), key=len)) if size else sub
    connected = core.number_of_nodes() > 1
    return {
        "proteins": size,
        # MCODE cluster score = density * number of nodes
        "mcode_score": round(density * size, 3),
        "density": round(density, 3),
        "avg_clustering_coefficient": round(nx.average_clustering(sub), 3) if size else 0.0,
        "diameter": nx.diameter(core) if connected else 0,
    }


def main() -> pd.DataFrame:
    graph = nx.read_graphml(OUTPUT_GRAPHML_FILE).to_undirected()
    graph.remove_edges_from(list(nx.selfloop_edges(graph)))
    find_mcode_column(graph)

    groups = nodes_by_cluster(graph)
    clustered = {c: ids for c, ids in groups.items() if c is not None}
    clustered_count = len({n for ids in clustered.values() for n in ids})
    total = clustered_count + len(groups.get(None, []))

    table = pd.DataFrame(
        [{"cluster": c, "cluster_name": cluster_name(c), **cluster_metrics(graph, ids)}
         for c, ids in clustered.items()]
    )

    print(f"Total proteins: {total}")
    print(f"Clustered proteins: {clustered_count}\n")
    print(table.to_string(index=False))

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT_PATH, index=False)
    print(f"\nExported to {OUT_PATH}")
    return table


if __name__ == "__main__":
    main()
