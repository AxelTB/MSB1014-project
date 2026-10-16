"""Graph analysis, shortest path calculations, MCODE cluster parsing, and rankings."""

import statistics
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import networkx as nx

from config import (
    CTL_TYPE_KEYS,
    QUERY_TERM_KEYS,
    TOP_N_DRUGS,
)

MCODE_PREFIX = "MCODE::Clusters"


def get_non_empty_attr(data: Dict[str, Any], keys: Tuple[str, ...]) -> Optional[str]:
    """Return the first non-empty value found in node attribute data for any given key."""
    for key in keys:
        value = data.get(key)
        if value is None:
            continue
        text = str(value).strip()
        if text and text.lower() not in ("nan", "none", "null"):
            return text
    return None


def get_node_label(graph: nx.Graph, node_id: str) -> str:
    """Return human-readable node label ('name', 'shared name', or node id)."""
    data = graph.nodes[node_id]
    return str(data.get("name") or data.get("shared name") or node_id)


def classify_nodes(graph: nx.Graph) -> Tuple[List[str], List[str], List[str]]:
    """Classify graph nodes into (snp_gene_proteins, interacting_proteins, drugs).

    - snp_gene_proteins: nodes with a non-empty query term attribute
    - interacting_proteins: nodes with CTL.Type == 'initial' and no query term
    - drugs: nodes with CTL.Type == 'drug'
    """
    snp_genes, interactors, drugs = [], [], []

    for node_id, data in graph.nodes(data=True):
        label = get_node_label(graph, node_id)
        query_term = get_non_empty_attr(data, QUERY_TERM_KEYS)
        ctl_type = get_non_empty_attr(data, CTL_TYPE_KEYS)
        ctl_type_clean = ctl_type.lower() if ctl_type else None

        if query_term:
            snp_genes.append(label)
        elif ctl_type_clean == "drug":
            drugs.append(label)
        elif ctl_type_clean == "initial":
            interactors.append(label)

    return snp_genes, interactors, drugs


def find_mcode_column(graph: nx.Graph) -> str:
    """Locate the MCODE cluster column name (e.g., 'MCODE::Clusters (1)') in node attributes."""
    all_keys = {k for _, data in graph.nodes(data=True) for k in data.keys()}
    matches = sorted(k for k in all_keys if k.lower().startswith(MCODE_PREFIX.lower()))
    if not matches:
        raise KeyError(
            f"No node attribute starting with {MCODE_PREFIX!r} found. Available attributes: {sorted(all_keys)}"
        )
    return matches[0]


def parse_cluster_list(value: Any) -> List[str]:
    """Parse cluster specification ('3', '[1, 2]', 'Cluster 1') into clean string IDs."""
    if value is None:
        return []
    text = str(value).strip().strip("[]")
    if not text or text.lower() in ("nan", "none", "null"):
        return []
    # Replace semicolons with commas and strip quotations
    cleaned = [p.strip().strip("'\"") for p in text.replace(";", ",").split(",") if p.strip()]
    return cleaned


def nodes_by_cluster(graph: nx.Graph) -> Dict[Optional[str], List[str]]:
    """Group protein nodes (SNP genes and interactors, not drugs) by MCODE cluster ID."""
    mcode_col = find_mcode_column(graph)
    clusters: Dict[Optional[str], List[str]] = {}

    for node_id, data in graph.nodes(data=True):
        query_term = get_non_empty_attr(data, QUERY_TERM_KEYS)
        ctl_type = (get_non_empty_attr(data, CTL_TYPE_KEYS) or "").lower()
        is_protein = bool(query_term) or (ctl_type == "initial")
        if not is_protein:
            continue

        raw_clusters = parse_cluster_list(data.get(mcode_col))
        assigned = raw_clusters if raw_clusters else [None]

        for cluster in assigned:
            clusters.setdefault(cluster, []).append(node_id)

    def sort_key(item: Tuple[Optional[str], List[str]]) -> Tuple[int, float, str]:
        c = item[0]
        if c is None:
            return (2, 0.0, "")
        try:
            return (0, float(c), "")
        except ValueError:
            return (1, 0.0, str(c))

    return dict(sorted(clusters.items(), key=sort_key))


def proteins_by_cluster(graph: nx.Graph) -> Dict[Optional[str], List[str]]:
    """Map MCODE clusters to human-readable protein node labels."""
    return {
        c: [get_node_label(graph, nid) for nid in node_ids]
        for c, node_ids in nodes_by_cluster(graph).items()
    }


def compute_drug_cluster_distances(
    graph: nx.Graph,
) -> Dict[str, Dict[str, Dict[str, Optional[int]]]]:
    """Compute shortest-path distances (edge count, undirected) between each drug and cluster proteins.

    Returns:
        {drug_label: {cluster: {protein_label: distance or None}}}
    """
    undirected = graph.to_undirected()
    clusters = {c: nids for c, nids in nodes_by_cluster(graph).items() if c is not None}

    drug_nodes = [
        node_id
        for node_id, data in graph.nodes(data=True)
        if not get_non_empty_attr(data, QUERY_TERM_KEYS)
        and (get_non_empty_attr(data, CTL_TYPE_KEYS) or "").lower() == "drug"
    ]

    distances: Dict[str, Dict[str, Dict[str, Optional[int]]]] = {}
    for drug_id in drug_nodes:
        drug_label = get_node_label(graph, drug_id)
        lengths = nx.single_source_shortest_path_length(undirected, drug_id)
        distances[drug_label] = {
            c: {get_node_label(graph, p_id): lengths.get(p_id) for p_id in p_ids}
            for c, p_ids in clusters.items()
        }

    return distances


def build_cluster_distance_tables(
    distances: Dict[str, Dict[str, Dict[str, Optional[int]]]]
) -> Dict[str, List[Dict[str, Any]]]:
    """Organize distance dict into ranked table rows per MCODE cluster.

    Each row contains: drug, min_distance, max_distance, mean_distance.
    Connected drugs are sorted by min_distance, then mean, max, drug name.
    Unreachable drugs are listed last with 'NA'.
    """
    connected: Dict[str, List[Dict[str, Any]]] = {}
    unreachable: Dict[str, List[Dict[str, Any]]] = {}

    for drug, per_cluster in distances.items():
        for cluster, per_protein in per_cluster.items():
            connected.setdefault(cluster, [])
            unreachable.setdefault(cluster, [])
            reachable = [d for d in per_protein.values() if d is not None]
            if reachable:
                connected[cluster].append({
                    "drug": drug,
                    "min_distance": min(reachable),
                    "max_distance": max(reachable),
                    "mean_distance": round(statistics.mean(reachable), 2),
                })
            else:
                unreachable[cluster].append({
                    "drug": drug,
                    "min_distance": "NA",
                    "max_distance": "NA",
                    "mean_distance": "NA",
                })

    tables: Dict[str, List[Dict[str, Any]]] = {}
    for cluster in connected:
        connected[cluster].sort(
            key=lambda r: (
                r["min_distance"],
                r["mean_distance"],
                r["max_distance"],
                str(r["drug"]),
            )
        )
        unreachable[cluster].sort(key=lambda r: str(r["drug"]))
        tables[cluster] = connected[cluster] + unreachable[cluster]

    return tables


def top_ranked_rows(candidates: List[Dict[str, Any]], top_n: int = TOP_N_DRUGS) -> List[Dict[str, Any]]:
    """Pick top N ranked rows including all tied rows at the cutoff distance."""
    if not candidates:
        return []
    cutoff = candidates[min(top_n, len(candidates)) - 1]["min_distance"]
    return [r for r in candidates if r["min_distance"] <= cutoff]


def get_cluster_top_drugs(
    tables: Dict[str, List[Dict[str, Any]]], top_n: int = TOP_N_DRUGS
) -> Dict[str, List[Dict[str, Any]]]:
    """Return top N closest drugs (plus ties) for every cluster."""
    return {
        cluster: top_ranked_rows(
            [r for r in rows if r["min_distance"] != "NA"], top_n=top_n
        )
        for cluster, rows in tables.items()
    }


def assign_drug_border_clusters(
    tables: Dict[str, List[Dict[str, Any]]]
) -> Dict[str, Dict[str, Any]]:
    """Assign each reachable drug to the single MCODE cluster where it is ranked highest.

    Resolves ties by minimum distance, then by cluster with fewer assigned drugs.
    """
    sizes: Dict[str, int] = {}
    appearances: Dict[str, List[Tuple[int, int, str]]] = {}

    for cluster, rows in tables.items():
        connected = [r for r in rows if r["min_distance"] != "NA"]
        sizes[cluster] = len(connected)
        for rank, row in enumerate(connected, start=1):
            appearances.setdefault(row["drug"], []).append(
                (rank, row["min_distance"], cluster)
            )

    assignment: Dict[str, Dict[str, Any]] = {}
    skipped = []
    counts = {cluster: 0 for cluster in tables}

    def assign(drug: str, found: List[Tuple[int, int, str]], best: Tuple[int, int, str]) -> None:
        assignment[drug] = {
            "cluster": best[2],
            "rank": best[0],
            "others": [c for _, _, c in found if c != best[2]],
        }
        counts[best[2]] += 1

    for drug, found in sorted(appearances.items(), key=lambda item: str(item[0])):
        best_key = min(item[:2] for item in found)
        tied = [item for item in found if item[:2] == best_key]
        if len(tied) == 1:
            assign(drug, found, tied[0])
        else:
            skipped.append((drug, found, tied))

    for drug, found, tied in skipped:
        best = min(tied, key=lambda item: (counts[item[2]], sizes[item[2]]))
        assign(drug, found, best)

    return assignment


def compute_snp_closest_drugs(
    graph: nx.Graph,
    tables: Dict[str, List[Dict[str, Any]]],
    top_n: int = TOP_N_DRUGS,
    only_node: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Identify closest candidate drugs for each disease-associated SNP protein."""
    mcode_col = find_mcode_column(graph)
    rows: List[Dict[str, Any]] = []

    snps = [
        (nid, data)
        for nid, data in graph.nodes(data=True)
        if get_non_empty_attr(data, QUERY_TERM_KEYS)
        and (only_node is None or nid == only_node)
    ]
    snps.sort(key=lambda item: str(get_node_label(graph, item[0])))

    for nid, data in snps:
        label = get_node_label(graph, nid)
        clusters = parse_cluster_list(data.get(mcode_col))
        if not clusters:
            rows.append({
                "snp_protein": label,
                "cluster": "NA",
                "drug": "NA",
                "cluster_drug_min_distance": "NA",
            })
            continue

        for cluster in clusters:
            candidates = [r for r in tables.get(cluster, []) if r["min_distance"] != "NA"]
            ranked = top_ranked_rows(candidates, top_n)
            if not ranked:
                rows.append({
                    "snp_protein": label,
                    "cluster": cluster,
                    "drug": "NA",
                    "cluster_drug_min_distance": "NA",
                })
                continue

            for r in ranked:
                rows.append({
                    "snp_protein": label,
                    "cluster": cluster,
                    "drug": r["drug"],
                    "cluster_drug_min_distance": r["min_distance"],
                })

    return rows


def get_snp_drug_edges(graph: nx.Graph, snp_node: str) -> List[Tuple[str, str]]:
    """Find all edges directly linking an SNP protein node to a drug node."""
    edges = []
    for u, v in graph.edges():
        if snp_node not in (u, v) or u == v:
            continue
        other = v if u == snp_node else u
        other_data = graph.nodes[other]
        if (
            not get_non_empty_attr(other_data, QUERY_TERM_KEYS)
            and (get_non_empty_attr(other_data, CTL_TYPE_KEYS) or "").lower() == "drug"
        ):
            edges.append((u, v))
    return edges


def compute_edge_removal_snp_closest_drugs(
    graph: nx.Graph,
    distances: Dict[str, Dict[str, Dict[str, Optional[int]]]],
    top_n: int = TOP_N_DRUGS,
) -> List[Dict[str, Any]]:
    """Perform sensitivity analysis: recompute closest drugs with direct SNP-drug edges removed."""
    original_tables = build_cluster_distance_tables(distances)
    snps = sorted(
        [
            nid
            for nid, data in graph.nodes(data=True)
            if get_non_empty_attr(data, QUERY_TERM_KEYS)
        ],
        key=lambda n: str(get_node_label(graph, n)),
    )

    rows: List[Dict[str, Any]] = []
    for snp_node in snps:
        direct_edges = get_snp_drug_edges(graph, snp_node)
        if direct_edges:
            mod_graph = graph.copy()
            mod_graph.remove_edges_from(direct_edges)
            mod_distances = compute_drug_cluster_distances(mod_graph)
            mod_tables = build_cluster_distance_tables(mod_distances)
        else:
            mod_graph = graph
            mod_tables = original_tables

        rows.extend(
            compute_snp_closest_drugs(
                mod_graph, mod_tables, top_n=top_n, only_node=snp_node
            )
        )

    return rows
