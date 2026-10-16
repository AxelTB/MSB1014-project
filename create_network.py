"""Step 2: Build STRING PPI network, cluster with MCODE, and extend with CyTargetLinker."""

import glob
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import py4cytoscape as p4c

from config import (
    ADDITIONAL_INTERACTORS,
    CTL_DIRECTION,
    CTL_ID_ATTRIBUTE,
    DARK_GRAY,
    DELETE_EXISTING_NETWORKS,
    DRUG_SHAPE,
    DRUG_SIZE,
    GENE_SHAPE,
    LINKSETS_DIR,
    MCODE_FLUFF,
    MCODE_HAIRCUT,
    MCODE_K_CORE,
    MCODE_MAX_DEPTH,
    MCODE_NODE_SCORE_CUTOFF,
    NETWORK_TITLE,
    NODE_COLOR_COLUMN,
    NODE_SHAPE_COLUMN,
    NODE_SIZE_COLUMN,
    OTHER_SIZE,
    SNP_GENE_SIZE,
    SPECIES,
    STRING_CONFIDENCE_THRESHOLD,
    UNIQUE_GENES_FILE,
    cluster_color,
    cluster_number,
)
from utils.cytoscape_utils import check_cytoscape
from utils.data_utils import load_unique_genes


def build_string_network(
    genes: List[str],
    confidence_threshold: float = STRING_CONFIDENCE_THRESHOLD,
    additional_interactors: int = ADDITIONAL_INTERACTORS,
    species: str = SPECIES,
    title: str = NETWORK_TITLE,
) -> int:
    """Query the STRING database directly in Cytoscape for the given genes."""
    query = ",".join(genes)
    p4c.commands.commands_post(
        f'string protein query query="{query}" species="{species}" '
        f'cutoff={confidence_threshold} limit={additional_interactors}'
    )

    network_suid = p4c.get_network_suid()
    if network_suid is None:
        raise RuntimeError("STRING query did not create a network in Cytoscape.")

    p4c.set_current_network(network_suid)
    p4c.rename_network(title, network=network_suid)

    nodes = p4c.get_table_columns("node", columns=["query term"], network=network_suid)
    n_query = nodes["query term"].replace("", pd.NA).notna().sum()
    print(
        f"Network '{title}' (SUID {network_suid}): {len(nodes)} nodes, "
        f"{n_query} of {len(genes)} input genes mapped, "
        f"{len(nodes) - n_query} added interactors."
    )
    return network_suid


def run_mcode(
    network_suid: int,
    node_score_cutoff: float = MCODE_NODE_SCORE_CUTOFF,
    k_core: int = MCODE_K_CORE,
    max_depth: int = MCODE_MAX_DEPTH,
    haircut: bool = MCODE_HAIRCUT,
    fluff: bool = MCODE_FLUFF,
) -> Optional[str]:
    """Cluster the network with the standalone MCODE app and record cluster membership."""
    p4c.set_current_network(network_suid)

    result = p4c.commands.commands_post(
        f"mcode cluster network=SUID:{network_suid} scope=NETWORK includeLoops=false "
        f"nodeScoreCutoff={node_score_cutoff} kCore={k_core} "
        f"maxDepthFromStart={max_depth} haircut={str(haircut).lower()} "
        f"fluff={str(fluff).lower()}"
    )

    clusters = (result or {}).get("clusters", [])
    if not clusters:
        print("MCODE ran, but found no clusters.")
        return None

    result_id = result.get("id", 1)
    col = f"MCODE::Cluster Name ({result_id})"
    membership: Dict[int, List[str]] = {}
    for cluster in clusters:
        for suid in cluster["nodes"]:
            membership.setdefault(suid, []).append(cluster["name"])

    all_suids = p4c.get_table_columns("node", columns=["name"], network=network_suid).index
    data = pd.DataFrame({
        "SUID": list(all_suids),
        col: [", ".join(membership.get(s, [])) for s in all_suids],
    })
    p4c.load_table_data(
        data,
        data_key_column="SUID",
        table="node",
        table_key_column="SUID",
        network=network_suid,
    )

    print(f"MCODE identified {len(clusters)} cluster(s) (recorded in column '{col}'):")
    for cluster in clusters:
        print(f"  {cluster['name']}: {len(cluster['nodes'])} nodes, score {cluster['score']:.2f}")
    return col


def _parse_cluster_ids(value: object) -> List[str]:
    """Parse cluster specification cell into a list of cluster id strings."""
    if isinstance(value, (list, tuple, set)):
        items = list(value)
    else:
        items = str(value).strip("[] ").split(",")
    return [str(i).strip() for i in items if str(i).strip() not in ("", "nan", "None")]


def add_cluster_colors(
    network_suid: int,
    cluster_col: Optional[str],
    color_column: str = NODE_COLOR_COLUMN,
) -> Dict[str, str]:
    """Assign node fill colors, shapes, and sizes based on MCODE cluster and query status."""
    if cluster_col is None:
        nodes = p4c.get_table_columns("node", columns=["name"], network=network_suid)
        first = pd.Series([None] * len(nodes), index=nodes.index)
    else:
        nodes = p4c.get_table_columns("node", columns=[cluster_col], network=network_suid)
        first = nodes[cluster_col].map(lambda v: (_parse_cluster_ids(v) or [None])[0])

    names = sorted({n for n in first.dropna()}, key=lambda n: cluster_number(n) or 0)
    color_of = {n: cluster_color(n) for n in names}

    color_data = pd.DataFrame({
        "SUID": list(nodes.index),
        color_column: [cluster_color(n) for n in first],
    })
    p4c.load_table_data(
        color_data,
        data_key_column="SUID",
        table="node",
        table_key_column="SUID",
        network=network_suid,
    )

    shape_data = pd.DataFrame({
        "SUID": list(nodes.index),
        NODE_SHAPE_COLUMN: [GENE_SHAPE] * len(nodes),
    })
    p4c.load_table_data(
        shape_data,
        data_key_column="SUID",
        table="node",
        table_key_column="SUID",
        network=network_suid,
    )

    queries = p4c.get_table_columns("node", columns=["query term"], network=network_suid)["query term"]
    is_snp = queries.map(lambda v: isinstance(v, str) and v.strip() != "")
    size_data = pd.DataFrame({
        "SUID": list(queries.index),
        NODE_SIZE_COLUMN: [SNP_GENE_SIZE if x else OTHER_SIZE for x in is_snp],
    })
    p4c.load_table_data(
        size_data,
        data_key_column="SUID",
        table="node",
        table_key_column="SUID",
        network=network_suid,
    )

    print(f"Added visual attributes to '{NETWORK_TITLE}' (SUID {network_suid}):")
    for n, c in color_of.items():
        print(f"  {n}: {c}")
    return color_of


def extend_with_cytargetlinker(
    network_suid: int,
    linkset_directory: Path = LINKSETS_DIR,
    id_attribute: str = CTL_ID_ATTRIBUTE,
    direction: str = CTL_DIRECTION,
) -> int:
    """Extend STRING network with drug interactions using CyTargetLinker."""
    linksets = sorted(glob.glob(os.path.join(str(linkset_directory), "*.xgmml")))
    if not linksets:
        raise FileNotFoundError(f"No .xgmml linkset files found in: {linkset_directory}")
    print(f"Extending with {len(linksets)} linkset(s): {[os.path.basename(f) for f in linksets]}")

    p4c.set_current_network(network_suid)
    networks_before = set(p4c.get_network_list())

    p4c.commands.commands_post(
        f'cytargetlinker extend idAttribute="{id_attribute}" '
        f'linkSetFiles="{",".join(linksets)}" '
        f"direction={direction} network=current"
    )

    new_networks = set(p4c.get_network_list()) - networks_before
    if new_networks:
        print(f"New network created: {sorted(new_networks)}")
    extended_suid = p4c.get_network_suid()
    print(
        f"Extended network SUID: {extended_suid}, "
        f"{p4c.get_node_count()} nodes, {p4c.get_edge_count()} edges."
    )
    return extended_suid


def color_drugs(
    string_suid: int,
    extended_suid: int,
    color_column: str = NODE_COLOR_COLUMN,
    shape_column: str = NODE_SHAPE_COLUMN,
    size_column: str = NODE_SIZE_COLUMN,
) -> None:
    """Propagate colors from STRING network and configure drug nodes (DIAMOND, white, size 60)."""
    src = p4c.get_table_columns("node", columns=["name", color_column, size_column], network=string_suid)
    color_by_name = dict(zip(src["name"], src[color_column]))
    size_by_name = dict(zip(src["name"], src[size_column]))

    ext = p4c.get_table_columns("node", columns=["name"], network=extended_suid)
    ext_data = pd.DataFrame({
        "SUID": list(ext.index),
        color_column: [color_by_name.get(n, DARK_GRAY) for n in ext["name"]],
        shape_column: [GENE_SHAPE if n in color_by_name else DRUG_SHAPE for n in ext["name"]],
        size_column: [size_by_name.get(n, DRUG_SIZE) for n in ext["name"]],
    })
    p4c.load_table_data(
        ext_data,
        data_key_column="SUID",
        table="node",
        table_key_column="SUID",
        network=extended_suid,
    )
    n_drugs = sum(n not in color_by_name for n in ext["name"])
    print(f"Configured {len(ext)} nodes in extended network ({n_drugs} drug compounds).")


def main() -> None:
    """Execute Step 2: Build STRING network, MCODE cluster, and extend with CyTargetLinker."""
    print("=== Step 2: Creating & Extending Network in Cytoscape ===")
    check_cytoscape(delete_existing_networks=DELETE_EXISTING_NETWORKS)

    genes = load_unique_genes(UNIQUE_GENES_FILE)
    string_suid = build_string_network(genes)
    mcode_col = run_mcode(string_suid)
    add_cluster_colors(string_suid, mcode_col)
    extended_suid = extend_with_cytargetlinker(string_suid)
    color_drugs(string_suid, extended_suid)
    print("Network creation completed successfully.")


if __name__ == "__main__":
    main()
