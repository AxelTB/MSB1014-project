"""Cytoscape automation, export, and visual style helper functions using py4cytoscape."""

import glob
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd
import py4cytoscape as p4c
from py4cytoscape import commands

from config import (
    BORDER_COLOR_COLUMN,
    BORDER_WIDTH_COLUMN,
    LABEL_COLUMN,
    LEGACY_PATH_EDGE_TYPE,
    NODE_COLOR_COLUMN,
    NODE_SHAPE_COLUMN,
    NODE_SIZE_COLUMN,
    STYLE_NAME,
)


def check_cytoscape(delete_existing_networks: bool = False) -> None:
    """Verify that Cytoscape is running and CyREST is responsive on localhost."""
    try:
        ping = p4c.cytoscape_ping()
    except Exception as exc:
        sys.exit(f"ERROR: Cannot communicate with Cytoscape: {exc}")

    if not ping:
        sys.exit("ERROR: Cytoscape is not open or not reachable.")

    if delete_existing_networks:
        p4c.delete_all_networks()


def find_networks(titles: List[str]) -> Dict[str, int]:
    """Map network title to SUID for the titles that currently exist in Cytoscape."""
    existing = p4c.get_network_list()
    found = {t: p4c.get_network_suid(t) for t in titles if t in existing}
    if not found:
        sys.exit(f"ERROR: None of {titles} found in Cytoscape. Available: {existing}")
    return found


def export_network_graphml(output_file: Union[str, Path], network: Optional[Union[str, int]] = "current") -> Path:
    """Export the active (or specified) Cytoscape network to a GraphML file.

    Cytoscape needs an absolute path so file resolution does not depend on Cytoscape's CWD.
    """
    out_path = Path(output_file).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.is_file():
        out_path.unlink()

    net_param = "current" if network in ("current", None) else (
        f"SUID:{network}" if isinstance(network, int) else network
    )
    cmd = f'network export network={net_param} options=GraphML outputFile="{out_path}"'
    commands.commands_post(cmd)

    if not out_path.is_file():
        raise FileNotFoundError(f"Export failed, file not created: {out_path}")
    return out_path


def remove_legacy_path_edges(edge_type: str = LEGACY_PATH_EDGE_TYPE) -> None:
    """Delete legacy drug-path edges created by older script versions to avoid corrupting distances."""
    try:
        edge_table = p4c.get_table_columns("edge", columns=["interaction"])
    except Exception:
        return

    if edge_table is None or edge_table.empty or "interaction" not in edge_table:
        return

    legacy_edges = edge_table[edge_table["interaction"] == edge_type]
    if legacy_edges.empty:
        return

    suids = [int(s) for s in legacy_edges.index]
    p4c.delete_edges(suids)
    print(f"Removed {len(suids)} legacy path edge(s) of type '{edge_type}'.")


def apply_visual_style(network_suid: int, style_name: str = STYLE_NAME) -> None:
    """Apply visual style mappings (fill color, shape, size, label, border) to a network."""
    columns = p4c.get_table_column_names("node", network=network_suid)

    # Core required visual columns
    for col in (NODE_COLOR_COLUMN, NODE_SHAPE_COLUMN, NODE_SIZE_COLUMN, LABEL_COLUMN):
        if col not in columns:
            raise KeyError(
                f"Node column '{col}' missing in network {network_suid}. "
                f"Ensure network creation and node attribute assignment ran first."
            )

    if style_name not in p4c.get_visual_style_names():
        p4c.create_visual_style(style_name)
    p4c.set_visual_style(style_name, network=network_suid)

    # Lock width and height so NODE_SIZE drives uniform diameter
    p4c.lock_node_dimensions(True, style_name=style_name)

    # Passthrough property mappings
    p4c.update_style_mapping(
        style_name, p4c.map_visual_property("NODE_FILL_COLOR", NODE_COLOR_COLUMN, "p")
    )
    p4c.update_style_mapping(
        style_name, p4c.map_visual_property("NODE_SHAPE", NODE_SHAPE_COLUMN, "p")
    )
    p4c.update_style_mapping(
        style_name, p4c.map_visual_property("NODE_SIZE", NODE_SIZE_COLUMN, "p")
    )
    p4c.update_style_mapping(
        style_name, p4c.map_visual_property("NODE_LABEL", LABEL_COLUMN, "p")
    )

    # Border color and width mappings if present in the node table
    if BORDER_COLOR_COLUMN in columns:
        p4c.update_style_mapping(
            style_name,
            p4c.map_visual_property(
                "NODE_BORDER_PAINT", BORDER_COLOR_COLUMN, "p", network=network_suid
            ),
        )
    else:
        print(
            f"Node column '{BORDER_COLOR_COLUMN}' missing in network {network_suid}; border color skipped."
        )

    if BORDER_WIDTH_COLUMN in columns:
        p4c.update_style_mapping(
            style_name,
            p4c.map_visual_property(
                "NODE_BORDER_WIDTH", BORDER_WIDTH_COLUMN, "p", network=network_suid
            ),
        )
    else:
        print(
            f"Node column '{BORDER_WIDTH_COLUMN}' missing in network {network_suid}; border width skipped."
        )
