"""Step 4: Apply visual style mappings in Cytoscape.

Executed after `drug_target_distance.py` so that `border_color` and `border_width`
mappings apply to both the base and extended networks.
"""

from typing import Dict

from config import (
    BORDER_COLOR_COLUMN,
    BORDER_WIDTH_COLUMN,
    EXTENDED_NETWORK_TITLE,
    LABEL_COLUMN,
    NETWORK_TITLE,
    NODE_COLOR_COLUMN,
    NODE_SHAPE_COLUMN,
    NODE_SIZE_COLUMN,
    STYLE_NAME,
)
from utils.cytoscape_utils import apply_visual_style, check_cytoscape, find_networks

NETWORK_TITLES = [NETWORK_TITLE, EXTENDED_NETWORK_TITLE]


def main() -> None:
    """Apply the visual style to all active project networks in Cytoscape."""
    print("=== Step 4: Applying Cytoscape Visual Styles ===")
    check_cytoscape(delete_existing_networks=False)

    found_networks = find_networks(NETWORK_TITLES)
    for title, suid in found_networks.items():
        apply_visual_style(suid, style_name=STYLE_NAME)
        print(
            f"Styled '{title}' (SUID {suid}): "
            f"fill <- '{NODE_COLOR_COLUMN}', shape <- '{NODE_SHAPE_COLUMN}', "
            f"size <- '{NODE_SIZE_COLUMN}', label <- '{LABEL_COLUMN}', "
            f"border <- '{BORDER_COLOR_COLUMN}', border_width <- '{BORDER_WIDTH_COLUMN}'."
        )
    print("Visual styling applied successfully.")


if __name__ == "__main__":
    main()
