"""Settings and paths shared across all modules in the Network Biology pipeline."""

import os
from pathlib import Path
from typing import Optional

# Base project directories
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
DATA_INPUT_DIR = DATA_DIR / "input"
DATA_OUTPUT_DIR = DATA_DIR / "output"
REFERENCE_DIR = BASE_DIR / "reference"

# Input file paths
CLINVAR_FILTERED_FILE = DATA_INPUT_DIR / "clinvar_variant_filtered.tsv"
VARIANT_GENE_NETWORK_FILE = DATA_INPUT_DIR / "variant_gene_network.tsv"
UNIQUE_GENES_FILE = DATA_INPUT_DIR / "unique_genes.tsv"
LINKSETS_DIR = DATA_INPUT_DIR / "linksets"

# Output file paths
OUTPUT_GRAPHML_FILE = DATA_OUTPUT_DIR / "network.graphml"
COMPLETE_GRAPH_DIR = DATA_OUTPUT_DIR / "drug_cluster_selection"
CLUSTER_TABLES_DIR = COMPLETE_GRAPH_DIR / "cluster_tables"
CLUSTER_TOP_DRUGS_FILE = COMPLETE_GRAPH_DIR / "cluster_top_drugs.csv"
SNP_TABLE_FILE = COMPLETE_GRAPH_DIR / "snp_closest_drugs.csv"
SNP_SPECIFIC_DIR = DATA_OUTPUT_DIR / "snp_specific"
SNP_SPECIFIC_TABLE_FILE = SNP_SPECIFIC_DIR / "snp_closest_drugs.csv"
EDGE_REMOVAL_DIR = DATA_OUTPUT_DIR / "edge_removal"
EDGE_REMOVAL_SNP_TABLE_FILE = EDGE_REMOVAL_DIR / "snp_closest_drugs.csv"
ENRICHMENT_DIR = DATA_OUTPUT_DIR / "cluster_enrichment"
ENRICHMENT_PLOTS_DIR = ENRICHMENT_DIR / "plots"
CLUSTER_CHARACTERIZATION_FILE = DATA_OUTPUT_DIR / "cluster_characterization_summary.csv"
DRUG_THERAPEUTIC_ADVICE_FILE = DATA_OUTPUT_DIR / "drug_therapeutic_advice.csv"

# Reference graph path for validation
REFERENCE_GRAPHML_FILE = REFERENCE_DIR / "reference_network.graphml"

# Cytoscape Network Titles
NETWORK_TITLE = "Gene Network"
EXTENDED_NETWORK_TITLE = f"CTL_{NETWORK_TITLE}"   # CyTargetLinker network

# Visual Style Settings
STYLE_NAME = "Gene Network Style"

# Node table columns
LABEL_COLUMN = "display name"
NODE_COLOR_COLUMN = "node_color"
NODE_SHAPE_COLUMN = "node_shape"
NODE_SIZE_COLUMN = "node_size"
BORDER_COLOR_COLUMN = "border_color"
BORDER_WIDTH_COLUMN = "border_width"

# Colors: one per MCODE cluster, cluster N gets CLUSTER_COLORS[N - 1] (colorblind-friendly)
CLUSTER_COLORS = [
   # '#FFF176',  # 1: Lemon Yellow
    '#FF7043',  # 2: Vibrant Coral (Distinct from orange/vermillion)
    '#BA68C8',  # 3: Orchid Purple
    '#4DD0E1',  # 4: Bright Cyan
    '#81C784',  # 5: Spring Green
    '#7986CB',  # 6: Periwinkle Blue
]
CLUSTER_NAMES = [
    "Red Cluster",
    "Purple Cluster",
    "Cyan Cluster",
    "Green Cluster",
    "Blue Cluster",
]
LIGHT_GRAY = "#D9D9D9"   # nodes outside any cluster
DARK_GRAY = "#FFFFFF"    # drug nodes

# Shapes
GENE_SHAPE = "ELLIPSE"
DRUG_SHAPE = "DIAMOND"

# Sizes
SNP_GENE_SIZE = 70       # input (SNP) genes, i.e. STRING nodes with a query term
DRUG_SIZE = 60
OTHER_SIZE = 40          # added STRING interactors

# Borders
DEFAULT_BORDER_WIDTH = 2
HIGHLIGHTED_BORDER_WIDTH = 6

# Network Generation Parameters
DELETE_EXISTING_NETWORKS = True
SPECIES = "Homo sapiens"
STRING_CONFIDENCE_THRESHOLD = 0.7   # STRING cutoff (0.0 - 1.0)
ADDITIONAL_INTERACTORS = 30         # Maximum additional interactors from STRING

# MCODE Clustering Parameters
MCODE_NODE_SCORE_CUTOFF = 0.3
MCODE_K_CORE = 2
MCODE_MAX_DEPTH = 100
MCODE_HAIRCUT = True
MCODE_FLUFF = False

# CyTargetLinker Parameters
CTL_ID_ATTRIBUTE = "display name"
CTL_DIRECTION = "SOURCES"

# Drug Target Distance Parameters
TOP_N_DRUGS = 3
LEGACY_PATH_EDGE_TYPE = "cluster drug path"

# Standard attribute keys for robustness
QUERY_TERM_KEYS = ("query term", "query_term", "queryterm", "Query Term")
CTL_TYPE_KEYS = ("CTL.Type", "CTL_Type", "ctl.type")


def cluster_number(cluster: object) -> Optional[int]:
    """Number of an MCODE cluster name such as 'Cluster 3' (or 3); None if it has none."""
    last = str(cluster).strip().split()[-1] if str(cluster).strip() else ""
    try:
        return int(float(last))
    except ValueError:
        return None


def cluster_color(cluster: object) -> str:
    """The single color of an MCODE cluster; LIGHT_GRAY for no cluster."""
    number = cluster_number(cluster) if cluster is not None else None
    if number is None:
        return LIGHT_GRAY
    if not 1 <= number <= len(CLUSTER_COLORS):
        raise ValueError(f"Cluster {cluster!r} has no color: only {len(CLUSTER_COLORS)} in CLUSTER_COLORS.")
    return CLUSTER_COLORS[number - 1]



def cluster_name(cluster: object) -> str:
    """Color name of an MCODE cluster (e.g. 'Cluster 2' -> 'Purple Cluster'); falls back to the input."""
    number = cluster_number(cluster)
    if number is None or not 1 <= number <= len(CLUSTER_NAMES):
        return str(cluster)
    return CLUSTER_NAMES[number - 1]
