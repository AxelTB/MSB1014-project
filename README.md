# Network Biology: Disease Variant Network and Drug Target Proximity Analysis

## Overview
This project implements an automated, reproducible network biology workflow that investigates disease-associated genetic variants (specifically Familial Hypercholesterolemia and related lipid disorders), builds functional protein-protein interaction (PPI) networks, identifies modular protein complexes via MCODE clustering, expands the network with drug-target interactions using DrugBank linksets via CyTargetLinker, and evaluates candidate drug repurposing potential using topological shortest-path distance metrics.

---

## AI Usage Notice
> **Notice**: This codebase and documentation were refactored, modularized, and validated with the assistance of an AI programming assistant (Copilot SDK in VS Code). The underlying biological pipeline logic, MCODE clustering parameters, topological distance ranking methodology, and verification criteria were preserved and verified against baseline reference exports.

---

## Directory Structure
```text
Project_repo/
├── .venv/                         # Python virtual environment (Python 3.12)
├── config.py                      # Global configuration parameters, styling, and paths
├── export_genes.py                # Step 1: Extract unique disease genes from variant data
├── create_network.py              # Step 2: STRING query, MCODE clustering, CyTargetLinker extension
├── drug_target_distance.py        # Step 3: Shortest-path distance metrics, drug ranking, and border styling
├── set_graph_style.py             # Step 4: Apply visual styling mappings in Cytoscape
├── pathway_analysis.py            # Step 5: Retrieve clusters from Cytoscape, perform GO & pathway ORA
├── validate_pipeline.py           # Step 6: Download live graph and verify coincidence with reference
├── main.py                        # Master pipeline orchestrator (Steps 1 to 6)
├── stats/                         # Standalone analyses
│   ├── __init__.py
│   ├── drug_therapeutic_advice.py # Drug therapeutic advice and cluster characterization
│   ├── mutation_stats.py          # Per-gene mutation counts
│   └── cluster_stats.py           # Protein counts and per-cluster MCODE topology metrics
├── requirements.txt               # Pinned Python package dependencies
├── reference/
│   └── reference_network.graphml  # Ground-truth reference GraphML export
├── utils/                         # Reusable auxiliary modules
│   ├── __init__.py
│   ├── cytoscape_utils.py         # Py4Cytoscape automation, export, and style mappings
│   ├── graph_utils.py             # NetworkX analysis, node classification, shortest paths, rankings
│   └── data_utils.py              # TSV/CSV parsing and export helper functions
└── data/
    ├── input/                     # Immutable input files for the pipeline
    │   ├── clinvar_variant_filtered.tsv  # Filtered ClinVar variants
    │   ├── variant_gene_network.tsv      # Cleaned variant-to-gene associations
    │   ├── unique_genes.tsv              # Target gene list for STRING query
    │   └── linksets/                     # CyTargetLinker linkset files (.xgmml)
    │       └── drugbank4-2.xgmml
    └── output/                    # Generated pipeline outputs (created at runtime)
        ├── network.graphml        # Exported GraphML from Cytoscape
        ├── drug_cluster_selection/        # Distance tables and SNP closest drug rankings
        │   ├── cluster_tables/    # (legacy) per-cluster CSVs; replaced by cluster_top_drugs.csv
        │   └── snp_closest_drugs.csv
        ├── cluster_enrichment/    # GO and pathway over-representation analysis per cluster
        │   └── plots/             # Pathway vs gene count scatter plot graphs (.png)
        ├── cluster_characterization_summary.csv # High-level characterization summary table across all clusters
        └── edge_removal/          # Sensitivity analysis removing direct SNP-drug edges
            └── snp_closest_drugs.csv
```

---

## Pipeline Workflow

1. **Gene Extraction (`export_genes.py`)**:
   - Parses the ClinVar variant-to-gene table (`data/input/variant_gene_network.tsv`).
   - Handles multi-gene annotations (separated by `|`) and outputs a clean, deduplicated gene list to `data/input/unique_genes.tsv`.

2. **Network Creation & Extension (`create_network.py`)**:
   - Checks Cytoscape connection on `localhost:1234` via CyREST.
   - Queries STRING directly for *Homo sapiens* PPI interactions (confidence cutoff = 0.7, limit = 30 additional interactors).
   - Runs MCODE clustering (`nodeScoreCutoff=0.3`, `kCore=2`, `maxDepth=100`, `haircut=True`) and assigns cluster membership.
   - Assigns colorblind-friendly fill colors for clusters, larger ellipses for query genes, and medium ellipses for interactors.
   - Extends the network with drug-target interactions using CyTargetLinker and `data/input/linksets/drugbank4-2.xgmml`.
   - Formats drug nodes (white diamonds, size 60).

3. **Topological Distance Analysis & Drug Border Assignment (`drug_target_distance.py`)**:
   - **Crucial Sequence Requirement**: Must run **before** visual style application so that `border_color` and `border_width` columns are populated in Cytoscape.
   - Exports the live extended network to `data/output/network.graphml`.
   - Computes undirected shortest-path distances from every drug compound to all proteins in each MCODE cluster.
   - Generates one table (`data/output/drug_cluster_selection/cluster_top_drugs.csv`) with the top 5 closest drugs per cluster, their DrugBank category, and min/max/mean distances.
   - Identifies top candidate drugs for each disease-associated SNP protein and writes `snp_closest_drugs.csv`.
   - Performs an edge-removal sensitivity analysis (removing direct SNP-drug edges) in `data/output/edge_removal/snp_closest_drugs.csv`.
   - Computes drug-to-cluster assignment and writes `border_color` and `border_width=6` to Cytoscape.

4. **Visual Style Application (`set_graph_style.py`)**:
   - Applies the visual style (`Gene Network Style`) across Cytoscape networks.
   - Maps `node_color`, `node_shape`, `node_size`, `display name`, `border_color`, and `border_width`.

5. **Pathway & GO Over-Representation Analysis (`pathway_analysis.py`)**:
   - Retrieves protein clusters from the Cytoscape graph (or exported GraphML).
   - Performs Gene Ontology (Biological Process, Cellular Component, Molecular Function) and pathway over-representation analysis (ORA).
   - Prints formatted cluster characterization reports with top enriched terms and FDR values.
   - Generates and saves pathway vs. gene count scatter plots for each cluster in `data/output/cluster_enrichment/plots/`.
   - Exports high-level cluster characterization summary to `data/output/cluster_characterization_summary.csv` and detailed enrichment tables to `data/output/cluster_enrichment/`.

6. **Graph Verification (`validate_pipeline.py`)**:
   - Downloads the active network from Cytoscape.
   - Validates 100% coincidence against `reference/reference_network.graphml` (verifying exact node counts, edge counts, node identities, edge connectivity, and key attributes).

---

## Prerequisites and Setup

### 1. Cytoscape Environment
- Cytoscape 3.10+ running with CyREST enabled on `localhost:1234`.
- Installed Cytoscape Apps:
  - **stringApp** (for querying STRING PPI database)
  - **MCODE** (for cluster detection)
  - **CyTargetLinker** (for drug-target network extension)

### 2. Python Environment
Activate the project's Python virtual environment:
```bash
source .venv/bin/activate
pip install -r requirements.txt
```

---

## Running the Pipeline

### Option A: Run the Complete Pipeline (Recommended)
To execute all stages from start to finish with automated validation:
```bash
.venv/bin/python main.py
```

### Option B: Run Individual Pipeline Steps
```bash
# Step 1: Extract unique genes
.venv/bin/python export_genes.py

# Step 2: Create PPI network, cluster, and extend with drugs
.venv/bin/python create_network.py

# Step 3: Compute distances and set drug border attributes (run before styling)
.venv/bin/python drug_target_distance.py

# Step 4: Apply visual style in Cytoscape
.venv/bin/python set_graph_style.py

# Step 5: Perform pathway and GO over-representation analysis on clusters
.venv/bin/python pathway_analysis.py

# Step 6: Verify that the network coincides with the reference graph
.venv/bin/python validate_pipeline.py
```

## Standalone Statistics

Run these analyses from the project root:
```bash
# Export per-gene total and nonsense mutation counts
.venv/bin/python -m stats.mutation_stats

# Print protein/cluster counts and export per-cluster metrics (cluster_stats.csv)
.venv/bin/python -m stats.cluster_stats

# Export drug therapeutic advice for drugs in the analyzed network
.venv/bin/python -m stats.drug_therapeutic_advice
```
Both scripts write their CSV results to `data/output/`.
