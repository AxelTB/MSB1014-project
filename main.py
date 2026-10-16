"""Master execution script for the Network Biology pipeline.

Executes all pipeline stages in sequence:
  1. export_genes.py: Extract unique genes from ClinVar variants.
  2. create_network.py: Query STRING, cluster with MCODE, extend with CyTargetLinker.
  3. drug_target_distance.py: Compute shortest paths, export tables, and write drug borders.
  3b. snp_drug_distance.py: Find closest drugs per SNP-affected gene after node removal.
  4. set_graph_style.py: Apply visual style in Cytoscape (using node and border attributes).
  5. pathway_analysis.py: Retrieve clusters from graph, perform GO over-representation analysis.
  6. validate_pipeline.py: Download the new Cytoscape graph and verify coincidence with reference.
"""

import sys
from export_genes import main as export_genes_main
from create_network import main as create_network_main
from drug_target_distance import main as drug_target_distance_main
from snp_drug_distance import main as snp_drug_distance_main
from set_graph_style import main as set_graph_style_main
from string_functional_enrichment import main as pathway_analysis_main


def run_pipeline() -> None:
    """Run all pipeline stages end-to-end."""
    print("=================================================================")
    print("       STARTING COMPLETE NETWORK BIOLOGY PIPELINE               ")
    print("=================================================================\n")

    # Step 1: Gene extraction
    export_genes_main()
    print()

    # Step 2: Network creation, clustering, extension
    create_network_main()
    print()

    # Step 3: Drug distances & border attribute computation (MUST run before styling)
    drug_target_distance_main()
    print()

    # Step 3b: SNP-specific closest drugs after gene removal
    snp_drug_distance_main()
    print()

    # Step 4: Apply visual style (including fill color, shape, size, and borders)
    set_graph_style_main()
    print()

    # Step 5: STRING functional enrichment analysis on MCODE clusters
    pathway_analysis_main()
    print()


    print("=================================================================")
    print("       PIPELINE COMPLETED!                                       ")
    print("=================================================================")


if __name__ == "__main__":
    run_pipeline()
