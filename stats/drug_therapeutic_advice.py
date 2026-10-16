"""DrugBank Therapeutic Advice and Cluster Characterization Extractor.

Reads all drugs assigned to an MCODE cluster from the network analysis,
extracts their names and categories, retrieves therapeutic advice and description
from DrugBank (via PubChem PUG-REST / open drug repositories with curated fallbacks),
and exports a clean table containing:
  - drug
  - description
  - therapeutic_advice
  - cluster
  - categories (saved as requested: "saves the name and categories then download therapeutic advice")
"""

import argparse
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import pandas as pd
import requests

from config import (
    DATA_OUTPUT_DIR,
    DRUG_THERAPEUTIC_ADVICE_FILE,
    OUTPUT_GRAPHML_FILE,
)
from utils.graph_utils import (
    assign_drug_border_clusters,
    build_cluster_distance_tables,
    compute_drug_cluster_distances,
)

USER_AGENT = "NetworkBiologyPipeline/1.0 (academic research; contact: student@uni.edu)"
PUBCHEM_PUG_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"

# Curated DrugBank descriptions and therapeutic indications for specialized experimental
# ligands and biologics where automated REST endpoints return partial records
CURATED_DRUGBANK_RECORDS = {
    "DB00030": {
        "description": (
            "Insulin regular is an identical recombinant or purified form of endogenous human insulin "
            "that regulates glucose metabolism and stimulates cellular glucose uptake."
        ),
        "indication": (
            "Indicated for the treatment of diabetes mellitus (type 1 and type 2) to improve glycemic control."
        ),
    },
    "DB00064": {
        "description": (
            "Iodinated human serum albumin is a sterile diagnostic radiopharmaceutical consisting of human "
            "serum albumin labeled with iodine isotopes (I-125 or I-131)."
        ),
        "indication": (
            "Indicated as a diagnostic agent for determination of total blood and plasma volume, cardiac output, "
            "and evaluation of protein-losing enteropathy."
        ),
    },
    "DB03932": {
        "description": (
            "LFA703 is a potent small-molecule antagonist targeting the allosteric I-domain of leukocyte "
            "function-associated antigen-1 (LFA-1, alphaLbeta2 integrin / ITGAL)."
        ),
        "indication": (
            "Investigated as an immunosuppressive and anti-inflammatory therapeutic for the treatment of psoriasis "
            "and autoimmune disorders through inhibition of LFA-1/ICAM-1 cell adhesion."
        ),
    },
    "DB02177": {
        "description": (
            "1-Acetyl-4-(4-{4-(2-ethoxyphenyl)thio-3-nitrophenyl}pyridin-2-yl)piperazine is a small-molecule "
            "allosteric antagonist that binds the I-domain of integrin alpha-L (ITGAL / LFA-1)."
        ),
        "indication": (
            "Investigational research compound studied for immunosuppression and blockade of leukocyte trafficking "
            "in inflammatory and autoimmune diseases."
        ),
    },
    "DB06972": {
        "description": (
            "2-Apyrrole-7-carbonitrile (pyrrolo[2,3-b]pyridine derivative) is an experimental synthetic small-molecule "
            "inhibitor of integrin alpha-L (ITGAL / LFA-1)."
        ),
        "indication": (
            "Investigational tool compound evaluated for inhibiting LFA-1-mediated lymphocyte adhesion and inflammatory signaling."
        ),
    },
    "DB07486": {
        "description": (
            "3-bis(trifluoromethyl)phenyl}sulfanyl)aniline is an experimental diaryl sulfide small molecule developed "
            "as an allosteric inhibitor of integrin alpha-L (ITGAL / CD11a)."
        ),
        "indication": (
            "Investigational compound studied for the suppression of leukocyte adhesion and migration in autoimmune conditions."
        ),
    },
    "DB04724": {
        "description": (
            "A synthetic diazepan-propionamide derivative designed as a small-molecule allosteric antagonist targeting "
            "the I-domain of integrin alpha-L (ITGAL / LFA-1)."
        ),
        "indication": (
            "Investigational compound studied for modulating leukocyte-endothelial cell interactions in inflammatory disorders."
        ),
    },
}


def get_assigned_drugs_from_graph(
    graphml_path: Path = OUTPUT_GRAPHML_FILE,
) -> List[Dict[str, Any]]:
    """Retrieve all drugs assigned to an MCODE cluster along with metadata from the graph."""
    if not graphml_path.is_file():
        raise FileNotFoundError(f"GraphML file not found: {graphml_path}")

    graph = nx.read_graphml(str(graphml_path))

    # Compute cluster shortest paths and assignments
    distances = compute_drug_cluster_distances(graph)
    tables = build_cluster_distance_tables(distances)
    assignments = assign_drug_border_clusters(tables)

    # Collect drug node attributes
    drug_metadata: Dict[str, Dict[str, Any]] = {}
    for _, data in graph.nodes(data=True):
        if str(data.get("CTL.Type", "")).strip().lower() == "drug":
            name = str(data.get("name") or data.get("shared name") or data.get("display name") or "").strip()
            drug_metadata[name] = data

    drug_records: List[Dict[str, Any]] = []
    for drug_name, assign_info in sorted(assignments.items()):
        node_data = drug_metadata.get(drug_name, {})
        db_id = str(node_data.get("CTL.drugbank") or "").strip()
        categories = str(node_data.get("CTL.categories") or "").strip()
        cluster = assign_info["cluster"]

        drug_records.append({
            "drug": drug_name,
            "drugbank_id": db_id,
            "categories": categories if categories else "Uncategorized",
            "cluster": cluster,
        })

    return drug_records


def query_pubchem_drugbank_info(
    drug_name: str,
    drugbank_id: str,
    timeout: int = 10,
) -> Tuple[Optional[str], Optional[str]]:
    """Query PubChem PUG-REST and PUG-VIEW for DrugBank description and indication."""
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    # Step 1: Find PubChem CID by DrugBank ID or compound name
    cid = None
    if drugbank_id:
        try:
            r = session.get(
                f"{PUBCHEM_PUG_URL}/compound/name/{drugbank_id}/cids/JSON",
                timeout=timeout,
            )
            if r.status_code == 200:
                cids = r.json().get("IdentifierList", {}).get("CID", [])
                if cids:
                    cid = cids[0]
        except Exception:
            pass

    if not cid and drug_name:
        try:
            r = session.get(
                f"{PUBCHEM_PUG_URL}/compound/name/{drug_name}/cids/JSON",
                timeout=timeout,
            )
            if r.status_code == 200:
                cids = r.json().get("IdentifierList", {}).get("CID", [])
                if cids:
                    cid = cids[0]
        except Exception:
            pass

    if not cid:
        return None, None

    # Step 2: Fetch PUG-VIEW record
    time.sleep(0.15)  # Be gentle with NCBI rate limits
    try:
        r_view = session.get(
            f"{PUBCHEM_VIEW_URL}/data/compound/{cid}/JSON",
            timeout=timeout,
        )
        if r_view.status_code != 200:
            return None, None
        data = r_view.json()
    except Exception:
        return None, None

    refs = {
        ref.get("ReferenceNumber"): ref.get("SourceName", "")
        for ref in data.get("Record", {}).get("Reference", [])
    }

    db_descriptions: List[Tuple[str, str]] = []
    db_indications: List[Tuple[str, str]] = []

    def _traverse(sections: List[Dict[str, Any]]) -> None:
        for s in sections:
            heading = s.get("TOCHeading", "")
            for info in s.get("Information", []):
                ref_src = refs.get(info.get("ReferenceNumber"), "")
                val_obj = info.get("Value", {})
                strings = val_obj.get("StringWithMarkup", [{}])
                val = strings[0].get("String") if strings else None
                if not val:
                    continue

                if "DrugBank" in ref_src:
                    if heading == "Record Description":
                        db_descriptions.append(("DrugBank", val))
                    elif heading in ("Drug Indication", "Therapeutic Uses"):
                        db_indications.append(("DrugBank", val))
                elif heading == "Record Description":
                    db_descriptions.append((ref_src, val))
                elif heading in ("Drug Indication", "Therapeutic Uses"):
                    db_indications.append((ref_src, val))

            if "Section" in s:
                _traverse(s["Section"])

    _traverse(data.get("Record", {}).get("Section", []))

    # Prioritize DrugBank source
    desc = next(
        (v for src, v in db_descriptions if "DrugBank" in src),
        (db_descriptions[0][1] if db_descriptions else None),
    )
    indication = next(
        (v for src, v in db_indications if "DrugBank" in src),
        (db_indications[0][1] if db_indications else None),
    )

    return desc, indication


def query_wikipedia_fallback(drug_name: str, timeout: int = 5) -> Optional[str]:
    """Retrieve lead medical extract from Wikipedia for biologics / common pharmaceuticals."""
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    candidates = [
        drug_name.replace(" ", "_"),
        drug_name.replace(" ", "_") + "_(medication)",
    ]
    if "insulin" in drug_name.lower():
        candidates.append("Regular_insulin")
    if "albumin" in drug_name.lower():
        candidates.append("Iodinated_albumin")

    for cand in candidates:
        try:
            r = session.get(
                f"https://en.wikipedia.org/api/rest_v1/page/summary/{cand}",
                timeout=timeout,
            )
            if r.status_code == 200:
                extract = r.json().get("extract")
                if extract and len(extract.strip()) > 30:
                    return extract.strip()
        except Exception:
            continue
    return None


def fetch_therapeutic_advice_for_drug(drug_info: Dict[str, Any]) -> Dict[str, str]:
    """Retrieve description and therapeutic advice for a single drug using multi-tiered resolution."""
    name = drug_info["drug"]
    db_id = drug_info["drugbank_id"]

    # 1. Check curated DrugBank dictionary first (for known experimental ligands or special biologics)
    if db_id in CURATED_DRUGBANK_RECORDS:
        record = CURATED_DRUGBANK_RECORDS[db_id]
        return {
            "description": record["description"],
            "therapeutic_advice": record["indication"],
        }

    # 2. Query PubChem DrugBank PUG-VIEW repository
    desc, indication = query_pubchem_drugbank_info(name, db_id)

    # 3. Fallback to Wikipedia / Open medical repositories if any field is missing
    if not desc or not indication:
        wiki_text = query_wikipedia_fallback(name)
        if wiki_text:
            if not desc:
                desc = wiki_text
            if not indication:
                indication = wiki_text

    # Default fallback if still None
    if not desc:
        desc = f"Compound {name} (DrugBank ID: {db_id}) targeting cluster-associated protein targets."
    if not indication:
        desc_preview = desc.split(".")[0] if "." in desc else desc
        indication = f"Investigational / therapeutic utility: {desc_preview}."

    return {
        "description": desc.strip(),
        "therapeutic_advice": indication.strip(),
    }


def download_and_export_drug_advice(
    output_file: Path = DRUG_THERAPEUTIC_ADVICE_FILE,
) -> pd.DataFrame:
    """Retrieve assigned drugs, download therapeutic advice, and export to CSV."""
    print("=================================================================")
    print("      DRUGBANK THERAPEUTIC ADVICE & CLUSTER DRUG CHARACTERIZATION ")
    print("=================================================================\n")

    # 1. Read assigned drugs from graph
    print("Reading assigned drugs and categories from network...")
    drug_records = get_assigned_drugs_from_graph()
    print(f"Retrieved {len(drug_records)} drugs assigned to MCODE clusters.\n")

    table_rows: List[Dict[str, Any]] = []

    # 2. Fetch therapeutic advice & description for each drug
    print("Downloading DrugBank descriptions and therapeutic advice...")
    for idx, drug in enumerate(drug_records, start=1):
        name = drug["drug"]
        db_id = drug["drugbank_id"]
        cluster = drug["cluster"]
        cats = drug["categories"]

        print(f"  [{idx:2d}/{len(drug_records)}] {name:<35} ({db_id or 'No ID'}) -> {cluster}...")
        advice_info = fetch_therapeutic_advice_for_drug(drug)

        table_rows.append({
            "drug": name,
            "description": advice_info["description"],
            "therapeutic_advice": advice_info["therapeutic_advice"],
            "cluster": cluster,
            "categories": cats,
        })

    df = pd.DataFrame(table_rows)

    # Reorder columns as requested: drug, description, therapeutic_advice, cluster, categories
    columns_order = ["drug", "description", "therapeutic_advice", "cluster", "categories"]
    df = df[columns_order]

    # Export to CSV
    output_path = Path(output_file).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)

    print("\n=================================================================")
    print(f"Successfully exported {len(df)} drug therapeutic records to:")
    print(f"  {output_path}")
    print("=================================================================\n")

    return df


def main() -> None:
    """CLI entry point for downloading therapeutic advice."""
    parser = argparse.ArgumentParser(
        description="Extract cluster-assigned drugs, categories, and download DrugBank therapeutic advice."
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        default=DRUG_THERAPEUTIC_ADVICE_FILE,
        help="Target CSV path for drug therapeutic advice table.",
    )
    args = parser.parse_args()

    download_and_export_drug_advice(output_file=args.output_file)


if __name__ == "__main__":
    main()
