from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import sys

import cobra

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix

def prefix_merged_model(
    merged_path: Path,
    model1_path: Path,
    model2_path: Path,
    output_path: Path,
) -> None:
    merged = cobra.io.read_sbml_model(merged_path)
    model1 = cobra.io.read_sbml_model(model1_path)
    model2 = cobra.io.read_sbml_model(model2_path)

    model1_mets = {met.id for met in model1.metabolites}
    model2_mets = {met.id for met in model2.metabolites}
    model1_rxns = {rxn.id for rxn in model1.reactions}
    model2_rxns = {rxn.id for rxn in model2.reactions}

    for met in merged.metabolites:
        if met.id.startswith(("H_", "M_")):
            continue
        if met.compartment == "f":
            continue
        in_m1 = met.id in model1_mets
        in_m2 = met.id in model2_mets
        if in_m1 and not in_m2:
            met.id = f"H_{met.id}"
        elif in_m2 and not in_m1:
            met.id = f"M_{met.id}"
        else:
            if met.compartment in ("b", "p"):
                met.id = f"M_{met.id}"
            else:
                met.id = f"H_{met.id}"

    for rxn in merged.reactions:
        if rxn.id.startswith(("H_", "M_")):
            continue
        in_m1 = rxn.id in model1_rxns
        in_m2 = rxn.id in model2_rxns
        if in_m1 and not in_m2:
            rxn.id = f"H_{rxn.id}"
            continue
        if in_m2 and not in_m1:
            rxn.id = f"M_{rxn.id}"
            continue

        if any(met.compartment in ("b", "p") for met in rxn.metabolites):
            rxn.id = f"M_{rxn.id}"
        elif any(met.id.startswith("H_") for met in rxn.metabolites):
            rxn.id = f"H_{rxn.id}"

    cobra.io.write_sbml_model(merged, output_path)


def infer_manual_pairs(manual_path: Path, output_csv: Path) -> None:
    manual = cobra.io.read_sbml_model(manual_path)

    feed_set = {met.id for met in manual.metabolites if met.compartment == "f"}
    host_map: dict[str, set[str]] = defaultdict(set)
    micro_map: dict[str, set[str]] = defaultdict(set)

    for rxn in manual.reactions:
        if rxn.id.startswith("H_IEX_"):
            feed_in_rxn = [met.id for met in rxn.metabolites if met.compartment == "f"]
            host_mets = [
                met.id
                for met in rxn.metabolites
                if met.id.startswith("H_") and met.compartment == "e"
            ]
            for feed_id in feed_in_rxn:
                for host_id in host_mets:
                    host_map[feed_id].add(host_id)
        elif rxn.id.startswith("M_IEX_"):
            feed_in_rxn = [met.id for met in rxn.metabolites if met.compartment == "f"]
            if len(feed_in_rxn) > 1 and "h[f]" in feed_in_rxn:
                feed_in_rxn = [fid for fid in feed_in_rxn if fid != "h[f]"]
            micro_mets = [
                met.id
                for met in rxn.metabolites
                if met.id.startswith("M_") and met.compartment == "e"
            ]
            for feed_id in feed_in_rxn:
                for micro_id in micro_mets:
                    micro_map[feed_id].add(micro_id)

    rows = []
    for feed_id in sorted(feed_set):
        if not host_map.get(feed_id) or not micro_map.get(feed_id):
            continue
        for host_id in sorted(host_map[feed_id]):
            for micro_id in sorted(micro_map[feed_id]):
                rows.append((feed_id, host_id, micro_id))

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", encoding="utf-8") as handle:
        handle.write("feed_met,host_met,micro_met\n")
        for feed_id, host_id, micro_id in rows:
            handle.write(f"{feed_id},{host_id},{micro_id}\n")


def normalize_prefixed_model(
    prefixed_path: Path,
    model1_path: Path,
    model2_path: Path,
    output_path: Path,
) -> None:
    merged = cobra.io.read_sbml_model(prefixed_path)
    model1 = cobra.io.read_sbml_model(model1_path)
    model2 = cobra.io.read_sbml_model(model2_path)

    model1_base = {
        handle_metabolites_prefix_suffix(met.id)
        for met in model1.metabolites
        if handle_metabolites_prefix_suffix(met.id)
    }
    model2_base = {
        handle_metabolites_prefix_suffix(met.id)
        for met in model2.metabolites
        if handle_metabolites_prefix_suffix(met.id)
    }
    model1_rxns = {rxn.id for rxn in model1.reactions}
    model2_rxns = {rxn.id for rxn in model2.reactions}

    for met in merged.metabolites:
        if met.compartment == "f":
            continue
        base_id = handle_metabolites_prefix_suffix(met.id)
        if not base_id:
            continue
        if base_id.startswith(("H_", "M_")):
            base_id = base_id.split("_", 1)[1]
        in_m1 = base_id in model1_base
        in_m2 = base_id in model2_base
        if in_m1 and not in_m2:
            met.id = f"H_{base_id}[{met.compartment}]"
        elif in_m2 and not in_m1:
            met.id = f"M_{base_id}[{met.compartment}]"
        else:
            if met.compartment in ("b", "p"):
                met.id = f"M_{base_id}[{met.compartment}]"
            else:
                met.id = f"H_{base_id}[{met.compartment}]"

    for rxn in merged.reactions:
        base_id = rxn.id
        if base_id.startswith(("H_", "M_")):
            base_id = base_id.split("_", 1)[1]
        in_m1 = base_id in model1_rxns
        in_m2 = base_id in model2_rxns
        if in_m1 and not in_m2:
            rxn.id = f"H_{base_id}"
            continue
        if in_m2 and not in_m1:
            rxn.id = f"M_{base_id}"
            continue

        if any(met.compartment in ("b", "p") for met in rxn.metabolites):
            rxn.id = f"M_{base_id}"
        else:
            rxn.id = f"H_{base_id}"

    cobra.io.write_sbml_model(merged, output_path)

def main() -> None:
    base = Path("tests/dat/manually_merged_models/gapseq_recon3D")
    merged_path = base / "merged_model_2025.xml"
    model1_path = base / "M1_recon3D_301_modified.xml"
    model2_path = base / "M2_bacterial_model.xml"
    output_prefixed = base / "merged_model_2025_prefixed.xml"
    output_prefixed_normalized = base / "merged_model_2025_prefixed_normalized.xml"
    manual_path = base / "Host_Microbiome_metamodel_2025.xml"
    output_pairs = base / "output" / "manual_inferred_pairs.csv"

    prefix_merged_model(merged_path, model1_path, model2_path, output_prefixed)
    normalize_prefixed_model(
        output_prefixed, model1_path, model2_path, output_prefixed_normalized
    )
    infer_manual_pairs(manual_path, output_pairs)


if __name__ == "__main__":
    main()
