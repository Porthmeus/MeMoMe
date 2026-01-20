from __future__ import annotations

from pathlib import Path
from collections import defaultdict
import sys

import cobra

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix


def base_id(met_id: str) -> str | None:
    return handle_metabolites_prefix_suffix(met_id)


def strip_prefix(met_id: str, prefixes: tuple[str, ...]) -> str:
    for prefix in prefixes:
        if met_id.startswith(prefix):
            return met_id[len(prefix) :]
    return met_id


def collect_base_sets(model: cobra.Model, host_prefix: str, micro_prefix: str, feed_comp: str) -> dict[str, set[str]]:
    host = set()
    micro = set()
    feed = set()
    for met in model.metabolites:
        if met.compartment == feed_comp:
            bid = base_id(met.id)
            if bid:
                feed.add(bid)
            continue
        if met.id.startswith(host_prefix):
            bid = base_id(met.id[len(host_prefix) :])
            if bid:
                host.add(bid)
        elif met.id.startswith(micro_prefix):
            bid = base_id(met.id[len(micro_prefix) :])
            if bid:
                micro.add(bid)
    return {"host": host, "micro": micro, "feed": feed}


def infer_auto_pairs(model: cobra.Model) -> set[tuple[str, str]]:
    pairs = set()
    translation_mets = [m for m in model.metabolites if m.compartment == "t"]
    for met in translation_mets:
        host_ids = set()
        micro_ids = set()
        for rxn in met.reactions:
            if len(rxn.metabolites) != 2:
                continue
            other = next(m for m in rxn.metabolites if m.id != met.id)
            if other.id.startswith("model1_"):
                bid = base_id(other.id[len("model1_") :])
                if bid:
                    host_ids.add(bid)
            elif other.id.startswith("model2_"):
                bid = base_id(other.id[len("model2_") :])
                if bid:
                    micro_ids.add(bid)
        if host_ids and micro_ids:
            for host in host_ids:
                for micro in micro_ids:
                    pairs.add((host, micro))
    return pairs


def infer_manual_pairs_from_csv(csv_path: Path) -> set[tuple[str, str]]:
    pairs = set()
    with csv_path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        for line in handle:
            feed_id, host_id, micro_id = line.strip().split(",")
            host = strip_prefix(host_id, ("H_",))
            micro = strip_prefix(micro_id, ("M_",))
            host_base = base_id(host)
            micro_base = base_id(micro)
            if host_base and micro_base:
                pairs.add((host_base, micro_base))
    return pairs


def write_set_csv(path: Path, header: str, values: set[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"{header}\n")
        for value in sorted(values):
            handle.write(f"{value}\n")


def write_pairs_csv(path: Path, pairs: set[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("host_base,micro_base\n")
        for host, micro in sorted(pairs):
            handle.write(f"{host},{micro}\n")


def main() -> None:
    base = Path("tests/dat/manually_merged_models/gapseq_recon3D")
    output_dir = base / "output" / "compare_auto_manual"
    output_dir.mkdir(parents=True, exist_ok=True)

    auto_path = base / "output" / "automatically_merged_metamodel.xml"
    manual_path = base / "merged_model_2025_prefixed_normalized.xml"
    manual_pairs_path = base / "output" / "manual_inferred_pairs.csv"

    auto = cobra.io.read_sbml_model(auto_path)
    manual = cobra.io.read_sbml_model(manual_path)

    auto_sets = collect_base_sets(auto, "model1_", "model2_", "t")
    manual_sets = collect_base_sets(manual, "H_", "M_", "f")

    auto_pairs = infer_auto_pairs(auto)
    manual_pairs = infer_manual_pairs_from_csv(manual_pairs_path)

    # Set diffs
    host_auto_only = auto_sets["host"] - manual_sets["host"]
    host_manual_only = manual_sets["host"] - auto_sets["host"]
    micro_auto_only = auto_sets["micro"] - manual_sets["micro"]
    micro_manual_only = manual_sets["micro"] - auto_sets["micro"]
    feed_auto_only = auto_sets["feed"] - manual_sets["feed"]
    feed_manual_only = manual_sets["feed"] - auto_sets["feed"]

    pairs_auto_only = auto_pairs - manual_pairs
    pairs_manual_only = manual_pairs - auto_pairs
    pairs_overlap = auto_pairs & manual_pairs

    # Write CSVs
    write_set_csv(output_dir / "host_only_auto.csv", "host_base", host_auto_only)
    write_set_csv(output_dir / "host_only_manual.csv", "host_base", host_manual_only)
    write_set_csv(output_dir / "micro_only_auto.csv", "micro_base", micro_auto_only)
    write_set_csv(output_dir / "micro_only_manual.csv", "micro_base", micro_manual_only)
    write_set_csv(output_dir / "feed_only_auto.csv", "feed_base", feed_auto_only)
    write_set_csv(output_dir / "feed_only_manual.csv", "feed_base", feed_manual_only)

    write_pairs_csv(output_dir / "pairs_only_auto.csv", pairs_auto_only)
    write_pairs_csv(output_dir / "pairs_only_manual.csv", pairs_manual_only)
    write_pairs_csv(output_dir / "pairs_overlap.csv", pairs_overlap)

    # Summary
    summary_path = output_dir / "summary.txt"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("Auto vs Manual Comparison Summary\n")
        handle.write(f"auto metabolites: {len(auto.metabolites)}\n")
        handle.write(f"manual metabolites: {len(manual.metabolites)}\n")
        handle.write(f"auto reactions: {len(auto.reactions)}\n")
        handle.write(f"manual reactions: {len(manual.reactions)}\n")
        handle.write("\n")
        handle.write(f"host base auto: {len(auto_sets['host'])}\n")
        handle.write(f"host base manual: {len(manual_sets['host'])}\n")
        handle.write(f"host auto-only: {len(host_auto_only)}\n")
        handle.write(f"host manual-only: {len(host_manual_only)}\n")
        handle.write("\n")
        handle.write(f"micro base auto: {len(auto_sets['micro'])}\n")
        handle.write(f"micro base manual: {len(manual_sets['micro'])}\n")
        handle.write(f"micro auto-only: {len(micro_auto_only)}\n")
        handle.write(f"micro manual-only: {len(micro_manual_only)}\n")
        handle.write("\n")
        handle.write(f"feed base auto: {len(auto_sets['feed'])}\n")
        handle.write(f"feed base manual: {len(manual_sets['feed'])}\n")
        handle.write(f"feed auto-only: {len(feed_auto_only)}\n")
        handle.write(f"feed manual-only: {len(feed_manual_only)}\n")
        handle.write("\n")
        handle.write(f"pairs auto: {len(auto_pairs)}\n")
        handle.write(f"pairs manual: {len(manual_pairs)}\n")
        handle.write(f"pairs overlap: {len(pairs_overlap)}\n")
        handle.write(f"pairs auto-only: {len(pairs_auto_only)}\n")
        handle.write(f"pairs manual-only: {len(pairs_manual_only)}\n")

    print("Wrote comparison artifacts to", output_dir)


if __name__ == "__main__":
    main()
