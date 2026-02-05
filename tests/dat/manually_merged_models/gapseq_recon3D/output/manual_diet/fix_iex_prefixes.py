#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import cobra

ROOT = Path("tests/dat/manually_merged_models/gapseq_recon3D")
INPUT = ROOT / "output/manual_diet/merged_model_2025_prefixed_normalized_diet.xml"
OUTPUT = ROOT / "output/manual_diet/merged_model_2025_prefixed_normalized_diet_fixIEX.xml"
REPORT = ROOT / "output/manual_diet/merged_model_2025_prefixed_normalized_diet_fixIEX_report.txt"

PREFIXES = ("H_", "M_")


def _infer_prefix_from_mets(rxn: cobra.Reaction) -> str | None:
    prefixes = set()
    for met in rxn.metabolites:
        for pref in PREFIXES:
            if met.id.startswith(pref):
                prefixes.add(pref)
                break
    if len(prefixes) == 1:
        return prefixes.pop()
    return None


def _strip_existing_iex_prefix(rxn_id: str) -> tuple[str | None, str]:
    if rxn_id.startswith("H_IEX_"):
        return "H_", rxn_id[len("H_"):]
    if rxn_id.startswith("M_IEX_"):
        return "M_", rxn_id[len("M_"):]
    return None, rxn_id


def main() -> None:
    model = cobra.io.read_sbml_model(str(INPUT))

    changed = []
    skipped = []
    unchanged = 0

    for rxn in list(model.reactions):
        if "IEX_" not in rxn.id:
            continue

        old_id = rxn.id
        existing_prefix, base_id = _strip_existing_iex_prefix(old_id)
        inferred = _infer_prefix_from_mets(rxn)

        if inferred is None:
            skipped.append(old_id)
            continue

        new_id = f"{inferred}{base_id}" if not base_id.startswith(inferred) else base_id
        if new_id == old_id:
            unchanged += 1
            continue

        if new_id in model.reactions:
            skipped.append(old_id)
            continue

        rxn.id = new_id
        changed.append((old_id, new_id))

    cobra.io.write_sbml_model(model, str(OUTPUT))

    report_lines = [
        f"input: {INPUT}",
        f"output: {OUTPUT}",
        f"changed: {len(changed)}",
        f"unchanged: {unchanged}",
        f"skipped: {len(skipped)}",
        "",
        "changed_ids:",
    ]
    report_lines += [f"{old} -> {new}" for old, new in changed]
    report_lines += ["", "skipped_ids:"]
    report_lines += skipped

    REPORT.write_text("\n".join(report_lines))


if __name__ == "__main__":
    main()
