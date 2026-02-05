import pandas as pd
import numpy as np
from pathlib import Path


ROOT = Path("tests/dat/manually_merged_models/gapseq_recon3D/output")
AUTO_MODEL_PATH = ROOT / "automatically_merged_metamodel.xml"
MANUAL_MODEL_PATH = ROOT / "merged_model_2025_prefixed_normalized_diet_fixIEX.xml"
IN_DIR = ROOT / "fva_compare"

# Load FVA results
auto_host_fva_baseline_path = IN_DIR / "auto_host_fva_baseline.csv"
auto_host_fva_baseline = pd.read_csv(auto_host_fva_baseline_path, index_col=0)
manual_host_fva_baseline_path = IN_DIR / "manual_host_fva_baseline.csv"
manual_host_fva_baseline = pd.read_csv(manual_host_fva_baseline_path, index_col=0)
auto_host_fva_closed_path = IN_DIR / "auto_host_fva_closed.csv"
auto_host_fva_closed = pd.read_csv(auto_host_fva_closed_path, index_col=0)
manual_host_fva_closed_path = IN_DIR / "manual_host_fva_closed.csv"
manual_host_fva_closed = pd.read_csv(manual_host_fva_closed_path, index_col=0)

# Canonicalize reaction IDs to compare auto vs manual
def canonicalize_rxn_id(rxn_id: str) -> str:
    base = rxn_id
    for pref in ("model1_", "model2_", "H_", "M_"):
        if base.startswith(pref):
            base = base[len(pref):]
            break
    for pref in ("DM_model1_", "DM_model2_"):
        if base.startswith(pref):
            base = "DM_" + base[len(pref):]
            break
    return base

# Flag/remove loop-participating reactions:
#    - A reaction is considered a loop if abs(min) >= 999 and abs(max) >= 999
#      in the baseline (open microbiome) FVA.
#    - Exclude union of auto/manual loop sets.
loop_threshold = 999
auto_loops = (auto_host_fva_baseline["minimum"].abs() >= loop_threshold) & (auto_host_fva_baseline["maximum"].abs() >= loop_threshold)
manual_loops = (manual_host_fva_baseline["minimum"].abs() >= loop_threshold) & (manual_host_fva_baseline["maximum"].abs() >= loop_threshold)
auto_loop_ids = set(auto_host_fva_baseline.index[auto_loops].tolist())
manual_loop_ids = set(manual_host_fva_baseline.index[manual_loops].tolist())
loops_union = auto_loop_ids | manual_loop_ids
# Exclude loop reactions from all datasets
auto_host_fva_baseline = auto_host_fva_baseline.drop(index=loops_union, errors="ignore")
manual_host_fva_baseline = manual_host_fva_baseline.drop(index=loops_union, errors="ignore")
auto_host_fva_closed = auto_host_fva_closed.drop(index=loops_union, errors="ignore")
manual_host_fva_closed = manual_host_fva_closed.drop(index=loops_union, errors="ignore")
# Save list of looping reactions
with open(IN_DIR / "excluded_loop_reactions.txt", "w") as f:
    for rxn_id in sorted(loops_union):
        f.write(f"{rxn_id}\n")

# Identify blocked reactions (both ranges equal to 0) to avoid division by zero
auto_blocked = (auto_host_fva_baseline["range"] == 0) & (auto_host_fva_closed["range"] == 0)
manual_blocked = (manual_host_fva_baseline["range"] == 0) & (manual_host_fva_closed["range"] == 0)
#exclude blocked reactions from dataframes
auto_host_fva_baseline = auto_host_fva_baseline[~auto_blocked]
manual_host_fva_baseline = manual_host_fva_baseline[~manual_blocked]
auto_host_fva_closed = auto_host_fva_closed[~auto_blocked]
manual_host_fva_closed = manual_host_fva_closed[~manual_blocked]
# Compute relative change, where relative_change = (range_open − range_closed) / range_open)
auto_host_fva_relative_change = (auto_host_fva_baseline["range"] - auto_host_fva_closed["range"]) / auto_host_fva_baseline["range"]
manual_host_fva_relative_change = (manual_host_fva_baseline["range"] - manual_host_fva_closed["range"]) / manual_host_fva_baseline["range"] 

# - Save comparison CSVs + summary
auto_host_fva_comparison = pd.DataFrame({
    "baseline_range": auto_host_fva_baseline["range"],
    "closed_range": auto_host_fva_closed["range"],
    "relative_change": auto_host_fva_relative_change
})
manual_host_fva_comparison = pd.DataFrame({
    "baseline_range": manual_host_fva_baseline["range"],
    "closed_range": manual_host_fva_closed["range"],
    "relative_change": manual_host_fva_relative_change
})

# Build a single comparison dataframe on canonicalized IDs
auto_host_fva_comparison["canonical_id"] = [canonicalize_rxn_id(r) for r in auto_host_fva_comparison.index]
manual_host_fva_comparison["canonical_id"] = [canonicalize_rxn_id(r) for r in manual_host_fva_comparison.index]

auto_host_fva_comparison = auto_host_fva_comparison.set_index("canonical_id", drop=True)
manual_host_fva_comparison = manual_host_fva_comparison.set_index("canonical_id", drop=True)

common_ids = auto_host_fva_comparison.index.intersection(manual_host_fva_comparison.index)
auto_common = auto_host_fva_comparison.loc[common_ids]
manual_common = manual_host_fva_comparison.loc[common_ids]

comparison_df = pd.DataFrame({
    "baseline_range_auto": auto_common["baseline_range"],
    "closed_range_auto": auto_common["closed_range"],
    "relative_change_auto": auto_common["relative_change"],
    "baseline_range_manual": manual_common["baseline_range"],
    "closed_range_manual": manual_common["closed_range"],
    "relative_change_manual": manual_common["relative_change"],
})

comparison_path = IN_DIR / "auto_manual_host_fva_comparison.csv"
comparison_df.to_csv(comparison_path)
