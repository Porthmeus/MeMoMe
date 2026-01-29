#!/usr/bin/env python3
"""Over-representation analysis (ORA) for subsystem enrichment."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_DROP = ("Transport, extracellular", "Exchange/demand reaction")


def _load_dependencies(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "reaction_id" not in df.columns:
        first = df.columns[0]
        df = df.rename(columns={first: "reaction_id"})
    df["reaction_id"] = df["reaction_id"].astype(str)
    if "subsystem" not in df.columns:
        raise ValueError("Missing subsystem column in dependencies file.")
    df["subsystem"] = df["subsystem"].fillna("NA").astype(str)
    return df


def _load_list(path: Path) -> set[str]:
    df = pd.read_csv(path)
    col = df.columns[0]
    return set(df[col].astype(str))


def _benjamini_hochberg(pvals: list[float]) -> list[float]:
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    qvals = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = n - rank
        q = min(prev, pvals[idx] * n / (i + 1))
        qvals[idx] = q
        prev = q
    return qvals


def _fisher_exact(overlap: int, in_list: int, in_bg: int, total: int) -> float:
    # Right-tailed Fisher exact using hypergeometric survival function.
    try:
        from scipy.stats import hypergeom
    except ImportError as exc:
        raise ImportError("Install: pip install scipy") from exc
    return float(hypergeom.sf(overlap - 1, total, in_bg, in_list))


def _ora(
    list_ids: set[str],
    background: pd.DataFrame,
    out_path: Path,
) -> None:
    background = background.copy()
    background["in_list"] = background["reaction_id"].isin(list_ids)
    total = len(background)
    in_list = int(background["in_list"].sum())

    rows = []
    for subsystem, group in background.groupby("subsystem"):
        in_bg = len(group)
        overlap = int(group["in_list"].sum())
        pval = _fisher_exact(overlap, in_list, in_bg, total) if in_list > 0 else 1.0
        rows.append(
            {
                "subsystem": subsystem,
                "overlap": overlap,
                "list_size": in_list,
                "subsystem_size": in_bg,
                "background_size": total,
                "p_value": pval,
            }
        )
    res = pd.DataFrame(rows)
    res["fdr_bh"] = _benjamini_hochberg(res["p_value"].tolist())
    res = res.sort_values(["fdr_bh", "p_value", "overlap"], ascending=[True, True, False])
    res.to_csv(out_path, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="ORA subsystem enrichment using high-dependency lists.")
    parser.add_argument(
        "--dependencies",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare/"
            "microbiome_dependency/all_dependencies_with_subsystems.csv"
        ),
    )
    parser.add_argument(
        "--auto-list",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare/"
            "microbiome_dependency/auto_high_dependency.csv"
        ),
    )
    parser.add_argument(
        "--manual-list",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare/"
            "microbiome_dependency/manual_high_dependency.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare/"
            "microbiome_dependency/ora_subsystems"
        ),
    )
    parser.add_argument(
        "--suffix",
        type=str,
        default="",
        help="Optional suffix to append to output filenames.",
    )
    parser.add_argument(
        "--drop-subsystems",
        nargs="*",
        default=list(DEFAULT_DROP),
        help="Subsystems to exclude from analysis.",
    )
    parser.add_argument(
        "--loop-mode",
        choices=["include", "exclude"],
        default="exclude",
        help=(
            "include: add looping reactions to the high-dependency lists; "
            "exclude: remove looping reactions from both background and lists."
        ),
    )
    args = parser.parse_args()

    deps = _load_dependencies(args.dependencies)
    if args.drop_subsystems:
        deps = deps[~deps["subsystem"].isin(set(args.drop_subsystems))]

    auto_ids = _load_list(args.auto_list)
    manual_ids = _load_list(args.manual_list)

    if "part_of_loop" in deps.columns:
        loop_ids = set(deps.loc[deps["part_of_loop"].astype(bool), "reaction_id"])
        if args.loop_mode == "exclude":
            deps = deps[~deps["part_of_loop"].astype(bool)]
            auto_ids = auto_ids - loop_ids
            manual_ids = manual_ids - loop_ids
        else:
            rel_cols = {
                "auto": "auto_rel_change",
                "manual": "manual_rel_change",
            }
            rel_map = deps.set_index("reaction_id")[list(rel_cols.values())]
            auto_loop = set(rel_map[rel_map[rel_cols["auto"]] <= -0.5].index) & loop_ids
            manual_loop = set(rel_map[rel_map[rel_cols["manual"]] <= -0.5].index) & loop_ids
            auto_ids = auto_ids | auto_loop
            manual_ids = manual_ids | manual_loop

    args.output_dir.mkdir(parents=True, exist_ok=True)
    suffix = f"_{args.suffix}" if args.suffix else ""
    _ora(auto_ids, deps, args.output_dir / f"auto_ora_results{suffix}.csv")
    _ora(manual_ids, deps, args.output_dir / f"manual_ora_results{suffix}.csv")


if __name__ == "__main__":
    main()
