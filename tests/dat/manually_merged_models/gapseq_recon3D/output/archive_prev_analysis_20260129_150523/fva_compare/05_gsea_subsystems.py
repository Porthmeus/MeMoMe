#!/usr/bin/env python3
"""Run subsystem enrichment (preranked GSEA) on microbiome-dependency scores."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


DEFAULT_DROP = ("Transport, extracellular", "Exchange/demand reaction")


def _load_dependencies(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "reaction_id" not in df.columns:
        first = df.columns[0]
        if first.lower().startswith("unnamed"):
            df = df.rename(columns={first: "reaction_id"})
        else:
            df = df.rename(columns={first: "reaction_id"})
    df["reaction_id"] = df["reaction_id"].astype(str)
    return df


def _load_mapping(df: pd.DataFrame, drop_subsystems: tuple[str, ...]) -> pd.Series:
    if "subsystem" not in df.columns:
        return pd.Series(dtype=str)
    mapping = df[["reaction_id", "subsystem"]].copy()
    mapping["subsystem"] = mapping["subsystem"].fillna("NA").astype(str)
    if drop_subsystems:
        mapping = mapping[~mapping["subsystem"].isin(set(drop_subsystems))]
    return mapping.set_index("reaction_id")["subsystem"]


def _run_prerank(
    scores: pd.Series,
    subsystem_map: pd.Series,
    outdir: Path,
    min_size: int,
    max_size: int,
    permutation_num: int,
    seed: int,
) -> pd.DataFrame:
    try:
        import gseapy as gp
    except ImportError as exc:
        raise ImportError("Install: pip install gseapy") from exc

    used = [rid for rid in scores.index if rid in subsystem_map.index]
    if len(used) < 2000:
        raise ValueError(
            f"Too few reactions after intersection (n={len(used)}). "
            "Check reaction IDs / ensure you exported ALL reactions."
        )
    rnk = scores.loc[used].sort_index().sort_values(ascending=False)

    subs: dict[str, list[str]] = defaultdict(list)
    for rid in rnk.index:
        subs[subsystem_map[rid]].append(rid)

    gene_sets = {k: v for k, v in subs.items() if min_size <= len(v) <= max_size}
    if not gene_sets:
        raise ValueError("No subsystems passed size filters.")

    outdir.mkdir(parents=True, exist_ok=True)
    pre = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        outdir=str(outdir),
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,
        seed=seed,
        verbose=False,
    )
    res = pre.res2d.sort_values(["FDR q-val", "NES"], ascending=[True, False])
    return res


def main() -> None:
    parser = argparse.ArgumentParser(description="Subsystem GSEA on dependency scores.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare/"
            "microbiome_dependency/all_dependencies_with_subsystems.csv"
        ),
        help="CSV with reaction_id, auto_rel_change, manual_rel_change, subsystem, part_of_loop.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare/"
            "microbiome_dependency/gsea_subsystems"
        ),
    )
    parser.add_argument("--min-size", type=int, default=10)
    parser.add_argument("--max-size", type=int, default=3000)
    parser.add_argument("--permutation-num", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument(
        "--include-loops",
        action="store_true",
        help="Include loop reactions; default excludes part_of_loop==True.",
    )
    parser.add_argument(
        "--drop-subsystems",
        nargs="*",
        default=list(DEFAULT_DROP),
        help="Subsystems to exclude from analysis.",
    )
    args = parser.parse_args()

    df = _load_dependencies(args.input)
    if "part_of_loop" in df.columns and not args.include_loops:
        df = df[~df["part_of_loop"].astype(bool)]

    subsystem_map = _load_mapping(df, tuple(args.drop_subsystems))
    if subsystem_map.empty:
        raise ValueError("No subsystem mapping found in input.")

    # Define dependency score as -relative_change (so higher = more dependent).
    for label, col in (("auto", "auto_rel_change"), ("manual", "manual_rel_change")):
        if col not in df.columns:
            raise ValueError(f"Missing column: {col}")
        scores = (-df.set_index("reaction_id")[col].astype(float)).dropna()
        outdir = args.output_dir / f"{label}_prerank"
        res = _run_prerank(
            scores=scores,
            subsystem_map=subsystem_map,
            outdir=outdir,
            min_size=args.min_size,
            max_size=args.max_size,
            permutation_num=args.permutation_num,
            seed=args.seed,
        )
        res.to_csv(args.output_dir / f"{label}_gsea_results.csv")


if __name__ == "__main__":
    main()
