#!/usr/bin/env python3
"""
Benchmark a MeMoMe matching table against a reference (e.g. manually curated)
matching table over a grid of InChI / DB / name cutoffs.

For every cutoff combination the MeMoMe table is filtered with
``filter_matching_table`` (the same filter main.py uses) and the resulting
pairs are compared to the reference pairs.

How a cutoff combination selects pairs (filter_matching_table):
    1. A candidate pair is kept if
           inchi_score >= InChI cutoff                              (rule A)
           OR (Name_score >= name cutoff AND DB_score >= DB cutoff) (rule B)
       inchi_score is always 0 or 1, so an InChI cutoff of 1.0 lets every InChI
       match through on its own, and any cutoff > 1 (default grid: 1.01) switches
       rule A off. InChI then only counts through total_score in step 2.
    2. Per met_id1 the kept candidate with the highest --score-type (default
       total_score = mean of inchi, DB, name and formula score) wins. Ties are
       broken by the ids: the alphabetically first met_id2 wins.
    As a consequence every met_id1 gets at most one partner, so a reference that
    maps one met_id1 to several met_id2 can never reach a recall of 1.

How reference ids are matched to MeMoMe ids:
    MeMoMe stores ids without prefix/compartment ('M_glc__D_e' -> 'glc__D').
    A reference id is used as-is if the MeMoMe table contains it; otherwise it
    is stripped the same way ('tre_c' -> 'tre'), and if that is still unknown
    the SBML escaping of special characters is applied ('ala-L' -> 'ala__45__L').
    The start of the output reports how many reference pairs were found in the
    MeMoMe table - pairs that are not found count as FN for every cutoff.

Scoring:
    TP  predicted pair that is in the reference
    FP  predicted pair that is not in the reference, but at least one of its
        metabolites appears in the reference (the curator matched it to
        something else). Predicted pairs where neither metabolite appears in
        the reference cannot be verified and are reported as "unverifiable";
        pass --count-unverifiable-as-fp to count them as FP instead.
    FN  reference pair that was not predicted

Exchange-only mode (--exchanges-only, needs --model1/--model2):
    Only metabolites that are exchanged with the environment are evaluated -
    these are the metabolites the ModelMerger translates between the models.
    The filter still runs on the full table (exactly as in main.py); afterwards
    both the predicted and the reference pairs are restricted to pairs whose
    two metabolites are exchange metabolites of their model.
    Exchange metabolites are the participants of reactions with a single
    (non-boundary) metabolite that sit in the external compartment, taken as the
    compartment holding most of these reactions (sinks/demands are excluded).
    Reference pairs outside the scope are dropped and listed in the output.

Picking the best cutoffs:
    All combinations are sorted by --metric (default f1). Ties are broken by the
    remaining two metrics, then by the lower InChI cutoff (rule A on, as in the
    main.py defaults), then by the higher DB and name cutoffs (the strictest
    combination that reaches the same result).

Example:
    python3 utils/benchmark_matching.py \
        --memome matching_tables/gapseq_recon3D_matching_table.csv \
        --reference tests/dat/manually_merged_models/gapseq_recon3D/met_matches.csv \
        --output benchmark_grid.csv

    # exchanges only
    python3 utils/benchmark_matching.py \
        --memome matching_tables/gapseq_recon3D_matching_table.csv \
        --reference tests/dat/manually_merged_models/gapseq_recon3D/met_matches.csv \
        --exchanges-only \
        --model1 tests/dat/manually_merged_models/gapseq_recon3D/M1_recon3D_301_modified.xml \
        --model2 tests/dat/manually_merged_models/gapseq_recon3D/M2_bacterial_model.xml
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import libsbml
import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix
from src.matchMets import filter_matching_table


def load_memome_table(path: Path, score_type: str) -> pd.DataFrame:
    """
    Load a matching table produced by MeMoModel.match().
    Parameters
    ----------
        path : Path - the matching_table.csv written by MeMoMe
        score_type : str - score column used by filter_matching_table to rank candidates
    Returns the table reduced to the columns needed for filtering, with missing scores set to 0.
    """
    columns = ["met_id1", "met_id2", "inchi_score", "DB_score", "Name_score", score_type]
    # only read the needed columns - full tables can be several hundred MB
    table = pd.read_csv(path, usecols=lambda c: c in columns)
    missing = set(columns) - set(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    # rows added by --keep-unmatched have only one id and can never be a pair
    table = table.dropna(subset=["met_id1", "met_id2"])
    # a missing score means that matching method found nothing for this pair
    return table.fillna({"inchi_score": 0, "DB_score": 0, "Name_score": 0, score_type: 0})


def _detect_column(columns: list[str], prefix: str) -> str:
    """
    Find the reference column holding the ids of one model. The met_matches.csv
    convention names them M1_<namespace> and M2_<namespace>, but the column order
    differs between the examples, so we look them up by prefix.
    """
    hits = [c for c in columns if c.startswith(prefix)]
    if len(hits) != 1:
        raise ValueError(f"Could not detect a unique '{prefix}*' column in reference ({hits}); "
                         f"set it explicitly with --ref-m1-col/--ref-m2-col")
    return hits[0]


def _resolve_id(met_id: str, known_ids: set[str]) -> str:
    """
    Translate a reference id into the id MeMoMe uses in its matching table.
    Keep the id if MeMoMe knows it, otherwise strip prefix/compartment the way MeMoMe does
    (e.g. 'tre_c' -> 'tre') and/or apply the SBML escaping of special characters
    (e.g. 'ala-L' -> 'ala__45__L'). Ids are only changed when needed, because stripping a
    valid id can corrupt it (e.g. an id that legitimately ends in '_<lowercase letter>').
    """
    met_id = str(met_id).strip()
    if met_id in known_ids:
        return met_id
    # handle_metabolites_prefix_suffix returns None for invalid cpd ids - keep the raw id then
    stripped = handle_metabolites_prefix_suffix(met_id) or met_id
    # SBML ids may only contain letters, digits and '_', so exporters (e.g. cobra) write other
    # characters as __<ascii code>__ - the reference usually has the readable form
    for candidate in (stripped, _sbml_escape(met_id), _sbml_escape(stripped)):
        if candidate in known_ids:
            return candidate
    return stripped


def _sbml_escape(met_id: str) -> str:
    """Escape characters that are not allowed in SBML ids the way cobra does ('-' -> '__45__')."""
    return re.sub(r"[^A-Za-z0-9_]", lambda m: f"__{ord(m.group())}__", met_id)


def exchange_metabolite_ids(sbml_path: Path) -> set[str]:
    """
    Return the MeMoMe ids of all metabolites that are exchanged with the environment.
    Parameters
    ----------
        sbml_path : Path - the SBML model the matching table was created from
    The file is read with libsbml (not cobra) and ids go through
    handle_metabolites_prefix_suffix, exactly as in parseMetaboliteInfoFromSBML - cobra would
    decode ids like 'ala__45__L' to 'ala-L', which would not match the matching table.
    """
    if not sbml_path.exists():
        raise FileNotFoundError(f"{sbml_path} does not exist")
    model = libsbml.SBMLReader().readSBMLFromFile(str(sbml_path)).getModel()
    # species with boundaryCondition=True are old-style boundary metabolites (e.g. 'glc_b'),
    # they are not part of the network and are ignored when counting reaction participants
    boundary_species = {sp.getId() for sp in model.getListOfSpecies() if sp.getBoundaryCondition()}
    compartment_of = {sp.getId(): sp.getCompartment() for sp in model.getListOfSpecies()}

    # exchange reactions have exactly one (non-boundary) participant, e.g. 'glc_e <=>'
    exchanged = []
    for rxn in model.getListOfReactions():
        participants = {ref.getSpecies() for ref in list(rxn.getListOfReactants()) + list(rxn.getListOfProducts())}
        participants -= boundary_species
        if len(participants) == 1:
            exchanged.append(participants.pop())

    # sinks and demands have the same form but sit in internal compartments; like
    # cobra.medium.find_external_compartment, take the compartment holding most of these
    # reactions as the external one and keep only its metabolites
    if not exchanged:
        return set()
    external = Counter(compartment_of[sp] for sp in exchanged).most_common(1)[0][0]
    ids = {handle_metabolites_prefix_suffix(sp) for sp in exchanged if compartment_of[sp] == external}
    # handle_metabolites_prefix_suffix returns None for invalid cpd ids
    ids.discard(None)
    return ids


def load_reference(path: Path, memome: pd.DataFrame,
                   m1_col: str | None = None, m2_col: str | None = None) -> set[tuple[str, str]]:
    """
    Load the reference matching table (e.g. a manually curated met_matches.csv).
    Parameters
    ----------
        path : Path - the reference csv
        memome : pd.DataFrame - the MeMoMe matching table, used to translate the reference ids
        m1_col, m2_col : str|None - reference columns with the model1/model2 ids, detected if None
    Returns the reference as a set of (met_id1, met_id2) pairs in MeMoMe's id convention.
    """
    # read as str so numeric-looking ids are not converted
    ref = pd.read_csv(path, dtype=str)
    m1_col = m1_col or _detect_column(list(ref.columns), "M1_")
    m2_col = m2_col or _detect_column(list(ref.columns), "M2_")
    ref = ref[[m1_col, m2_col]].dropna()
    ids1 = set(memome["met_id1"])
    ids2 = set(memome["met_id2"])
    return {(_resolve_id(a, ids1), _resolve_id(b, ids2)) for a, b in zip(ref[m1_col], ref[m2_col])}


def score(predicted: set[tuple[str, str]], reference: set[tuple[str, str]],
          count_unverifiable_as_fp: bool = False) -> dict:
    """
    Compare the predicted pairs to the reference pairs (see the module docstring for
    the TP/FP/FN definitions) and return counts plus precision, recall and F1.
    """
    ref_ids1 = {a for a, _ in reference}
    ref_ids2 = {b for _, b in reference}
    tp = len(predicted & reference)
    wrong = predicted - reference
    # a wrong pair is only a verifiable error if the curator looked at one of its
    # metabolites - the reference usually covers only the shared exchange metabolites
    verifiable = {p for p in wrong if p[0] in ref_ids1 or p[1] in ref_ids2}
    unverifiable = len(wrong) - len(verifiable)
    fp = len(wrong) if count_unverifiable_as_fp else len(verifiable)
    fn = len(reference - predicted)
    # guard against division by zero for cutoffs that select nothing
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {"predicted": len(predicted), "TP": tp, "FP": fp, "FN": fn,
            "unverifiable": unverifiable, "precision": precision, "recall": recall, "f1": f1}


def run_grid(memome: pd.DataFrame, reference: set[tuple[str, str]],
             inchi_cutoffs: list[float], db_cutoffs: list[float], name_cutoffs: list[float],
             score_type: str = "total_score", count_unverifiable_as_fp: bool = False,
             scope: tuple[set[str], set[str]] | None = None) -> pd.DataFrame:
    """
    Filter the MeMoMe table with every combination of the given cutoffs and score each result.
    Parameters
    ----------
        scope : tuple of (model1 ids, model2 ids) or None - if given, only predicted pairs whose
            met_id1 and met_id2 are in these sets are scored (used for --exchanges-only). The
            reference is expected to be restricted to the same scope already.
    Returns one row per cutoff combination with the cutoffs and the output of score().
    """
    # rows that fail even the loosest cutoffs can never be selected - drop them once
    loosest = ((memome["inchi_score"] >= min(inchi_cutoffs))
               | ((memome["Name_score"] >= min(name_cutoffs)) & (memome["DB_score"] >= min(db_cutoffs))))
    candidates = memome.loc[loosest]

    # full cartesian product of the cutoffs - filtering the pruned table is fast,
    # so a dense grid is fine
    rows = []
    for inchi in inchi_cutoffs:
        for db in db_cutoffs:
            for name in name_cutoffs:
                filtered = filter_matching_table(candidates, Inchi_threshold=inchi, DB_threshold=db,
                                                 Name_threshold=name, score_type=score_type)
                # filter_matching_table keeps at most one partner per met_id1
                predicted = set(zip(filtered["met_id1"], filtered["met_id2"]))
                # restrict only after filtering, so the chosen partner per metabolite is the
                # same as in a real main.py run
                if scope is not None:
                    predicted = {(a, b) for a, b in predicted if a in scope[0] and b in scope[1]}
                rows.append({"InChI_threshold": inchi, "DB_threshold": db, "Name_threshold": name,
                             **score(predicted, reference, count_unverifiable_as_fp)})
    return pd.DataFrame(rows)


def cutoff_range(spec: str) -> list[float]:
    """Parse 'start:stop:step' (inclusive) or a comma separated list."""
    if ":" in spec:
        start, stop, step = (float(x) for x in spec.split(":"))
        # add half a step so 'stop' is included despite float rounding, and round
        # so values like 0.30000000000000004 print and compare cleanly
        return [round(x, 6) for x in np.arange(start, stop + step / 2, step)]
    return [float(x) for x in spec.split(",")]


def main() -> None:
    """Parse the command line, run the grid and print/save the results."""
    # show the module docstring in --help, it explains the selection and scoring rules
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--memome", required=True, type=Path, help="Matching table generated by MeMoMe (matching_table.csv)")
    parser.add_argument("--reference", required=True, type=Path, help="Reference matching table (e.g. met_matches.csv)")
    parser.add_argument("--ref-m1-col", default=None, help="Reference column with model1 ids (default: the column starting with 'M1_')")
    parser.add_argument("--ref-m2-col", default=None, help="Reference column with model2 ids (default: the column starting with 'M2_')")
    parser.add_argument("--inchi", default="1.0,1.01", help="InChI cutoffs; the score is 0/1, so 1.01 disables InChI-only matches (default: 1.0,1.01)")
    parser.add_argument("--db", default="0:1:0.05", help="DB cutoffs as start:stop:step or comma list (default: 0:1:0.05)")
    parser.add_argument("--name", default="0.5:1:0.05", help="Name cutoffs as start:stop:step or comma list (default: 0.5:1:0.05)")
    parser.add_argument("--score-type", default="total_score", help="Score column used to pick the best candidate per metabolite (default: total_score)")
    parser.add_argument("--metric", default="f1", choices=["f1", "precision", "recall"], help="Metric used to pick the best cutoffs (default: f1)")
    parser.add_argument("--count-unverifiable-as-fp", action="store_true", help="Count predicted pairs whose metabolites are absent from the reference as FP")
    parser.add_argument("--exchanges-only", action="store_true", help="Only evaluate metabolites that are exchanged with the environment in both models (needs --model1/--model2)")
    parser.add_argument("--model1", type=Path, default=None, help="SBML of model1 the matching table was created from (for --exchanges-only)")
    parser.add_argument("--model2", type=Path, default=None, help="SBML of model2 the matching table was created from (for --exchanges-only)")
    parser.add_argument("--top", default=10, type=int, help="Number of best cutoff combinations to print (default: 10)")
    parser.add_argument("--output", default=None, type=Path, help="Optional CSV to write the full grid results to")
    args = parser.parse_args()
    if args.exchanges_only and (args.model1 is None or args.model2 is None):
        parser.error("--exchanges-only needs --model1 and --model2")

    memome = load_memome_table(args.memome, args.score_type)
    reference = load_reference(args.reference, memome, args.ref_m1_col, args.ref_m2_col)
    # reference pairs absent from the MeMoMe table are FN for every cutoff - report how many
    in_table = set(zip(memome["met_id1"], memome["met_id2"]))
    print(f"MeMoMe table: {len(memome)} rows | reference: {len(reference)} pairs, "
          f"{len(reference & in_table)} of them present in the MeMoMe table")

    scope = None
    if args.exchanges_only:
        scope = (exchange_metabolite_ids(args.model1), exchange_metabolite_ids(args.model2))
        # reference pairs outside the scope could never be predicted in scope - drop them
        # instead of counting them as FN, and say which ones so they can be checked
        out_of_scope = sorted(p for p in reference if p[0] not in scope[0] or p[1] not in scope[1])
        reference = reference - set(out_of_scope)
        print(f"Exchanges only: {len(scope[0])} exchange metabolites in model1, {len(scope[1])} in model2 | "
              f"{len(reference)} reference pairs in scope")
        if out_of_scope:
            print(f"  {len(out_of_scope)} reference pairs dropped (not exchanged in both models): {out_of_scope}")

    results = run_grid(memome, reference, cutoff_range(args.inchi), cutoff_range(args.db),
                       cutoff_range(args.name), args.score_type, args.count_unverifiable_as_fp, scope)
    # ties: prefer the other metrics, then keep the InChI-only rule enabled (lower InChI
    # cutoff, as in main.py's defaults), then prefer the stricter DB and name cutoffs
    tie_breakers = [m for m in ["f1", "precision", "recall"] if m != args.metric]
    results = results.sort_values(
        by=[args.metric, *tie_breakers, "InChI_threshold", "DB_threshold", "Name_threshold"],
        ascending=[False, False, False, True, False, False],
    ).reset_index(drop=True)

    if args.output:
        results.to_csv(args.output, index=False)
        print(f"Wrote {len(results)} grid results to {args.output}")

    with pd.option_context("display.width", 200, "display.max_columns", None, "display.float_format", "{:.3f}".format):
        print(f"\nTop {args.top} by {args.metric}:")
        print(results.head(args.top).to_string(index=False))
    # print the winner as flags that can be pasted into main.py
    best = results.iloc[0]
    print(f"\nBest: --InChI-threshold {best.InChI_threshold} --DB-threshold {best.DB_threshold} "
          f"--Name-threshold {best.Name_threshold}  ->  {args.metric}={best[args.metric]:.3f} "
          f"(precision={best.precision:.3f}, recall={best.recall:.3f})")


if __name__ == "__main__":
    main()
