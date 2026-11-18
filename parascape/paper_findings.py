"""The analyses behind the paper, declared once so the program can run them.

Instruction 101, items B2 and B3: every analysis in the paper has to be runnable inside ParaScape
by somebody who did not write it, and has to save its result beside the record of what it was made
from. The analyses are corpus-level rather than gene-level -- they are about what the literature
can and cannot say -- so they belong in the Data band beside Datasets and Information, not in a
module that takes a gene.

ONE TABLE. The key, the figure it belongs to, the claim it ends on, the script that produces it
and where its result lands. Nothing else in the program may carry a second copy of this: the
module, the tests and any future report are all computed from it.

No Qt here, so the catalogue can be listed, tested and run headless.
"""
from __future__ import annotations

import glob
import json
import os
import subprocess
import sys
from dataclasses import dataclass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@dataclass(frozen=True)
class Finding:
    """One analysis of the paper.

    :param key: stable identifier.
    :param figure: which figure it is drawn as.
    :param name: what it is called, in a few words.
    :param claim: the sentence it ends on, about the parasite. Not a description of the method:
        a finding whose claim is about the data is a finding that has not finished.
    :param script: the script under `scripts/` that produces it.
    :param results: a glob under `results/` matching the directories it writes.
    :param needs_omics: whether it reads the omics store rather than the release alone.
    """

    key: str
    figure: str
    name: str
    claim: str
    script: str
    results: str
    needs_omics: bool = False


FINDINGS = (
    Finding("corpus", "1b", "What the corpus reaches",
            "Transcript reaches 99.9 per cent of genes and host response 5.3 per cent.",
            "build_corpus_inventory.py", "corpus_inventory_*"),
    Finding("comparability", "2", "Which studies can be compared",
            "Removing the cell-cycle axis raises agreement between untreated arms from 0.832 to "
            "0.887, three times what a random axis gives.",
            "build_control_comparability.py", "control_comparability_*", needs_omics=True),
    Finding("modules", "3", "The co-expression modules",
            "14 of 23 modules are a compartment, and the secretory machinery splits into two.",
            "build_coexpression_modules.py", "coexpression_modules_*"),
    Finding("structure", "4", "What structures the corpus",
            "Two contrasts from one laboratory agree 3.3 times better than two from different "
            "laboratories.",
            "build_what_structures_contrasts.py", "what_structures_contrasts_*"),
    Finding("reach", "5", "Which subtype a treatment reaches",
            "Bradyzoite subtype C is reached by 4 of 201 published contrasts, against 86 for E.",
            "build_bradyzoite_reach.py", "bradyzoite_reach_*"),
    Finding("protein", "6b", "What sets protein abundance",
            "Transcript explains 48.6 per cent of the variance and every other layer adds 1.8.",
            "build_what_sets_protein_abundance.py", "what_sets_protein_abundance_*"),
    Finding("explains", "6c", "Which data explains which",
            "Protein abundance explains fitness at R2 = 0.373: what the parasite makes most of is "
            "what it can least afford to lose.",
            "build_which_data_explains_which.py", "which_data_explains_which_*"),
    Finding("network", "7", "The whole network against one layer",
            "Every layer at once reaches 24 times the genes of the best single layer, at two "
            "thirds the accuracy.",
            "build_network_layers.py", "network_layers_*"),
    Finding("strains", "8", "What differs between strains",
            "Dense granules vary at 2.06 against a length-matched 1.36; everything internal is "
            "at or below it.",
            "build_strain_variation_by_class.py", "strain_variation_by_class_*"),
    Finding("species", "9", "What transfers between species",
            "Where a protein is transfers between Toxoplasma and Plasmodium at 0.630 against a "
            "0.374 baseline; how its message is written does not transfer at all (0.027).",
            "build_what_transfers_between_species.py", "species_transfer_*"),
    Finding("host", "10a", "The host response per effector",
            "3,435 host genes are moved by exactly one of 21 effectors and 68 by half or more.",
            "build_host_response_shared.py", "host_response_shared_*"),
    Finding("single_cell", "10b", "The parasite in the same cells",
            "One of 85 effector knockouts changes the parasite total output, so the pseudobulk "
            "assumption holds.",
            "build_dual_perturb_single_cell.py", "dual_perturb_single_cell_*"),
    Finding("filling", "11", "What the map can fill",
            "5 of 12 targets are called better than guessing; is_exported is worse than guessing.",
            "build_filling_the_map.py", "filling_the_map_*"),
    Finding("absence", "12", "What nothing measures",
            "A host-response assay would give 7,707 genes their first measurement of that kind.",
            "build_measured_absence.py", "measured_absence_*"),
)


def get(key: str) -> Finding:
    """One finding by key.

    :raises KeyError: for a key that is not a finding.
    """
    for finding in FINDINGS:
        if finding.key == key:
            return finding
    raise KeyError(f"no finding named {key!r}")


def script_path(finding: Finding) -> str:
    return os.path.join(ROOT, "scripts", finding.script)


def latest(finding: Finding) -> str:
    """The newest result directory this finding has written, or "" when it has never run."""
    found = sorted(glob.glob(os.path.join(ROOT, "results", finding.results)))
    return found[-1] if found else ""


def record(finding: Finding) -> dict:
    """What the newest run was made from: its summary, with the directory it came from.

    This is the B3 half of the contract. A result with no record of its inputs is a number
    somebody has to take on trust, and the whole point of the paper is that nobody should have to.
    """
    directory = latest(finding)
    if not directory:
        return {}
    path = os.path.join(directory, "summary.json")
    if not os.path.exists(path):
        # some findings write one summary per sample rather than one per run
        other = sorted(glob.glob(os.path.join(directory, "*summary.json")))
        if not other:
            return {"result_directory": directory}
        path = other[-1]
    try:
        out = json.load(open(path))
    except (OSError, ValueError) as error:
        return {"result_directory": directory, "unreadable": str(error)}
    out["result_directory"] = directory
    return out


def run(key: str, cache: str = "", omics: str = "", timeout: int = 7200):
    """Run one finding, in its own process, and return the completed process.

    A separate process and not an import: these analyses read gigabytes and one of them draws, so
    a failure inside one must not take the window with it.
    """
    finding = get(key)
    env = dict(os.environ)
    if cache:
        env["PARASCAPE_CACHE"] = cache
    if omics:
        env["PARASCAPE_OMICS"] = omics
    return subprocess.run([sys.executable, script_path(finding)], cwd=ROOT, env=env,
                          capture_output=True, text=True, timeout=timeout)
