"""Instruction 65.04: the biological eligibility read.

The selection rule in instruction 65 ranks datasets by citation rate, newness and
comprehensiveness -- but it says plainly that **biological admission precedes ranking**:

    organism, host/cell type, stage, assay, quantity, units, measured versus inferred
    evidence, mapping, quality checks and data access must fit the question.

Ranking a source that does not answer the question is worse than not ranking at all, because the
number it produces looks like a comparison. This script performs the admission step for every
candidate that currently fills a slot, and it does it by comparing **two declarations** rather than
by my reading of a paper:

  * what the SLOT asks -- ``Slot.axis`` (the quantity), ``Slot.unit`` (what a row is),
    ``Slot.context`` and the three facet paths (``evidence_path``, ``biology_path``,
    ``context_path``), and ``Slot.role``;
  * what the SOURCE declares -- ``organism``, ``kind`` (the assay), ``level`` (the quantity class),
    ``coverage``, ``accession``, ``columns``.

Where the two declarations disagree, that is a finding regardless of what any paper says. Where the
declarations are silent, the verdict is **unknown** and the row is queued for reading -- never
guessed, never counted as a pass. That is the rule instruction 65 states for dates and counts and it
applies with more force to biology.

Offline and deterministic: it reads the shipped registry, the slot declarations and the audit tables
already in ``results/dataset_selection_2026_11_08/``. It makes no network call, so it can be re-run
to reproduce the table exactly.

Writes ``results/dataset_selection_2026_11_08/biological_eligibility.tsv``.
"""
from __future__ import annotations

import json
import os
import re
import sys
from collections import Counter, defaultdict

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# A script run by path puts its OWN directory on sys.path, not the working directory, so without
# this it would import whichever parascape is installed rather than the checkout being audited.
sys.path.insert(0, ROOT)

from parascape import datasets, slots as S  # noqa: E402

OUT = os.path.join(ROOT, "results", "dataset_selection_2026_11_08")
CANDIDATES = os.path.join(OUT, "slot_candidates_scored.tsv")


# --------------------------------------------------------------------------- assay vocabulary

def normalise_kind(kind: str) -> str:
    """Fold the spellings of one assay together.

    The registry spells two assays more than one way -- ``IP-MS``/``IPMS`` and
    ``Ribo-seq``/``RiboSeq``/``riboseq`` -- across 13 entries. Anything that groups or filters by
    assay without folding them treats one assay as two or three, which is how an audit silently
    reports a comparison it never made. Folding here, and reporting the drift separately, keeps the
    registry readable while the audit stays correct.
    """
    return re.sub(r"[^a-z0-9]", "", kind.lower())


# What each assay can actually measure. Stated explicitly rather than inferred from current usage,
# because inferring it from what is already installed would make every incumbent eligible by
# construction -- the circularity instruction 65 exists to prevent.
#
# A first version of this table was too NARROW, and it reported 28 candidates as biologically
# ineligible when every one of them was a correct pairing and the table was wrong. Each widening
# below therefore names the source that justified it, so the next reader can check the claim rather
# than trust it. The lesson is recorded here because a too-strict admission table produces exactly
# the failure instruction 65 is trying to prevent, only in the opposite direction: it rejects good
# evidence with a confident verdict.
MEASURES = {
    # An RNA-seq experiment is not confined to steady-state abundance. Decay after transcriptional
    # arrest is regulation (`mrna_decay_gse329845`, actinomycin D), and a reannotation built on long
    # reads gives the architecture of the untranslated region (`gse302108_utr5`), which is sequence.
    "rnaseq": {"transcription", "regulation", "sequence"},
    "scrnaseq": {"transcription"},
    "microarray": {"transcription"},
    "longread": {"transcription", "sequence", "regulation"},
    # Long-read models classify splicing, and which isoform is used is regulation (`pf_isoforms`).
    "nanopore": {"transcription", "regulation", "sequence"},
    "transcription": {"transcription"},
    "ripseq": {"transcription", "regulation"},
    "merip": {"transcription", "PTM"},
    # Ribosome profiling ships the footprint channel AND the matched RNA channel; the channel is
    # chosen by the column, not by the assay label, which `_channel_axis` below handles.
    "riboseq": {"translation"},
    # Thermal proteome profiling is a proteomics readout of compound binding: target engagement is
    # chemistry (`pf_target_engagement`).
    "proteomics": {"protein abundance", "translation", "chemistry"},
    "lopit": {"localization"},
    "rdeep": {"localization", "relation"},
    "phosphoproteomics": {"PTM"},
    "acetylation": {"PTM"},
    "lactylation": {"PTM"},
    "myristoylation": {"PTM"},
    "palmitoylation": {"PTM"},
    "ipms": {"relation"},
    # A bait anchored in an outer membrane reports what is NEAR it, which is localization evidence
    # about the outside of an organelle (`organelle_surface_turboid`).
    "bioid": {"relation", "localization"},
    "proximitylabelling": {"relation", "localization"},
    "crosslinkms": {"relation"},
    "xlms": {"relation"},
    "crisprscreen": {"fitness"},
    "insertionscreen": {"fitness"},
    # Clones evolved against a compound and sequenced whole-genome give resistance alleles, which
    # are sequence as much as they are chemistry (`pf_resistome`).
    "invitroevolution": {"fitness", "chemistry", "sequence"},
    "imagingscreen": {"phenotype", "localization"},
    "chipseq": {"regulation"},
    "chipchip": {"regulation"},
    "cutandtag": {"regulation"},
    "atacseq": {"regulation"},
    "metabolomics": {"metabolism"},
    "lipidomics": {"metabolism"},
    "immunity": {"immunity"},
    "variation": {"sequence"},
    "sequence": {"sequence"},
    "domains": {"sequence"},
    "structure": {"sequence", "localization"},
    "computedstructure": {"sequence", "localization"},
    "computedsequencefeatures": {"sequence", "localization"},
    "orthology": {"relation", "sequence"},
    "graph": {"relation"},
    "identity": {"sequence"},
    "annotation": {"sequence", "localization", "relation"},
    # Curated readings of papers: the only route to a resistance allele in Toxoplasma, which has no
    # resistome to parse, and to per-gene drug sensitivity, for which no genome-wide screen exists.
    # Admitted as evidence, always marked INFERRED rather than measured.
    "literature": {"relation", "chemistry", "fitness"},
}

# A ribosome-profiling deposit carries two channels. The RNA channel answers a transcription slot
# and the footprint channel answers a translation slot, and the registry tells them apart by the
# column name rather than by the assay. Five candidates sit on this distinction -- `rna99395_`,
# `rna129869_`, `rna245775_parent_tachy_`, `rna245775_parent_prebrady_` and `riboseq_mrna_` -- and
# reading the assay alone called all five a mismatch.
RNA_CHANNEL = re.compile(r"(^|_)(rna|mrna)", re.I)


def _channel_axis(pattern: str) -> str:
    """The axis a COLUMN answers, where the assay carries more than one channel."""
    return "transcription" if RNA_CHANNEL.search(str(pattern or "")) else ""

# Evidence that is derived or asserted rather than measured on the organism. Instruction 65 requires
# measured and inferred evidence to be told apart; a derived source may still be eligible, but it
# must never be admitted as a measurement.
INFERRED = {"computedsequencefeatures", "computedstructure", "annotation", "orthology",
            "literature", "identity", "graph", "domains", "sequence", "structure"}

HOST_ORGANISMS = {"Hs", "Mm"}

# The quantity class the registry records, against the axes a slot can ask for. Two levels are
# deliberately permissive: `reference` describes a resource rather than a quantity, and
# `post_translation` covers both modifications and the abundances measured alongside them.
#
# `level` records where in the central dogma the deposit was MADE, which is not always the question
# it answers: a whole-genome resistome is DNA-level and answers a chemistry question, and a thermal
# shift is post-translational and answers one too. The widenings here follow the same verified
# examples as the assay table above.
LEVEL_AXES = {
    "DNA": {"regulation", "sequence", "fitness", "chemistry", "host effect"},
    "transcription": {"transcription", "regulation", "sequence"},
    "translation": {"translation", "protein abundance", "transcription"},
    "post_translation": {"PTM", "protein abundance", "localization", "relation", "translation",
                         "chemistry"},
    "reference": None,   # None = imposes no constraint; the assay check still applies
}

# Slots that are never used as features, and the metadata slots, are excluded from the model by
# construction. Asking whether the quantity of a source fits an axis that is literally called
# "NEVER a feature" is a category error, not a finding.
EXEMPT_AXES = {"NEVER a feature"}
EXEMPT_ROLES = {"never", "metadata"}


def stage_words(slot) -> set:
    """The life stage the slot asks about, from the facets rather than from the slot name."""
    words = set()
    for path in (slot.evidence_path, slot.biology_path, slot.context_path):
        for part in path:
            p = str(part).lower()
            for stage in ("tachyzoite", "bradyzoite", "sporozoite", "oocyst", "merozoite",
                          "gametocyte", "ring", "trophozoite", "schizont", "liver", "cyst"):
                if stage in p:
                    words.add(stage)
    return words


# A column name is a declaration of what was measured, and it is usually abbreviated: `invivo_TZ_`,
# `rna51780_tachy_`, `expr_cyst`, `expr_sporulated`. Matching the full word "tachyzoite" against
# those finds nothing and reports "the source declares no stage", which is false -- the source
# declares it in the only place it had room to.
#
# `cyst` is the trap: "oocyst" contains it, so a bradyzoite stem that matched `cyst` would claim
# every oocyst column as a tissue cyst. Each stem therefore carries what must NOT be around it.
STAGE_STEMS = {
    "tachyzoite": ([r"tachy", r"\btz\b", r"_tz_", r"_tz$"], []),
    "bradyzoite": ([r"brady", r"\bbz\b", r"_bz_", r"tissue[_ ]?cyst", r"(?<!oo)cyst"], []),
    "sporozoite": ([r"sporo", r"\bspz\b", r"sporulated"], []),
    "oocyst": ([r"oocyst", r"sporulated", r"unsporulated"], []),
    "merozoite": ([r"merozoite", r"\bmero"], []),
    "gametocyte": ([r"gametocyte", r"\bgam\b", r"gameto"], []),
    "ring": ([r"\bring\b"], []),
    "trophozoite": ([r"troph"], []),
    "schizont": ([r"schizont"], []),
    "liver": ([r"liver", r"hepatocyte"], []),
    "cyst": ([r"tissue[_ ]?cyst", r"(?<!oo)cyst"], []),
}


def _matches(stage: str, text: str) -> bool:
    stems, forbidden = STAGE_STEMS.get(stage, ([re.escape(stage)], []))
    if any(re.search(f, text) for f in forbidden):
        return False
    return any(re.search(s, text) for s in stems)


def declares_stage(entry, stage: str, columns: tuple, pattern: str = "") -> str:
    """How strongly the source declares this stage: 'column', 'prose' or ''.

    The distinction is not pedantry, and one row proves it. `gse245775` matched "bradyzoite" from
    its note -- but the note exists to say the opposite: the arms are 48 h in RPMI pH 8.3 that the
    submitters label `cell type: pre-bradyzoites`, "named for that and NOT for bradyzoite". A word
    search over prose found the word inside the sentence written to deny it.

    So a COLUMN name is treated as a declaration of what was measured, because that is what the
    column is; free prose is treated as a mention that has to be read by a person before it counts.
    Reporting a mention as a fit is how an audit passes a source it never checked.
    """
    # ONLY the column this candidate actually uses. Searching every column of the deposit let a
    # candidate inherit a stage from a SIBLING: `rna108740_Day` was credited with "bradyzoite"
    # because `rna108740_Tissue_cysts` sits in the same source, and `rpf245775_eif12ko_` was
    # credited with it because a `prebrady` arm does. A neighbouring column is not evidence about
    # this one -- that is the whole reason slots are filled by pattern rather than by deposit.
    if _matches(stage, str(pattern or "").lower()):
        return "column"
    prose = " ".join(str(x) for x in (entry.name, entry.note, entry.coverage, entry.key)).lower()
    prose += " " + " ".join(str(c).lower() for c in columns)
    return "prose" if _matches(stage, prose) else ""


def check(slot, entry, columns: tuple, pattern: str = "") -> dict:
    """Run every admission criterion for one slot/source pair, each with its evidence."""
    kind = normalise_kind(entry.kind)
    verdicts, evidence = {}, {}
    exempt = slot.axis in EXEMPT_AXES or slot.role in EXEMPT_ROLES
    # The column decides which channel of a multi-channel deposit is being used.
    channel = _channel_axis(pattern)

    # -- organism, including the host case, which is a fit rather than a mismatch ---------------
    if slot.unit == "host_gene":
        ok = entry.organism in HOST_ORGANISMS
        verdicts["organism"] = "fits" if ok else "mismatch"
        evidence["organism"] = (f"slot asks for host genes; source organism {entry.organism}"
                                + ("" if ok else " is the parasite, not the host"))
    else:
        ok = entry.organism == slot.organism
        verdicts["organism"] = "fits" if ok else "mismatch"
        evidence["organism"] = f"slot {slot.organism} · source {entry.organism}"

    # -- quantity: the level the registry records against the axis the slot asks for ------------
    allowed = LEVEL_AXES.get(entry.level, set())
    if exempt:
        verdicts["quantity"] = "not applicable"
        evidence["quantity"] = f"slot role '{slot.role}', axis '{slot.axis}': never a feature"
    elif allowed is None:
        # `level: reference` says the source is a resource rather than a point on the central
        # dogma, so the level field does not speak to the quantity -- but the ASSAY check does, and
        # for all 28 such candidates it passes. Calling these "unknown" queued 28 rows for a reading
        # with nothing specific to read, which is the false-ineligible mistake pointing the other
        # way: a queue padded with rows that need no work hides the rows that do.
        can = MEASURES.get(kind)
        if can is not None and (slot.axis in can or slot.axis == "host effect"
                                or (channel and channel == slot.axis)):
            verdicts["quantity"] = "carried by the assay"
            evidence["quantity"] = (f"level 'reference' states a resource, not a quantity; the "
                                    f"quantity is established by the assay, and '{entry.kind}' "
                                    f"does measure '{slot.axis}'")
        else:
            verdicts["quantity"] = "unknown"
            evidence["quantity"] = ("level 'reference' states a resource, not a quantity, and the "
                                    "assay does not establish it either")
    else:
        ok = slot.axis in allowed or (channel and channel == slot.axis)
        verdicts["quantity"] = "fits" if ok else "mismatch"
        evidence["quantity"] = (f"slot axis '{slot.axis}' · source level '{entry.level}'"
                                + (f" · column is the {channel} channel" if channel else ""))

    # -- assay: can this measurement produce that quantity at all? ------------------------------
    can = MEASURES.get(kind)
    if exempt:
        verdicts["assay"] = "not applicable"
        evidence["assay"] = f"slot role '{slot.role}': excluded from features by declaration"
    elif can is None:
        verdicts["assay"] = "unknown"
        evidence["assay"] = f"assay '{entry.kind}' is not in the capability table"
    else:
        ok = (slot.axis in can or slot.axis == "host effect"
              or (channel and channel == slot.axis))
        verdicts["assay"] = "fits" if ok else "mismatch"
        evidence["assay"] = (f"'{entry.kind}' measures {sorted(can)}; slot asks '{slot.axis}'"
                             + (f"; the column is the {channel} channel" if channel else ""))

    # -- measured against inferred ---------------------------------------------------------------
    inferred = kind in INFERRED
    verdicts["evidence_grade"] = "inferred" if inferred else "measured"
    evidence["evidence_grade"] = (f"'{entry.kind}' is derived or asserted, not measured on the "
                                  f"organism" if inferred else f"'{entry.kind}' is a measurement")

    # -- unit: what one row is -------------------------------------------------------------------
    if slot.unit in ("gene", "host_gene"):
        verdicts["unit"] = "fits"
        evidence["unit"] = f"slot rows are {slot.unit}; the registry is gene-keyed"
    elif slot.unit == "pair":
        ok = slot.axis == "relation" or kind in ("ipms", "bioid", "xlms", "crosslinkms",
                                                 "proximitylabelling", "graph", "orthology",
                                                 "literature", "rdeep")
        verdicts["unit"] = "fits" if ok else "mismatch"
        evidence["unit"] = f"slot rows are gene PAIRS; '{entry.kind}' reports {'pairs' if ok else 'single genes'}"
    elif slot.unit == "metabolite":
        ok = kind in ("metabolomics", "lipidomics")
        verdicts["unit"] = "fits" if ok else "mismatch"
        evidence["unit"] = f"slot rows are metabolites; '{entry.kind}'"
    else:
        verdicts["unit"] = "unknown"
        evidence["unit"] = f"unit '{slot.unit}' not classified"

    # -- the life stage ---------------------------------------------------------------------------
    stages = stage_words(slot)
    if slot.unit == "host_gene":
        verdicts["stage"] = "not applicable"
        evidence["stage"] = ("the slot asks about the HOST; the parasite stage is the context of "
                             "the experiment, not the thing measured")
    elif not stages:
        verdicts["stage"] = "not applicable"
        evidence["stage"] = "the slot names no life stage"
    else:
        strength = {s: declares_stage(entry, s, columns, pattern) for s in stages}
        in_column = sorted(s for s, w in strength.items() if w == "column")
        in_prose = sorted(s for s, w in strength.items() if w == "prose")
        if in_column:
            verdicts["stage"] = "fits"
            evidence["stage"] = f"the column names {in_column}"
        elif in_prose:
            verdicts["stage"] = "unknown"
            evidence["stage"] = (f"{in_prose} appears only in the prose of the source, never in a "
                                 f"column name -- prose can mention a stage in order to exclude it, "
                                 f"so this needs reading")
        else:
            verdicts["stage"] = "unknown"
            evidence["stage"] = (f"slot asks about {sorted(stages)}; the source declares no stage "
                                 f"in its columns, name, note or coverage -- needs reading")

    # -- access and provenance --------------------------------------------------------------------
    if entry.accession:
        verdicts["access"] = "fits"
        evidence["access"] = f"accession {entry.accession}"
    elif entry.level == "reference" or kind in INFERRED:
        verdicts["access"] = "fits"
        evidence["access"] = "derived or local resource; no external accession is expected"
    else:
        verdicts["access"] = "unknown"
        evidence["access"] = "a measurement with no recorded accession"

    # -- coverage, which is the mapping denominator ------------------------------------------------
    cov = str(entry.coverage or "").strip()
    if cov:
        verdicts["mapping"] = "fits"
        evidence["mapping"] = cov[:120]
    else:
        verdicts["mapping"] = "unknown"
        evidence["mapping"] = "no coverage recorded, so the denominator is unknown"

    return verdicts, evidence


def admission(verdicts: dict) -> str:
    """One word for the row: a mismatch blocks, an unknown queues, otherwise admitted."""
    if any(v == "mismatch" for v in verdicts.values()):
        return "INELIGIBLE"
    if any(v == "unknown" for k, v in verdicts.items() if k != "evidence_grade"):
        return "needs reading"
    return "admitted"


def main() -> int:
    by_key = {e.key: e for e in datasets.REGISTRY}
    by_slot = {s.key: s for s in S.all_slots()}
    cand = pd.read_csv(CANDIDATES, sep="\t")

    rows = []
    for r in cand.itertuples():
        slot = by_slot.get(r.slot)
        entry = by_key.get(r.source_id)
        if slot is None or entry is None:
            rows.append({"slot": r.slot, "dataset_id": r.dataset_id, "source_id": r.source_id,
                         "admission": "unresolved",
                         "note": ("slot not declared" if slot is None else "source not registered")})
            continue
        columns = tuple(getattr(entry, "columns", ()) or ())
        verdicts, evidence = check(slot, entry, columns, str(getattr(r, "pattern", "") or ""))
        row = {"slot": r.slot, "organism": slot.organism, "axis": slot.axis, "unit": slot.unit,
               "slot_context": slot.context, "role": slot.role,
               "dataset_id": r.dataset_id, "source_id": r.source_id,
               "source_organism": entry.organism, "assay": entry.kind,
               "assay_normalised": normalise_kind(entry.kind), "level": entry.level,
               "installed_default": r.installed_default}
        row.update({f"{k}_verdict": v for k, v in verdicts.items()})
        row.update({f"{k}_evidence": evidence[k] for k in verdicts})
        row["admission"] = admission(verdicts)
        rows.append(row)

    frame = pd.DataFrame(rows)
    frame.to_csv(os.path.join(OUT, "biological_eligibility.tsv"), sep="\t", index=False)

    # -- the assay spelling drift, reported rather than silently folded --------------------------
    spell = defaultdict(set)
    for e in datasets.REGISTRY:
        spell[normalise_kind(e.kind)].add(e.kind)
    drift = {k: sorted(v) for k, v in spell.items() if len(v) > 1}
    with open(os.path.join(OUT, "assay_spelling_drift.json"), "w") as fh:
        json.dump(drift, fh, indent=2)

    print(f"candidates examined: {len(frame)}")
    print(frame["admission"].value_counts().to_string())
    print()
    for col in ("organism", "quantity", "assay", "unit", "stage", "access", "mapping"):
        c = Counter(frame.get(f"{col}_verdict", pd.Series(dtype=str)).dropna())
        print(f"  {col:14s} {dict(c)}")
    print()
    print("assay spellings that are one assay:", drift)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
