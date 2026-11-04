"""The information map (instruction 72): gene x trait values with an aligned probability and source table.

Everything is built from one long **evidence** table -- one row per fact, from any source -- and a registry of
**traits** that says how the rows of one trait become one matrix cell. Adding a column is therefore two steps:

1. register it once::

       information_map.register(Trait('host_facing', 'flag', 'Protein is found beyond the parasite, facing the host'))

2. give evidence rows with ``trait='host_facing'`` from any source (a loader, a notebook cell, a new dataset)::

       evidence.append(dict(gene_id='TGME49_270240', trait='host_facing', value='yes', probability=0.95,
                            source_tier='observed-direct', reference='PMID:...', dataset='...'))

``build`` then writes the column into ``values``, ``probability`` and ``source`` with the same empty cells in all
three; other columns are untouched. Predicted (filled) cells use the same evidence format with a
``filled:<strategy>@<version>`` tier, so observed and inferred facts never need a different code path.

Rules kept here: missing stays empty (never zero, never "no"); unlabeled is not negative; when several sources
speak about one cell the highest-probability one is shown and all remain in the evidence table.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

#: Provisional probabilities of observed facts by tier (v0); v1 replaces them by measured replication.
TIER_P = {'observed-direct': 0.95, 'observed-curated': 0.90, 'observed-screen': 0.80, 'observed-proximity': 0.70}
#: Order used to break probability ties: stronger evidence first, then filled, then exploratory.
TIER_RANK = {'observed-direct': 0, 'observed-curated': 1, 'observed-screen': 2, 'observed-proximity': 3}
#: Einar 2026-11-04: a predicted entry is filled only at >= 90 % held-out precision on direct-evidence genes.
FILL_PRECISION = 0.90
EVIDENCE_COLUMNS = ('gene_id', 'trait', 'value', 'probability', 'source_tier', 'reference', 'dataset')
PAIR_COLUMNS = ('gene_a', 'gene_b', 'partner_side', 'probability', 'source_tier', 'reference')
#: Einar 2026-11-04: tier 1 = conclusions that depend on many factors (host-facing, partners, function,
#: compartment) -- what is predicted; tier 2 = measurements (relative expression, protein abundance, fitness,
#: sequence features) -- what predictions are learned from.
TIERS = {1: 'conclusion (predicted from evidence)', 2: 'measurement (evidence)'}
#: Einar 2026-11-04: "relative values that are internally normalized (like housekeeping genes) or normalized to a
#: control are always better". Best first; learners prefer the earliest form available for a quantity.
NORMALIZATIONS = ('relative_to_control', 'internal_reference', 'library_size', 'raw', 'not_applicable', 'unknown')
#: measurement = observed value; prediction = another model's output (never used as an observed fact or as a
#: training label); identity = names and text (never a feature).
ROLES = ('measurement', 'prediction', 'identity')
KINDS = {
    'flag': 'True where a source shows it, False only where a source shows the opposite, empty when untested',
    'category': 'one label per gene (e.g. a compartment); the highest-probability label is shown',
    'terms': 'several terms per gene (e.g. GO terms); the best-supported term is shown and <name>_n counts terms',
    'measure': 'a measured result that may be negative (e.g. a screen phenotype True/False, a number)',
}


@dataclass(frozen=True)
class Trait:
    """One column of the map: its name, how its evidence rows collapse to a cell (KINDS), and what it means."""

    name: str
    kind: str
    description: str
    tier: int = 1
    family: str = ''                    # source family, held out together (a dataset key, 'computed', ...)
    normalization: str = 'not_applicable'
    role: str = 'measurement'
    depends_on: tuple = ()              # tier-2 families a tier-1 trait is derived from (excluded when predicting it)

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError(f'{self.name}: kind must be one of {sorted(KINDS)}')
        if self.tier not in TIERS or self.normalization not in NORMALIZATIONS or self.role not in ROLES:
            raise ValueError(f'{self.name}: tier, normalization or role outside {TIERS}, {NORMALIZATIONS}, {ROLES}')


TRAITS: dict[str, Trait] = {}


def register(trait: Trait) -> Trait:
    """Add a column; registering the same name again with a different definition is an error."""
    if trait.name in TRAITS and TRAITS[trait.name] != trait:
        raise ValueError(f'{trait.name} is already registered differently')
    TRAITS[trait.name] = trait
    return trait


for _trait in (
        Trait('host_facing', 'flag', 'Protein (or part of it) is found beyond the parasite plasma membrane, facing the host',
              family='host_facing_truth', depends_on=('pvm_proximity',)),
        Trait('compartment_hyperlopit', 'category', 'hyperLOPIT marker niche (Barylyuk 2020, curated)',
              family='lopit_tgon', depends_on=('lopit_tgon',)),
        Trait('compartment_go_cc', 'terms', 'Experimental GO cellular component', family='quickgo'),
        Trait('function_go_bp', 'terms', 'Experimental GO biological process', family='quickgo'),
        Trait('function_go_mf', 'terms', 'Experimental GO molecular function', family='quickgo')):
    register(_trait)


def evidence_frame(rows) -> pd.DataFrame:
    """Validate evidence rows (dicts, tuples in EVIDENCE_COLUMNS order, or a DataFrame)."""
    rows = rows if isinstance(rows, pd.DataFrame) else list(rows)
    if not len(rows):
        return pd.DataFrame({c: pd.Series(dtype=float if c == 'probability' else object) for c in EVIDENCE_COLUMNS})
    frame = rows.copy() if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows, columns=None)
    if not isinstance(rows, pd.DataFrame) and len(frame.columns) == len(EVIDENCE_COLUMNS) and \
            list(frame.columns) == list(range(len(EVIDENCE_COLUMNS))):
        frame.columns = list(EVIDENCE_COLUMNS)
    missing = set(EVIDENCE_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(f'evidence rows lack {sorted(missing)}')
    frame = frame[list(EVIDENCE_COLUMNS)].copy()
    frame['probability'] = frame['probability'].astype(float)
    if ((frame.probability < 0) | (frame.probability > 1)).any():
        raise ValueError('probabilities must lie in [0, 1]')
    unknown = set(frame.trait) - set(TRAITS)
    if unknown:
        raise ValueError(f'evidence for unregistered traits {sorted(unknown)}; register() them first')
    frame['value'] = frame['value'].astype(str)
    return frame


def _flag(text: str):
    """'yes'/'no' and their synonyms; anything else leaves the cell empty, because untested is not negative."""
    lowered = str(text).strip().lower()
    if lowered in ('yes', 'true', '1', 'host_facing'):
        return True
    if lowered in ('no', 'false', '0', 'not_exported'):
        return False
    return pd.NA


def conflicts(evidence) -> pd.DataFrame:
    """Genes whose flag traits have both a positive and a negative source; shown, never averaged away."""
    ev = evidence_frame(evidence)
    flags = [name for name, trait in TRAITS.items() if trait.kind == 'flag']
    ev = ev[ev.trait.isin(flags)].assign(flag=lambda d: d.value.map(_flag))
    both = ev.groupby(['gene_id', 'trait']).flag.agg(lambda s: {True, False} <= set(s.dropna()))
    return ev.merge(both[both].rename('conflict').reset_index(), on=['gene_id', 'trait'])


def _strongest(rows: pd.DataFrame) -> pd.DataFrame:
    rank = rows.source_tier.map(TIER_RANK).fillna(len(TIER_RANK))
    return (rows.assign(_rank=rank).sort_values(['gene_id', 'probability', '_rank', 'value'],
                                                ascending=[True, False, True, True])
            .drop_duplicates('gene_id').set_index('gene_id'))


def _column(trait: Trait, rows: pd.DataFrame, genes: pd.Index) -> dict:
    """{column name: (values, probability, source)} for one trait."""
    if rows.empty:
        empty = pd.Series(np.nan, index=genes)
        return {trait.name: (empty.astype(object), empty, empty.astype(object))}
    if trait.kind == 'terms':
        support = (rows.groupby(['gene_id', 'value'])
                   .agg(n=('reference', 'nunique'), probability=('probability', 'max'),
                        source_tier=('source_tier', 'first')).reset_index()
                   .sort_values(['gene_id', 'n', 'value'], ascending=[True, False, True]))
        top = support.drop_duplicates('gene_id').set_index('gene_id')
        count = support.groupby('gene_id').size()
        return {trait.name: (top.value.reindex(genes), top.probability.reindex(genes), top.source_tier.reindex(genes)),
                trait.name + '_n': (count.reindex(genes).astype('Int64'), top.probability.reindex(genes),
                                    top.source_tier.reindex(genes))}
    best = _strongest(rows)
    value = best.value
    if trait.kind == 'flag':
        value = value.map(_flag).astype('boolean')
    elif trait.kind == 'measure' and set(value.unique()) <= {'True', 'False'}:
        value = value.map({'True': True, 'False': False}).astype('boolean')
    return {trait.name: (value.reindex(genes), best.probability.reindex(genes), best.source_tier.reindex(genes))}


def partner_columns(pairs: pd.DataFrame, genes: pd.Index) -> dict:
    """Per-gene summaries of the pairs table: best partner, its probability, partners >= 0.95 and >= 0.5."""
    if pairs is None or pairs.empty:
        return {}
    pa = pairs.assign(partner=np.where(pairs.partner_side == 'host', 'host:' + pairs.gene_b.astype(str), pairs.gene_b))
    fwd = pa[['gene_a', 'partner', 'probability', 'source_tier']].rename(columns={'gene_a': 'gene'})
    rev = pa[pa.partner_side == 'parasite'][['gene_b', 'gene_a', 'probability', 'source_tier']].rename(
        columns={'gene_b': 'gene', 'gene_a': 'partner'})
    long = (pd.concat([fwd, rev], ignore_index=True)
            .sort_values(['gene', 'partner', 'probability'], ascending=[True, True, False])
            .drop_duplicates(['gene', 'partner']))
    best = long.sort_values(['gene', 'probability', 'partner'], ascending=[True, False, True]).drop_duplicates('gene').set_index('gene')
    has = genes.isin(best.index)
    p, s = best.probability.reindex(genes), best.source_tier.reindex(genes)
    out = {'partner_best': (best.partner.reindex(genes), p, s), 'partner_best_prob': (p, p, s)}
    for name, threshold in (('partner_n_ge95', 0.95), ('partner_n_ge50', 0.5)):
        n = long[long.probability >= threshold].groupby('gene').size().reindex(genes).fillna(0).astype('Int64')
        out[name] = (n.where(has), p, s)
    return out


def build(evidence, genes, pairs: pd.DataFrame | None = None, order=None):
    """Assemble ``values``, ``probability`` and ``source`` (same shape, same empty cells) from evidence rows.

    Columns follow registration order (or ``order``), each trait's helper columns next to it, partner summaries
    after the trait named ``partners`` in ``order`` or at the end.
    """
    genes = pd.Index(genes, name='gene_id')
    ev = evidence_frame(evidence)
    ev = ev[ev.gene_id.isin(genes)]
    names = list(order) if order is not None else [t for t in TRAITS if t in set(ev.trait)]
    columns = {}
    for name in names:
        if name == 'partners':
            columns.update(partner_columns(pairs, genes)); continue
        columns.update(_column(TRAITS[name], ev[ev.trait == name], genes))
    if 'partners' not in names:
        columns.update(partner_columns(pairs, genes))
    values = pd.DataFrame({c: v for c, (v, _, _) in columns.items()}, index=genes)
    probability = pd.DataFrame({c: p for c, (_, p, _) in columns.items()}, index=genes).astype(float)
    source = pd.DataFrame({c: s for c, (_, _, s) in columns.items()}, index=genes)
    empty = values.isna()
    probability, source = probability.mask(empty), source.mask(empty)
    assert (probability.isna() == empty).all().all() and (source.isna() == empty).all().all()
    return values, probability, source
