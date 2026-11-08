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

The second half of the module reads the **published** tables back (``load_published``), so the application can
show what is known about a gene without rebuilding anything. It resolves the release directory, loads the three
aligned tables plus the long evidence and pairs tables and the tier-2 measurement tables, and reports in
``PublishedMap.missing`` whatever it could not find rather than failing: a release without, say, the domain
table is still worth browsing. Nothing is inferred on the way in -- every cell it hands back is an observed
fact, because as of the 2026-11-04 release no strategy has reached the 90 % fill bar (``FILL_PRECISION``), and
``PublishedMap.predicted_cells`` is the check that this is still true rather than an assumption.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

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
              family='host_facing_truth',
              # Declared 2026-11-04 after the first scoring run: these assays are so close to the conclusion
              # that using them as evidence is near-circular, so they are held out when host_facing is predicted.
              depends_on=('pvm_proximity', 'hosttx_effectors', 'organelle_surface_turboid', 'pubmed')),
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


# --------------------------------------------------------------------------- reading a published release
#: Override for the published information-map directory, honoured absolutely (the `paths` convention).
PUBLISHED_ENV = 'PARASCAPE_INFORMATION_MAP'
#: Files of a published release. Only ``values`` is required; the rest are reported when absent.
MATRIX_FILES = ('values', 'probability', 'source', 'evidence', 'pairs')
#: Tier-2 measurement tables: name -> (directory glob, file name). One row per gene, no probability column,
#: so their cells carry the dataset as their source and no probability: these are measurements, not claims.
TIER2_TABLES = {
    'sequence': ('sequence_traits_*', 'sequence_traits.parquet'),
    'domain': ('domain_traits_*', 'domain_traits.parquet'),
}
#: A source tier that starts with one of these is a prediction, not an observed fact.
PREDICTED_PREFIXES = ('filled:', 'exploratory:')
#: What an empty cell says. The distinction this whole map exists for: nothing is known, which is not "no".
NOT_MEASURED = 'not measured'
_REPO = Path(__file__).resolve().parent.parent


def _newest_directory(parent: Path, pattern: str) -> Path | None:
    """The last directory matching ``pattern`` by name, so a dated release sorts to the newest."""
    if not parent.is_dir():
        return None
    hits = sorted((p for p in parent.glob(pattern) if p.is_dir()), key=lambda p: p.name)
    return hits[-1] if hits else None


def published_candidates(root=None) -> list:
    """Every directory a published information map could live in, in priority order.

    Listed rather than resolved silently, because a resolver that picks the wrong release is harder to
    debug than one that says where it looked -- the rule `paths` was written for.
    """
    if root is not None:
        return [Path(root)]
    out = []
    override = os.environ.get(PUBLISHED_ENV)
    if override:
        out.append(Path(override))
    try:
        from . import paths
        out.append(Path(paths.data_dir()) / 'information_map')
    except Exception:                       # noqa: BLE001 -- the loader must work without the cache resolver
        pass
    newest = _newest_directory(_REPO / 'results', 'information_map_v0_*')
    if newest is not None:
        out.append(newest)
    return out


def published_root(root=None) -> Path | None:
    """The first candidate directory that actually holds ``values.parquet``; None when there is none."""
    for candidate in published_candidates(root):
        if (candidate / 'values.parquet').is_file():
            return candidate
    return None


def tier2_path(name: str, root=None) -> Path | None:
    """Where one tier-2 measurement table lives: beside the published release, or in ``results/``."""
    glob, filename = TIER2_TABLES[name]
    # The published release first, then the repository's results. An installed copy has no `results`
    # directory, so without the release directory here the tier-2 rows would be missing for every user
    # who did not clone the repository.
    parents = [Path(root) if root is not None else published_root(), _REPO / 'results']
    for parent in [p for p in parents if p is not None]:
        direct = parent / filename
        if direct.is_file():
            return direct
        found = _newest_directory(parent, glob)
        if found is not None and (found / filename).is_file():
            return found / filename
    return None


def _domain_descriptions(path: Path) -> dict:
    lookup = path.parent / 'domain_lookup.tsv'
    if not lookup.is_file():
        return {}
    table = pd.read_csv(lookup, sep='\t', dtype=str).fillna('')
    return dict(zip(table['column'], (table['description'] + ' (' + table['accession'] + ')').str.strip()))


def column_catalogue(frames: dict, descriptions=None) -> pd.DataFrame:
    """One row per column of the loaded tables: its tier, kind, table and description.

    The tier decides how the application groups a column -- tier 1 conclusions first, then the tier-2
    measurements they would be learned from (instruction 72). It is taken from the trait registry where the
    column is registered, and otherwise from the table it came from, so the grouping never depends on whether
    a build notebook happened to run first in this process.
    """
    descriptions = descriptions or {}
    rows = []
    for table, frame in frames.items():
        for column in frame.columns:
            base = column[:-2] if column.endswith('_n') and column[:-2] in TRAITS else column
            trait = TRAITS.get(base)
            if table == 'information map':
                # Screen results are measurements (instruction 72 tier 2) whatever a build script registered.
                tier = 2 if column.startswith('screen_') else (trait.tier if trait is not None else 1)
            else:
                tier = 2
            description = descriptions.get(column) or (trait.description if trait is not None else '')
            if column.endswith('_n') and trait is not None:
                description = f'Number of distinct terms recorded for {base}'
            rows.append(dict(column=column, tier=tier, kind=trait.kind if trait is not None else '',
                             family=trait.family if trait is not None else table, table=table,
                             description=description))
    return pd.DataFrame(rows, columns=['column', 'tier', 'kind', 'family', 'table', 'description'])


@dataclass
class PublishedMap:
    """A published information map, read back for display: observed facts, their probability and their source.

    ``values``, ``probability`` and ``source`` are the three aligned tables; ``evidence`` and ``pairs`` are the
    long-form record behind them; ``tier2`` holds the measurement tables keyed by name. ``missing`` names every
    file that was not found, so a partial release is usable and says what it lacks.
    """

    root: Path | None = None
    values: pd.DataFrame = field(default_factory=pd.DataFrame)
    probability: pd.DataFrame = field(default_factory=pd.DataFrame)
    source: pd.DataFrame = field(default_factory=pd.DataFrame)
    evidence: pd.DataFrame = field(default_factory=pd.DataFrame)
    pairs: pd.DataFrame = field(default_factory=pd.DataFrame)
    tier2: dict = field(default_factory=dict)
    catalogue: pd.DataFrame = field(default_factory=pd.DataFrame)
    missing: tuple = ()

    @property
    def genes(self) -> pd.Index:
        """Every gene of the map (the rows of ``values``), tier-2 only genes included."""
        index = self.values.index
        for frame in self.tier2.values():
            index = index.union(frame.index)
        return pd.Index(index, name='gene_id')

    def has_gene(self, gene_id: str) -> bool:
        """Whether the gene is a row of this map at all -- absent from the map is not absent from biology."""
        return str(gene_id) in set(self.genes)

    def columns(self, tier=None) -> list:
        """Column names, optionally of one tier only, in the order the published tables carry them."""
        if self.catalogue.empty:
            return []
        rows = self.catalogue if tier is None else self.catalogue[self.catalogue.tier == tier]
        return rows.column.tolist()

    def gene_facts(self, gene_id: str, include_absent: bool = False) -> pd.DataFrame:
        """Every trait of one gene: value, probability, source tier, table and tier.

        Cells with no observed fact are left out unless ``include_absent``, and when they are shown their value
        is ``NOT_MEASURED`` with no probability: an empty cell is the absence of a measurement, never a "no".
        A measured negative (``False`` in a screen or a curated non-export) is a fact and is always shown.
        """
        gene_id = str(gene_id)
        rows = []
        for record in self.catalogue.itertuples(index=False):
            frame = self.values if record.table == 'information map' else self.tier2.get(record.table)
            if frame is None or gene_id not in frame.index:
                continue
            value = frame.at[gene_id, record.column]
            observed = not _is_blank(value)
            if not observed and not include_absent:
                continue
            probability = source_tier = None
            if record.table == 'information map' and observed:
                if gene_id in self.probability.index and record.column in self.probability.columns:
                    probability = self.probability.at[gene_id, record.column]
                if gene_id in self.source.index and record.column in self.source.columns:
                    source_tier = self.source.at[gene_id, record.column]
            rows.append(dict(column=record.column, tier=record.tier, table=record.table,
                             description=record.description, observed=observed,
                             value=value if observed else NOT_MEASURED,
                             probability=float(probability) if probability is not None and not pd.isna(probability)
                             else np.nan,
                             source_tier=str(source_tier) if source_tier is not None and not pd.isna(source_tier)
                             else ('measured value (no probability recorded)' if observed and record.tier == 2
                                   else '')))
        return pd.DataFrame(rows, columns=['column', 'tier', 'table', 'description', 'observed', 'value',
                                           'probability', 'source_tier'])

    def gene_evidence(self, gene_id: str, trait=None) -> pd.DataFrame:
        """The long-form evidence rows behind one gene: every source, including the ones that lost the cell."""
        if self.evidence.empty:
            return self.evidence
        rows = self.evidence[self.evidence.gene_id.astype(str) == str(gene_id)]
        if trait is not None:
            base = trait[:-2] if trait.endswith('_n') else trait
            rows = rows[rows.trait.isin({trait, base})]
        return rows.sort_values(['trait', 'probability'], ascending=[True, False])

    def gene_partners(self, gene_id: str, limit: int | None = None) -> pd.DataFrame:
        """One gene's recorded partners from the pairs table, both directions, strongest first."""
        if self.pairs.empty:
            return pd.DataFrame(columns=['partner', 'side', 'probability', 'source_tier', 'reference'])
        gene_id = str(gene_id)
        forward = self.pairs[self.pairs.gene_a.astype(str) == gene_id].assign(
            partner=lambda d: d.gene_b.astype(str), side=lambda d: d.partner_side)
        backward = self.pairs[(self.pairs.gene_b.astype(str) == gene_id) & (self.pairs.partner_side == 'parasite')].assign(
            partner=lambda d: d.gene_a.astype(str), side='parasite')
        both = (pd.concat([forward, backward], ignore_index=True)[
                    ['partner', 'side', 'probability', 'source_tier', 'reference']]
                .sort_values(['probability', 'partner'], ascending=[False, True])
                .drop_duplicates(['partner', 'reference']))
        return both.head(limit) if limit else both

    def coverage(self, column: str) -> pd.DataFrame:
        """How many genes carry an observed fact for one column, split by the source tier of the cell.

        The tier counted is the one that won the cell (the strongest source); the sources that lost it stay in
        ``evidence``, so these counts are of cells, not of evidence rows.

        Tier-2 tables have no source column, so their coverage is reported as one row naming the table: those
        cells are measurements of one dataset, and inventing a tier for them would be a claim of its own.
        """
        table = self._table_of(column)
        if table is None:
            return pd.DataFrame(columns=['source_tier', 'genes'])
        if table != 'information map':
            observed = ~self.tier2[table][column].map(_is_blank)
            return pd.DataFrame([{'source_tier': f'{table} measurement', 'genes': int(observed.sum())}])
        observed = ~self.values[column].map(_is_blank)
        tiers = self.source[column][observed] if column in self.source.columns else pd.Series(dtype=object)
        out = (tiers.fillna('source not recorded').value_counts().rename_axis('source_tier')
               .rename('genes').reset_index())
        return pd.concat([out, pd.DataFrame([{'source_tier': 'any tier', 'genes': int(observed.sum())}])],
                         ignore_index=True)

    def distribution(self, column: str, limit: int = 25) -> pd.DataFrame:
        """The values one column takes and how many genes hold each, commonest first."""
        series = self.column(column)
        if series is None:
            return pd.DataFrame(columns=['value', 'genes'])
        observed = series[~series.map(_is_blank)]
        counts = observed.astype(str).value_counts()
        out = counts.rename_axis('value').rename('genes').reset_index().head(limit)
        return out

    def column(self, column: str):
        """One column of the map as a Series over genes, from whichever table holds it."""
        table = self._table_of(column)
        if table is None:
            return None
        return (self.values if table == 'information map' else self.tier2[table])[column]

    def _table_of(self, column: str):
        if self.catalogue.empty:
            return None
        hit = self.catalogue[self.catalogue.column == column]
        return None if hit.empty else hit.table.iloc[0]

    def predicted_cells(self) -> pd.DataFrame:
        """Cells whose source is a prediction (``filled:`` or ``exploratory:``) rather than an observation.

        Empty for every release so far, and the application says so in a banner. It is computed rather than
        assumed, so the banner stops being true the moment a filled release is loaded.
        """
        if self.source.empty:
            return pd.DataFrame(columns=['gene_id', 'column', 'source_tier'])
        text = self.source.astype(object).where(self.source.notna(), '')
        hit = text.apply(lambda col: col.astype(str).str.startswith(PREDICTED_PREFIXES))
        rows = [dict(gene_id=str(self.source.index[r]), column=str(self.source.columns[c]),
                     source_tier=str(text.iat[r, c]))
                for r, c in zip(*np.where(hit.to_numpy()))]
        return pd.DataFrame(rows, columns=['gene_id', 'column', 'source_tier'])

    def banner(self) -> str:
        """The sentence the panel shows above everything: what kind of cells these are, and what is absent."""
        predicted = len(self.predicted_cells())
        if predicted:
            return (f'{predicted:,} cells of this release are predicted (filled or exploratory). '
                    'Check their source tier before reading one as a fact.')
        return ('Every cell shown is an observed fact. No predicted values are present: no strategy has yet '
                f'reached the {FILL_PRECISION:.0%} held-out precision bar that lets a predicted entry be '
                'written into the map, so no gap has been filled.')

    def summary(self) -> str:
        """One line of what was loaded, for the panel's status label."""
        if self.values.empty and not self.tier2:
            return 'No published information map was found.'
        parts = [f'{len(self.genes):,} genes', f'{len(self.catalogue):,} columns '
                 f'({len(self.columns(1))} tier 1, {len(self.columns(2))} tier 2)',
                 f'{len(self.evidence):,} evidence rows', f'{len(self.pairs):,} recorded pairs']
        return ' · '.join(parts) + (f' · from {self.root}' if self.root else '')


def _is_blank(value) -> bool:
    """Whether a cell holds no fact: NA, or an empty string. ``False`` and ``0`` are facts, not blanks."""
    if value is None or value is pd.NA:
        return True
    try:
        if pd.isna(value):
            return True
    except (TypeError, ValueError):
        return False
    return isinstance(value, str) and not value.strip()


def load_published(root=None, tier2=True) -> PublishedMap:
    """Read a published information map back from disk, reporting rather than raising what is absent.

    ``root`` names a release directory; without it the override, the installed cache and the newest
    ``results/information_map_v0_*`` are tried in that order (``published_candidates``). Tier-2 measurement
    tables are looked for beside the release and in ``results/``. A file that is missing or unreadable is
    recorded in ``PublishedMap.missing`` and the rest is still returned, because a partial release is worth
    browsing and a panel that refuses to open explains nothing.
    """
    found = published_root(root)
    missing, frames = [], {}
    if found is None:
        looked = ', '.join(str(p) for p in published_candidates(root)) or '(no candidate directories)'
        return PublishedMap(missing=(f'No published information map (values.parquet) found. Looked at: {looked}.',))
    loaded = {}
    for name in MATRIX_FILES:
        path = found / f'{name}.parquet'
        if not path.is_file():
            missing.append(f'{path} is absent, so its part of the map is not shown.')
            continue
        try:
            loaded[name] = pd.read_parquet(path)
        except (OSError, ValueError) as error:
            missing.append(f'{path} could not be read ({error}).')
    values = loaded.get('values', pd.DataFrame())
    if not values.empty:
        values.index = values.index.astype(str)
        frames['information map'] = values
    tables = {}
    if tier2:
        for name in TIER2_TABLES:
            path = tier2_path(name, root)
            if path is None:
                missing.append(f'The tier-2 {name} table ({TIER2_TABLES[name][1]}) was not found, '
                               'so its measurements are not shown.')
                continue
            try:
                frame = pd.read_parquet(path)
            except (OSError, ValueError) as error:
                missing.append(f'{path} could not be read ({error}).')
                continue
            frame.index = frame.index.astype(str)
            tables[name] = frame
            frames[name] = frame
            if name == 'domain':
                loaded.setdefault('_descriptions', {}).update(_domain_descriptions(path))
    for name in ('probability', 'source'):
        if name in loaded:
            loaded[name].index = loaded[name].index.astype(str)
    # Filled columns (instruction 72) carry their question in words beside the release.
    described = found / 'column_descriptions.json'
    if described.is_file():
        import json
        loaded.setdefault('_descriptions', {}).update(json.loads(described.read_text()))
    return PublishedMap(root=found, values=values,
                        probability=loaded.get('probability', pd.DataFrame()),
                        source=loaded.get('source', pd.DataFrame()),
                        evidence=loaded.get('evidence', pd.DataFrame()),
                        pairs=loaded.get('pairs', pd.DataFrame()),
                        tier2=tables,
                        catalogue=column_catalogue(frames, loaded.get('_descriptions', {})),
                        missing=tuple(missing))
