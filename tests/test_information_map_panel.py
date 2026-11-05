"""The information-map panel shows observed facts, and an unmeasured trait reads as unmeasured.

The failure these tests exist to prevent is not a crash: it is a gene with no recorded host-facing experiment
appearing to have been shown NOT to face the host. So the assertions are about wording and colour as much as
about structure -- "not measured" rather than "no", a measured negative distinguishable from an empty cell, and
the banner stating that no predicted value is present while that is still true of the release.

Headless: conftest pins QT_QPA_PLATFORM=offscreen.
"""
import pandas as pd
import pytest
from PyQt6 import QtWidgets

from parascape import information_map as IM
from parascape.information_map_panel import InformationMapPanel, value_text

GENES = ['TGME49_000001', 'TGME49_000002', 'TGME49_000003']
CONTROLS = (QtWidgets.QAbstractButton, QtWidgets.QLineEdit, QtWidgets.QComboBox, QtWidgets.QSpinBox,
            QtWidgets.QDoubleSpinBox, QtWidgets.QListWidget, QtWidgets.QTableWidget)


def _ours(widget) -> bool:
    """Whether this widget is a control the panel put there, rather than one Qt supplies.

    Excluded: the text field inside a spin or combo box (part of that control, which carries the tooltip),
    a table's corner button and a tab bar's scroll arrows -- Qt builds those, so they are not ours to explain.
    """
    parent = widget.parent()
    if isinstance(parent, (QtWidgets.QAbstractSpinBox, QtWidgets.QComboBox)):
        return False
    if isinstance(widget, QtWidgets.QAbstractButton) and isinstance(parent, (QtWidgets.QAbstractScrollArea,
                                                                            QtWidgets.QTabBar)):
        return False
    return isinstance(widget, CONTROLS)


def _rows(*items):
    return [dict(zip(IM.EVIDENCE_COLUMNS, item)) for item in items]


@pytest.fixture
def release(tmp_path, monkeypatch):
    """A tiny published release on disk: gene 1 has facts of both tiers, gene 3 no tier-1 conclusion at all.

    The module's repository root is pointed at the temporary directory as well, so a table missing from this
    release is reported missing instead of being found in the checkout's own `results/`.
    """
    monkeypatch.setattr(IM, '_REPO', tmp_path)
    evidence = _rows(
        ('TGME49_000001', 'host_facing', 'yes', 0.95, 'observed-direct', 'PMID:1', 'curated host interactions'),
        ('TGME49_000001', 'compartment_hyperlopit', 'dense granules', 0.90, 'observed-curated', 'PMID:2', 'hyperLOPIT'),
        ('TGME49_000002', 'host_facing', 'no', 0.90, 'observed-curated', 'PMID:3', 'gold negatives'))
    pairs = pd.DataFrame([('TGME49_000001', 'TGME49_000002', 'parasite', 0.97, 'observed-direct', 'XL-MS 1'),
                          ('TGME49_000001', 'STAT6', 'host', 0.90, 'observed-curated', 'PMID:4')],
                         columns=IM.PAIR_COLUMNS)
    values, probability, source = IM.build(evidence, GENES, pairs=pairs)
    values.to_parquet(tmp_path / 'values.parquet')
    probability.to_parquet(tmp_path / 'probability.parquet')
    source.to_parquet(tmp_path / 'source.parquet')
    IM.evidence_frame(evidence).to_parquet(tmp_path / 'evidence.parquet')
    pairs.to_parquet(tmp_path / 'pairs.parquet')
    # Tier-2 measurement tables, written into the release directory so the repository's own are not used.
    pd.DataFrame({'n_tm': [0, 2, 1], 'has_signal_peptide': [True, False, False]},
                 index=pd.Index(GENES, name='gene_id')).to_parquet(tmp_path / 'sequence_traits.parquet')
    pd.DataFrame({'dom_IPR000008': [True, False, False], 'n_interpro_domains': [1, 0, 0]},
                 index=pd.Index(GENES, name='gene_id')).to_parquet(tmp_path / 'domain_traits.parquet')
    return tmp_path


@pytest.fixture
def panel(release):
    QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    widget = InformationMapPanel(published=IM.load_published(release),
                                 products={'TGME49_000001': 'dense granule protein GRA-test'})
    widget.resize(900, 900)
    yield widget
    widget.close()


# --------------------------------------------------------------------------- the loader
def test_the_loader_returns_what_exists_and_names_what_is_absent(release):
    (release / 'pairs.parquet').unlink()
    (release / 'domain_traits.parquet').unlink()
    published = IM.load_published(release)
    assert not published.values.empty and published.pairs.empty
    assert any('pairs.parquet' in m for m in published.missing)
    assert any('domain' in m for m in published.missing)
    # And the rest is still usable: the panel opens on a partial release.
    assert published.gene_facts('TGME49_000001').shape[0] >= 2


def test_an_absent_release_is_reported_rather_than_raised(tmp_path):
    published = IM.load_published(tmp_path / 'nothing here')
    assert published.values.empty and published.missing and 'Looked at' in published.missing[0]
    assert published.summary() == 'No published information map was found.'


def test_the_catalogue_separates_conclusions_from_measurements(release):
    published = IM.load_published(release)
    assert 'host_facing' in published.columns(1) and 'compartment_hyperlopit' in published.columns(1)
    assert {'n_tm', 'dom_IPR000008'} <= set(published.columns(2))


def test_no_cell_of_the_published_release_is_a_prediction(release):
    published = IM.load_published(release)
    assert published.predicted_cells().empty
    assert 'observed fact' in published.banner() and '90' in published.banner()


def test_a_filled_cell_would_turn_the_banner_round(release):
    """The banner is measured, not asserted: a release with a filled cell must say so."""
    source = pd.read_parquet(release / 'source.parquet')
    source.loc['TGME49_000003', 'host_facing'] = 'filled:gbt@v1'
    source.to_parquet(release / 'source.parquet')
    published = IM.load_published(release)
    assert len(published.predicted_cells()) == 1
    assert 'predicted' in published.banner()


def test_the_shipped_release_loads_if_it_is_present():
    """The real 2026-11-04 release, when the repository holds it: no predicted cell, and genes to show."""
    published = IM.load_published()
    if published.values.empty:
        pytest.skip('no published information map in this checkout')
    assert published.predicted_cells().empty
    assert len(published.genes) > 1000 and 'host_facing' in published.columns(1)


# --------------------------------------------------------------------------- the gene view
def test_the_panel_builds_and_shows_the_observed_facts_with_their_tiers(panel):
    assert 'observed fact' in panel.banner.text()
    panel.search.setText('TGME49_000001')
    panel.find_genes()
    assert panel.gene == 'TGME49_000001'
    facts = panel.gene_facts().set_index('column')
    assert facts.loc['host_facing', 'tier'] == 1 and facts.loc['host_facing', 'probability'] == 0.95
    assert facts.loc['host_facing', 'source_tier'] == 'observed-direct'
    assert facts.loc['compartment_hyperlopit', 'source_tier'] == 'observed-curated'
    assert facts.loc['n_tm', 'tier'] == 2                       # measurements are grouped as tier 2
    assert list(facts.tier) == sorted(facts.tier)               # conclusions first, then measurements
    shown = {panel.facts.item(i, 0).text(): i for i in range(panel.facts.rowCount())}
    row = shown['host_facing']
    assert 'yes' in panel.facts.item(row, 2).text()
    assert 'PMID:1' in panel.facts.item(row, 5).text()          # the evidence behind the cell
    assert 'conclusion' in panel.facts.item(row, 1).text()


def test_a_gene_can_be_found_by_its_product(panel):
    panel.search.setText('GRA-test')
    panel.find_genes()
    assert panel.gene == 'TGME49_000001' and panel.matches.count() == 1


def test_an_unmeasured_trait_reads_as_not_measured_and_never_as_a_negative(panel):
    panel.show_gene('TGME49_000003')                            # nothing tier 1 was ever measured on it
    assert 'host_facing' not in set(panel.gene_facts().column)  # an empty cell is not a row by default
    panel.absent.setChecked(True)
    facts = panel.gene_facts().set_index('column')
    assert not facts.loc['host_facing', 'observed']
    assert facts.loc['host_facing', 'value'] == IM.NOT_MEASURED
    assert pd.isna(facts.loc['host_facing', 'probability'])
    texts = {panel.facts.item(i, 0).text(): panel.facts.item(i, 2) for i in range(panel.facts.rowCount())}
    empty = texts['host_facing']
    assert empty.text() == 'not measured' and 'no' != empty.text()
    assert 'not a' in empty.toolTip() and 'negative' in empty.toolTip()
    # A measured negative is a different thing and must not look like the empty cell.
    negative = texts['dom_IPR000008']
    assert negative.text().startswith('no —') and 'measured negative' in negative.text()
    assert negative.foreground().color() != empty.foreground().color()
    assert empty.font().italic() and not negative.font().italic()
    assert 'not measured' in panel.gene_note.text()


def test_a_measured_negative_of_a_curated_source_keeps_its_probability(panel):
    panel.show_gene('TGME49_000002')
    facts = panel.gene_facts().set_index('column')
    assert facts.loc['host_facing', 'value'] is False or str(facts.loc['host_facing', 'value']) == 'False'
    assert facts.loc['host_facing', 'probability'] == 0.90
    assert facts.loc['host_facing', 'source_tier'] == 'observed-curated'


def test_value_text_separates_the_three_states():
    assert value_text({'observed': False, 'value': IM.NOT_MEASURED}) == 'not measured'
    assert value_text({'observed': True, 'value': False}).startswith('no —')
    assert value_text({'observed': True, 'value': True}).startswith('yes')
    assert value_text({'observed': True, 'value': 'dense granules'}) == 'dense granules'
    assert value_text({'observed': True, 'value': 0.974591747815646}) == '0.9746'
    assert len(value_text({'observed': True, 'value': 'M' * 500})) < 130       # a whole protein is truncated


def test_the_partners_of_a_gene_are_shown_with_their_probabilities(panel):
    panel.show_gene('TGME49_000001')
    partners = {panel.partners.item(i, 0).text(): i for i in range(panel.partners.rowCount())}
    assert set(partners) == {'TGME49_000002', 'STAT6'}
    assert panel.partners.item(partners['TGME49_000002'], 2).text().startswith('0.97')
    assert panel.partners.item(partners['STAT6'], 1).text() == 'host'
    panel.show_gene('TGME49_000003')
    assert panel.partners.rowCount() == 0 and 'not the same as having none' in panel.partner_note.text()


def test_a_gene_absent_from_the_map_says_so_rather_than_showing_blanks(panel):
    panel.search.setText('TGME49_999999')
    panel.find_genes()
    assert panel.matches.count() == 0 and 'No gene of the map matches' in panel.match_note.text()
    assert panel.facts.rowCount() == 0
    panel.show_gene('TGME49_999999')
    assert panel.facts.rowCount() == 0 and 'not a row of this map' in panel.gene_note.text()
    assert 'not a statement about the gene' in panel.gene_note.text()


# --------------------------------------------------------------------------- the trait view
def test_the_trait_view_shows_coverage_by_tier_and_the_value_distribution(panel):
    assert panel.show_trait('host_facing')
    coverage = {panel.coverage.item(i, 0).text(): panel.coverage.item(i, 1).text()
                for i in range(panel.coverage.rowCount())}
    assert coverage['observed-direct'] == '1' and coverage['observed-curated'] == '1'
    assert coverage['any tier'] == '2'
    distribution = {panel.distribution.item(i, 0).text(): panel.distribution.item(i, 1).text()
                    for i in range(panel.distribution.rowCount())}
    assert distribution == {'True': '1', 'False': '1'}
    assert 'not measured' in panel.trait_note.text() and 'tier 1' in panel.trait_note.text()
    assert panel.show_trait('n_tm') and 'tier 2' in panel.trait_note.text()
    assert not panel.show_trait('no_such_column')


# --------------------------------------------------------------------------- the house rule
def test_every_interactive_control_explains_itself(panel):
    panel.search.setText('TGME49_000001')
    panel.find_genes()
    panel.absent.setChecked(True)
    missing = [f'{type(w).__name__}: {getattr(w, "text", lambda: "")()!r}'
               for w in panel.findChildren(QtWidgets.QWidget) if _ours(w) and not w.toolTip().strip()]
    assert not missing, f'controls without a tooltip: {missing}'
    assert len([w for w in panel.findChildren(QtWidgets.QWidget) if _ours(w)]) >= 9
