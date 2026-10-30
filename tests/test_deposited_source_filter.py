"""Source filtering preserves provenance and cannot export measurements from an earlier selection or query."""
import hashlib
import json
import shutil

import pandas as pd
import pytest

from parascape.deposited_values_panel import DepositedValuesPanel
from test_deposited_values import pack  # noqa: F401


def digest(path):
    """Pin source fixture contents independently of the panel."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def sources(pack):
    """Two sources share feature labels but have different quantity semantics and source identities."""
    namespace = pack.parent / 'protein_extra'
    namespace.mkdir()
    other = namespace / 'ratio_source'
    shutil.copytree(pack, other)
    info = json.loads((other / 'manifest.json').read_text())
    info.update(study='SOURCE2', source_analysis='file.tsv :: Ratio H/L',
                assay='deposited protein-group quantities', units='Reported source ratio')
    samples = json.loads((other / 'samples.json').read_text())
    samples.update(study=info['study'], source_analysis=info['source_analysis'])
    (other / 'samples.json').write_text(json.dumps(samples))
    (other / 'source_metadata.json').write_text(json.dumps(dict(source_units=info['units'],
        reported_quantity_role='Source-reported ratio; sample assignment unreviewed')))
    info['outputs'] = {name: digest(other / name) for name in info['outputs']}
    (other / 'manifest.json').write_text(json.dumps(info))
    receipt = json.loads((other / 'validation.json').read_text())
    receipt['manifest_sha256'] = digest(other / 'manifest.json')
    (other / 'validation.json').write_text(json.dumps(receipt))
    (namespace / 'catalogue.json').write_text(json.dumps(dict(packs=[other.name])))
    return pack, other


def test_filter_changes_source_without_leaving_stale_values_exports_or_design_dialog(qtbot, sources, tmp_path):
    original, ratio = sources
    panel = DepositedValuesPanel(original.parent)
    qtbot.addWidget(panel)
    panel.features.setPlainText('001')
    panel.show_values()
    assert panel.result.units.eq('Original units').all() and panel.export_button.isEnabled()
    panel.show_column_context()
    assert panel.context_dialog is not None
    panel.source_filter.setText('rAtIo H/l')
    assert panel.pack.count() == 1 and panel.pack.currentData()['path'] == str(ratio)
    assert panel.result is None and panel.table.rowCount() == 0
    assert panel.context_dialog is None and not panel.export_button.isEnabled()
    assert panel.save(str(tmp_path / 'stale.csv')) == 0
    assert not (tmp_path / 'stale.csv').exists()
    assert '1 of 2' in panel.catalogue_status.text()
    panel.show_values()
    output = tmp_path / 'ratio.csv'
    assert panel.save(str(output)) == 4
    exported = pd.read_csv(output, dtype=str, keep_default_na=False)
    assert exported.units.eq('Reported source ratio').all()
    assert exported.study.eq('SOURCE2').all()
    assert exported.source_analysis.eq('file.tsv :: Ratio H/L').all()
    assert exported.source_token.tolist() == ['-0.000', '1.2500', '0', '-2.500']


def test_filter_searches_source_metadata_preserves_visible_selection_and_disables_empty_results(qtbot, sources, tmp_path):
    original, ratio = sources
    panel = DepositedValuesPanel(original.parent)
    qtbot.addWidget(panel)
    for query in ('source2', 'file.tsv', 'protein-group quantities', 'Reported source ratio'):
        panel.source_filter.setText(query)
        assert panel.pack.count() == 1 and panel.pack.currentData()['path'] == str(ratio)
    panel.refresh()
    assert panel.source_filter.text() == 'Reported source ratio'
    assert panel.pack.currentData()['path'] == str(ratio)
    panel.source_filter.clear()
    assert panel.pack.count() == 2 and panel.pack.currentData()['path'] == str(ratio)
    panel.source_filter.setText('no matching source')
    assert panel.pack.count() == 0 and panel.result is None
    assert not any(button.isEnabled() for button in (panel.show_button, panel.qc_button, panel.context_button, panel.export_button))
    assert panel.save_qc(str(tmp_path / 'absent-qc.csv')) == 0
    assert not (tmp_path / 'absent-qc.csv').exists()
    assert 'No source analyses match' in panel.summary.text()
    panel.source_filter.clear()
    assert panel.pack.count() == 2 and panel.table.rowCount() == 0


def test_editing_feature_query_invalidates_previous_export_until_new_values_are_shown(qtbot, sources, tmp_path):
    original, _ = sources
    panel = DepositedValuesPanel(original.parent)
    qtbot.addWidget(panel)
    panel.features.setPlainText('001')
    panel.show_values()
    assert panel.table.rowCount() == 4
    panel.features.setPlainText('CaseID')
    assert panel.result is None and panel.table.rowCount() == 0 and not panel.export_button.isEnabled()
    assert panel.save(str(tmp_path / 'stale-feature.csv')) == 0
    panel.show_values()
    assert panel.table.rowCount() == 2 and panel.result.feature_label.eq('CaseID').all()


def test_empty_installed_catalogue_keeps_source_actions_unavailable(qtbot, tmp_path):
    panel = DepositedValuesPanel(tmp_path)
    qtbot.addWidget(panel)
    assert '0 of 0' in panel.catalogue_status.text()
    assert 'No reconciled deposited measurement data is installed' in panel.summary.text()
    assert not any(button.isEnabled() for button in (panel.show_button, panel.qc_button, panel.context_button, panel.export_button))
