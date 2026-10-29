"""Sample/run links must not manufacture independent replication or control pairing."""
import pandas as pd

from parascape.rna_design import annotate, independent_runs


def test_same_biosample_links_runs_without_claiming_independent_replicates():
    rows = [{'study': 's', 'run_accession': 'r1', 'sample_accession': 'b1', 'sample_title': 'control_1', 'sample_alias': 'WT_Rep1'},
            {'study': 's', 'run_accession': 'r2', 'sample_accession': 'b1', 'sample_title': 'control_1', 'sample_alias': 'WT_Rep1'},
            {'study': 't', 'run_accession': 'r1', 'sample_accession': 'b1', 'sample_title': 'control_1'}]
    d = annotate(rows)
    assert set(d.linked_runs) == {'r1, r2'}
    assert 'Rep1' in d.replicate_label.iloc[0]
    assert d.control_role.iloc[0] == 'control (source label)'
    assert d.shared_across_studies.iloc[0]
    assert not independent_runs(d[d.study == 's'], ['r1', 'r2'])
    assert independent_runs(d[d.study == 's'], ['r1'])
    assert d.replicate_type.str.contains('not specified').all()
    assert d.control_link.str.contains('No verified pairing').all()


def test_shared_strain_attribute_does_not_merge_different_title_groups():
    rows = [{'study': 's', 'run_accession': 'r1', 'sample_accession': 'b1', 'sample_title': 'wild type control_1'},
            {'study': 's', 'run_accession': 'r2', 'sample_accession': 'b2', 'sample_title': 'mutant_1'}]
    metadata = {b: {'attributes': {'strain': 'example'}, 'status': 'retrieved'} for b in ['b1', 'b2']}
    d = annotate(rows, metadata)
    assert d.condition_group.nunique() == 2


def test_study_browser_exposes_sample_relationships_and_exports(qtbot, tmp_path):
    from parascape.rna_design_panel import StudyDesignPanel
    d = annotate([{'study': 's', 'run_accession': 'r1', 'sample_accession': 'b1',
                   'sample_title': 'control_1', 'processing_status': 'pending'}])
    d.to_csv(tmp_path / 'sample_design.tsv', sep='\t', index=False)
    panel = StudyDesignPanel(str(tmp_path))
    qtbot.addWidget(panel)
    assert panel.tree.topLevelItemCount() == 1
    top = panel.tree.topLevelItem(0)
    assert top.text(1) == '0/1 runs verified'
    panel.search.setText('control')
    assert not top.isHidden()
    assert panel.save(str(tmp_path / 'out.csv')) == 1
