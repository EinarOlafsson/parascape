"""Context proposals come only from narrow named rules; everything else stays for review."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('propose_context_mapping', ROOT / 'scripts' / 'propose_context_mapping.py')
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)


def test_the_ambiguous_strain_tag_is_split_by_value_not_by_name():
    assert P.propose('strain', 'ME49')[0] == 'parasite_strain'
    assert P.propose('strain', 'RH strain')[0] == 'parasite_strain'
    assert P.propose('strain', 'C57BL/6J')[0] == 'host_strain'
    assert P.propose('strain', 'Lewis')[0] == 'host_strain'
    assert P.propose('strain', 'wild type') == ('review', '', 'strain value matches neither list')


def test_stages_times_and_unknowns():
    assert P.propose('cell type', 'Tachyzoites')[0] == 'parasite_stage'
    assert P.propose('cell type', 'HFF')[0] == 'host_cell_or_tissue'
    assert P.propose('developmental stage', 'adult')[0] == 'review'
    field, hours, rule = P.propose('timepoint', '9 dpi')
    assert (field, hours) == ('time_declared_hours', '216') and 'after infection' in rule
    field, hours, rule = P.propose('time', '6h')
    assert hours == '6' and 'needs review' in rule                    # start event not stated
    assert P.propose('time', 'early')[0] == 'review'
    assert P.propose('strain', 'N/A')[0] == 'review'
    assert P.propose('gating strategy', 'CD45+')[0] == 'review'
