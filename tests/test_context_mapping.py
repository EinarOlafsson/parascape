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
    assert P.propose('strain', 'N/A')[0] == 'not_available'
    assert P.propose('gating strategy', 'CD45+') == ('not_context', '', 'first-pass tag review: "gating strategy"')
    assert P.propose('favourite colour', 'blue')[0] == 'review'                # unknown tags stay for review


def test_first_pass_tag_review_and_units_carried_by_the_tag():
    assert P.propose('develpomental stage', 'bradyzoite')[0] == 'parasite_stage'   # misspelt stage tag, stage rule
    assert P.propose('develpomental stage', 'adult')[0] == 'review'
    assert P.propose('toxoplasma parasite strain', 'Pru dku80')[0] == 'parasite_strain'
    assert P.propose('mouse id', 'M12')[0] == 'unit_or_replicate_identifier'
    field, hours, rule = P.propose('days post-infection', '7')
    assert (field, hours) == ('time_declared_hours', '168') and 'after infection' in rule
    assert P.propose('days post-infection', 'seven')[0] != 'time_declared_hours'


def test_value_level_first_pass():
    assert P.propose('treatment', 'none')[0] == 'condition_declared'            # "none" = untreated, kept
    assert P.propose('parasite strain', 'n/a')[0] == 'not_available'
    assert P.propose('strain', 'RHKU80_AP2-XII-1-KD-HA')[0] == 'parasite_strain'
    assert P.propose('strain', 'RH::TIR1::METTL3-AID')[0] == 'parasite_strain'
    assert P.propose('cell line/strain', 'TgMe49')[0] == 'parasite_strain'
    assert P.propose('strain', 'NOD/ShiLtJ')[0] == 'host_strain'
    assert P.propose('strain', 'C57/BL6 x B6.C-H2d/bByJ')[0] == 'host_strain'
    assert P.propose('type', 'II') == ('parasite_strain', 'Type II', 'type tag with a T. gondii lineage numeral')
    assert P.propose('developmental stage', 'in vitro bradyzoite pH 8')[0] == 'parasite_stage'
    assert P.propose('developmental stage', 'myoblast')[0] == 'review'
    field, hours, rule = P.propose('developmentalstage', '3 days post infection')
    assert (field, hours) == ('time_declared_hours', '72') and 'after infection' in rule
    field, hours, rule = P.propose('time point', 'Day10')
    assert (field, hours) == ('time_declared_hours', '240') and 'needs review' in rule
    assert P.propose('time', 'post inoculation day 25')[1] == '600'
    assert P.propose('growth time', 'three days')[1] == '72'
    assert P.propose('mouse comments', 'DOB 5/1/22')[0] == 'free_text_comment'


def test_every_reviewed_time_reference_cites_its_evidence():
    assert len(P.STUDY_TIME_REFERENCE) == 24
    for study, (event, evidence) in P.STUDY_TIME_REFERENCE.items():
        assert study.startswith('GSE') and event and (evidence.startswith(('GEO', 'PMC', 'PMID')))
    assert P.STUDY_TIME_REFERENCE['GSE303610'][0] == 'auxin-induced knockdown'      # "24h" there is not infection
