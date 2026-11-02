"""The literal sample registry keeps author annotations verbatim and records SuperSeries membership."""
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('build_sample_registry', ROOT / 'scripts' / 'build_sample_registry.py')
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)


def series(tmp, name, tags, reconciled=True):
    folder = tmp / name; folder.mkdir()
    (folder / 'validation.json').write_text(json.dumps(dict(exact_sample_ids_reconciled=reconciled)))
    (folder / 'geo-samples.json').write_text(json.dumps(dict(samples=[dict(
        geo_accession='GSM1', title='  L1Toxo\n', sample_type='SRA', platform_refs=['GPL1'],
        relations=[dict(type='BioSample', target='https://www.ncbi.nlm.nih.gov/biosample/SAMN1'),
                   dict(type='SRA', target='https://www.ncbi.nlm.nih.gov/sra?term=SRX2')],
        channels=[dict(position='1', source=' rat peritoneal cells ', organisms=['Rattus norvegicus'],
                       characteristics=[dict(tag=t, value=v) for t, v in tags])])])))
    return folder


def test_values_stay_verbatim_apart_from_surrounding_whitespace(tmp_path):
    samples, annotations, info = R.study_records(series(tmp_path, 'GSE1', [('strain', '\nLewis  rat\n   '), ('Strain', 'RH')]))
    assert info['status'] == 'included' and samples[0]['biosample'] == 'SAMN1' and samples[0]['sra_experiment'] == 'SRX2'
    values = {(a['tag'], a['value']) for a in annotations}
    assert ('strain', 'Lewis  rat') in values and ('Strain', 'RH') in values          # internal spacing and case kept
    assert ('organism', 'Rattus norvegicus') in values and ('source_name', 'rat peritoneal cells') in values


def test_unreconciled_series_are_reported_not_used(tmp_path):
    samples, annotations, info = R.study_records(series(tmp_path, 'GSE2', [('strain', 'RH')], reconciled=False))
    assert samples == [] and annotations == [] and info['status'] == 'sample_ids_not_reconciled'
