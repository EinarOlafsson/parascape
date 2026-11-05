"""Independent scalar-only audit of the frozen EC recovery benchmark contract.

These checks read immutable metadata and bounded held-out JSON. They never import
the application, fit an estimator, decode the large model state or acquire data.
They establish source-recovery bookkeeping, not biological truth or independence.
"""
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / 'results/functional_ec_knn_pilot_2026_10_08_v2'


def read(name):
    """Read a bounded immutable metadata file and reject nonfinite JSON constants."""
    path = PILOT / name
    assert path.stat().st_size < 300_000, 'This audit only reads bounded metadata'
    def invalid(value):
        raise ValueError('Nonfinite JSON constant: ' + value)
    return json.loads(path.read_text(), parse_constant=invalid)


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(',', ':')).encode()


def digest_file(path):
    """Stream even model-state hashes without decoding or retaining those payloads."""
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(64 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def members(value):
    """Interpret only explicit canonical nonempty software profile labels."""
    values = json.loads(value)
    assert values and isinstance(values, list)
    assert values == sorted(set(values))
    assert set(values) <= set('1234567')
    return set(values)


def test_every_frozen_output_and_artifact_payload_retains_its_byte_identity():
    manifest = read('manifest.json')
    for relative, expected in manifest['output_sha256'].items():
        path = PILOT / relative
        assert path.is_file() and not path.is_symlink()
        assert digest_file(path) == expected, relative
    artifact = read('held_out/manifest.json')
    for relative, receipt in artifact['files'].items():
        path = PILOT / 'held_out' / relative
        assert path.stat().st_size == receipt['bytes']
        assert digest_file(path) == receipt['sha256']
    body = {key: value for key, value in artifact.items() if key != 'identity'}
    assert hashlib.sha256(canonical(body)).hexdigest() == artifact['identity']
    assert hashlib.sha256(canonical(artifact['spec'])).hexdigest() == artifact['key']


def test_unknown_universe_entities_are_excluded_without_creating_negative_truth():
    entry = read('registry.json')['entries'][0]
    universe = read('universe.json')
    assignments = read('split.json')['split']['assignments']
    packed = bytes.fromhex(entry['eligibility_hex'])
    bits = [bool(byte & (1 << shift)) for byte in packed for shift in range(7, -1, -1)]
    mask = bits[:len(universe)]
    eligible = [entity for entity, selected in zip(universe, mask) if selected]
    assert len(universe) == len(set(universe)) == entry['stored_population']
    assert not any(bits[len(universe):])
    assert len(eligible) == entry['eligible_population']
    assert set(eligible) == {row['entity'] for row in assignments}
    cohort = json.dumps(list(zip(universe, mask)), separators=(',', ':')).encode()
    assert hashlib.sha256(cohort).hexdigest() == entry['cohort_sha256']
    assert entry['negative_values'] == []
    assert entry['measured_population'] is None
    assert entry['status'] == 'candidate' and entry['evidence_grade'] == 'unresolved'


def test_homology_groups_and_fit_entities_do_not_cross_nested_roles():
    document = read('split.json')
    assignments = document['split']['assignments']
    roles = {}
    entities = set()
    for row in assignments:
        assert row['entity'] not in entities
        entities.add(row['entity'])
        assert row['role'] in {'train', 'tune', 'calibration', 'test'}
        assert row['group'] not in roles or roles[row['group']] == row['role']
        roles[row['group']] = row['role']
    train = [row['entity'] for row in assignments if row['role'] == 'train']
    test = [row['entity'] for row in assignments if row['role'] == 'test']
    spec = read('held_out/manifest.json')['spec']
    assert spec['fit_entities'] == train == document['exclusions']['fit_entities']
    assert spec['entity_order'] == test
    assert spec['fit_role'] == 'train' and spec['evaluation_partition'] == 'test'
    assert spec['confidence_kind'] == 'method_support'
    assert not (set(train) & set(test))
    summary = read('summary.json')
    assert Counter(row['role'] for row in assignments) == {role: summary[role] for role in ('train', 'tune', 'calibration', 'test')}


def test_source_exclusions_and_truth_dependency_are_pinned_in_the_artifact():
    document = read('split.json')
    entry = read('registry.json')['entries'][0]
    exclusions = document['exclusions']
    assert exclusions['columns'] == entry['exclusions']
    assert {'ec_number', 'ec_number_orthology', 'ec_major_classes', 'has_ec',
            'interpro_id', 'interpro_ids', 'interpro_desc', 'pfam_id', 'pfam_ids',
            'orthogroup', 'paralog_number'} <= set(exclusions['columns'])
    dependencies = read('held_out/manifest.json')['spec']['dependencies']
    by_kind = {(item['kind'], item['name']): item['sha256'] for item in dependencies}
    assert by_kind['exclusions', 'training_source_closure'] == hashlib.sha256(canonical(exclusions)).hexdigest()
    assert by_kind['split', 'nested'] == exclusions['split_identity']
    manifest = read('manifest.json')
    profile_receipt = next(value for path, value in manifest['input_sha256'].items() if path.endswith('_ec_number_profiles.parquet'))
    assert by_kind['truth', 'complete_major_class_profiles'] == profile_receipt


def test_abstentions_and_unsupported_profiles_remain_in_accuracy_denominators():
    rows = read('held_out/rows.json')
    card = read('held_out/card.json')
    spec = read('held_out/manifest.json')['spec']
    assert [row['entity'] for row in rows] == spec['entity_order']
    answered = correct = 0
    for row in rows:
        members(row['truth'])
        if row['prediction'] is not None:
            members(row['prediction'])
            answered += 1
            correct += row['truth'] == row['prediction']
        assert row['abstained'] == (row['prediction'] is None)
        assert row['abstained'] == (row['support'] < .3)
        assert math.isfinite(row['support']) and 0 <= row['support'] <= 1
        assert row['calibrated_confidence'] is None
    expected = {'eligible': len(rows), 'answered': answered, 'correct': correct,
                'abstained': len(rows) - answered, 'wrong': answered - correct,
                'unique_biological_entities': len(rows)}
    assert card['counts'] == expected
    assert card['metrics']['accuracy'] == correct / len(rows)
    assert card['metrics']['precision_of_calls'] == correct / answered
    assert card['metrics']['coverage'] == answered / len(rows)
    scores = read('held_out/class_scores.json')
    unseen = Counter(row['truth'] for row in rows if row['truth'] not in scores['classes'])
    assert dict(unseen) == read('summary.json')['unsupported_test_profiles']
    assert len(scores['scores']) == len(rows)
    for profile in scores['classes']:
        members(profile)


@pytest.mark.parametrize('major', list('1234567'))
def test_member_class_cards_preserve_overlaps_unknown_biology_and_abstentions(major):
    rows = read('held_out/rows.json')
    actual = [members(row['truth']) for row in rows]
    predicted = [members(row['prediction']) if row['prediction'] is not None else set() for row in rows]
    card = next(card for card in read('held_out/major_class_cards.json') if card['class'] == major)
    tp = sum(major in truth and major in call for truth, call in zip(actual, predicted))
    fp = sum(major not in truth and major in call for truth, call in zip(actual, predicted))
    fn = sum(major in truth and major not in call for truth, call in zip(actual, predicted))
    assert (card['true_positive'], card['false_positive'], card['false_negative']) == (tp, fp, fn)
    assert card['known_positive_genes'] == tp + fn
    assert card['eligible'] == len(rows)
    assert card['precision'] == (tp / (tp + fp) if tp + fp else None)
    assert card['recall'] == (tp / (tp + fn) if tp + fn else None)
    assert card['coverage'] == sum(row['prediction'] is not None for row in rows) / len(rows)
    assert card['biological_precision'] is card['biological_recall'] is card['calibrated_confidence'] is None


def test_baseline_cards_have_identical_frozen_scope_and_test_population():
    card = read('held_out/card.json')
    baselines = read('held_out/baseline_cards.json')
    assert set(baselines) == {'majority', 'prevalence_call'}
    for baseline in baselines.values():
        assert baseline['scope'] == card['scope']
        assert baseline['counts']['eligible'] == card['counts']['eligible']
        assert baseline['counts']['abstained'] == 0
    assert read('summary.json')['biological_benchmark_admitted'] is False
