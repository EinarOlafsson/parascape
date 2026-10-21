"""Unattended orchestration tests with synthetic files and fake service control."""
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import shutil
import unittest
from unittest import mock

from scripts import omics_queue_supervisor as Q


class FakeServices:
    def __init__(self, active=True):
        self.states = {unit: 'active' if active else 'inactive' for unit in Q.SERVICES}
        self.calls = []

    def state(self, unit):
        return self.states.get(unit, 'inactive')

    def command(self, action, *units):
        self.calls.append((action, units))
        for unit in units:
            self.states[unit] = 'active' if action == 'start' else 'inactive'


class Supervisor(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.state = self.root / 'automation'
        self.state.mkdir()
        self.outputs = self.root / 'outputs'
        self.plan = self.state / 'plan.jsonl'
        self.row = {'run_accession': 'SRR1', 'study': 'fixture', 'references': ['reference']}
        self.plan.write_text(json.dumps(self.row) + '\n')
        self.config = {'plan_sha256': hashlib.sha256(self.plan.read_bytes()).hexdigest()}
        Q.save(self.state / 'config.json', self.config)

    def complete_run(self):
        folder = self.outputs / 'fixture' / 'SRR1'
        folder.mkdir(parents=True)
        Q.save(folder / 'done.json', {**self.row, 'salmon': 'fixture', 'kallisto': 'fixture',
                                     'finished': 'fixture'})
        for tool, filename in (('salmon', 'quant.sf.gz'), ('kallisto', 'abundance.tsv.gz')):
            output = folder / f'{tool}_reference' / filename
            output.parent.mkdir()
            with gzip.open(output, 'wb') as stream:
                stream.write(Q.HEADERS[tool] + b'\nrow\t1\t1\t0\t0\n')
        return folder

    def test_receipt_alone_is_not_completion_and_is_preserved_for_review(self):
        folder = self.outputs / 'fixture' / 'SRR1'
        folder.mkdir(parents=True)
        Q.save(folder / 'done.json', {**self.row, 'salmon': 'x', 'kallisto': 'x', 'finished': 'x'})
        accepted, pending = Q.audit([self.row], self.outputs, {}, None)
        self.assertFalse(accepted)
        self.assertEqual(len(pending), 1)
        self.assertFalse((folder / 'done.json').exists())
        self.assertEqual(len(list(folder.glob('done.invalid-*.json'))), 1)

    def test_truncated_gzip_is_requeued_and_raw_output_is_preserved(self):
        folder = self.complete_run()
        output = folder / 'salmon_reference' / 'quant.sf.gz'
        output.write_bytes(output.read_bytes()[:-8])
        accepted, pending = Q.audit([self.row], self.outputs, {}, None)
        self.assertFalse(accepted)
        self.assertEqual(len(pending), 1)
        self.assertTrue(output.exists())
        self.assertFalse((folder / 'done.json').exists())

    def test_validated_results_reuse_hashes_until_a_file_changes(self):
        folder = self.complete_run()
        first = Q.inspect_run(self.row, self.outputs, None)
        self.assertTrue(first['valid'])
        with mock.patch.object(Q, 'verify_table', side_effect=AssertionError('unnecessary reread')):
            second = Q.inspect_run(self.row, self.outputs, first['evidence'])
        self.assertTrue(second['reused'])
        output = folder / 'salmon_reference' / 'quant.sf.gz'
        output.write_bytes(b'corrupt')
        self.assertFalse(Q.inspect_run(self.row, self.outputs, first['evidence'])['valid'])

    def test_active_run_is_not_quarantined_during_final_writes(self):
        folder = self.complete_run()
        (folder / 'done.json').write_text('{')
        accepted, pending = Q.audit([self.row], self.outputs, {}, 'SRR1')
        self.assertFalse(accepted)
        self.assertEqual(len(pending), 1)
        self.assertTrue((folder / 'done.json').exists())

    def test_successful_but_incomplete_worker_exit_is_restarted(self):
        services = FakeServices(active=False)
        status = Q.tick(self.plan, self.outputs, self.state, services)
        self.assertEqual(status['remaining'], 1)
        self.assertFalse((self.state / 'completion.json').exists())
        self.assertIn(('start', (Q.SERVICES[2],)), services.calls)

    def test_verified_completion_stops_owned_workers_and_timer(self):
        self.complete_run()
        services = FakeServices()
        status = Q.tick(self.plan, self.outputs, self.state, services)
        self.assertEqual(status['phase'], 'complete')
        receipt = Q.read_json(self.state / 'completion.json')
        self.assertEqual(receipt['verified_complete'], 1)
        self.assertEqual(receipt['plan_sha256'], self.config['plan_sha256'])
        self.assertEqual(receipt['verified_manifest_sha256'], Q.digest(self.state / 'verified.json'))
        self.assertIn(('stop', Q.SERVICES), services.calls)
        self.assertIn(('stop', (Q.TIMER,)), services.calls)

    def test_false_completion_is_removed_when_results_go_missing(self):
        folder = self.complete_run()
        Q.tick(self.plan, self.outputs, self.state, FakeServices())
        (folder / 'kallisto_reference' / 'abundance.tsv.gz').unlink()
        status = Q.tick(self.plan, self.outputs, self.state, FakeServices())
        self.assertEqual(status['remaining'], 1)
        self.assertFalse((self.state / 'completion.json').exists())

    def test_finished_queue_does_not_start_workers_again_after_reboot(self):
        self.complete_run()
        Q.tick(self.plan, self.outputs, self.state, FakeServices())
        services = FakeServices(active=False)
        status = Q.tick(self.plan, self.outputs, self.state, services)
        self.assertEqual(status['phase'], 'complete')
        self.assertFalse(any(action == 'start' for action, units in services.calls))

    def test_retry_backoff_is_bounded_but_has_no_retry_count_limit(self):
        services = FakeServices()
        services.states[Q.SERVICES[2]] = 'inactive'
        first = Q.recover(services, {}, 0, 0)
        self.assertEqual(first['retry_at'], 60)
        services.states[Q.SERVICES[2]] = 'inactive'
        services.calls.clear()
        second = Q.recover(services, first, 0, 10)
        self.assertEqual(second['worker'], 'retry_wait')
        self.assertFalse(services.calls)
        for _ in range(20):
            services.states[Q.SERVICES[2]] = 'inactive'
            now = second['retry_at']
            second = Q.recover(services, second, 0, now)
            self.assertLessEqual(second['retry_at'] - now, 3600)
        self.assertEqual(second['attempt'], 21)
        self.assertEqual(second['worker'], 'restarted')

    def test_pause_marker_controls_no_services(self):
        (self.state / 'PAUSE').touch()
        services = FakeServices(active=False)
        self.assertEqual(Q.tick(self.plan, self.outputs, self.state, services)['phase'], 'paused')
        self.assertFalse(services.calls)

    def test_empty_changed_or_duplicate_plan_cannot_report_completion(self):
        for raw in (b'', self.plan.read_bytes() * 2):
            self.plan.write_bytes(raw)
            with self.assertRaises(ValueError):
                Q.load_plan(self.plan, hashlib.sha256(raw).hexdigest())
        self.plan.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            Q.load_plan(self.plan, self.config['plan_sha256'])

    def test_service_backend_rejects_unrelated_process_control(self):
        with self.assertRaises(ValueError):
            Q.Services().command('stop', 'unrelated.service')

    def test_shared_run_in_two_studies_requires_both_outputs(self):
        folder = self.complete_run()
        second_row = {**self.row, 'study': 'alias'}
        rows = [self.row, second_row]
        self.plan.write_text(''.join(json.dumps(row) + '\n' for row in rows))
        self.assertEqual(len(Q.load_plan(self.plan, Q.digest(self.plan))), 2)
        accepted, pending = Q.audit(rows, self.outputs, {}, None)
        self.assertEqual(len(accepted), 1)
        self.assertEqual(len(pending), 1)
        second_folder = self.outputs / 'alias' / 'SRR1'
        shutil.copytree(folder, second_folder)
        receipt = Q.read_json(second_folder / 'done.json')
        Q.save(second_folder / 'done.json', {**receipt, 'study': 'alias'})
        accepted, pending = Q.audit(rows, self.outputs, accepted, None)
        self.assertEqual(set(accepted), {'fixture/SRR1', 'alias/SRR1'})
        self.assertFalse(pending)


if __name__ == '__main__':
    unittest.main()
