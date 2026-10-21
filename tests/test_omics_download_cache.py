"""Transport-only tests with synthetic bytes and a local HTTP server (no Qt or datasets)."""
import contextlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib.util
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C = load('omics_download_cache')
Q = load('omics_quant')
PAYLOAD = b'fixture data\n' * 5000


@contextlib.contextmanager
def server(*, ranges=True, truncate=False):
    requests = []

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            value = self.headers.get('Range')
            requests.append(value)
            offset = int(value.split('=')[1].split('-')[0]) if value and ranges else 0
            self.send_response(206 if value and ranges else 200)
            self.send_header('Content-Length', str(len(PAYLOAD) - offset))
            if value and ranges:
                self.send_header('Content-Range', f'bytes {offset}-{len(PAYLOAD)-1}/{len(PAYLOAD)}')
            self.end_headers()
            self.wfile.write(PAYLOAD[offset:-7] if truncate else PAYLOAD[offset:])

        def log_message(self, *_):
            pass

    httpd = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    try:
        yield f'http://127.0.0.1:{httpd.server_port}/fixture', requests
    finally:
        httpd.shutdown()
        httpd.server_close()
        thread.join()


class Downloads(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.archive = self.folder / 'SRR1.sra'

    def transfer(self, url, **kwargs):
        return C.download(url, self.archive, len(PAYLOAD), rate=10**9,
                          reserve=0, budget=len(PAYLOAD), guard=lambda: '', **kwargs)

    def test_resume_verifies_exact_bytes_and_receipt(self):
        self.archive.with_suffix('.part').write_bytes(PAYLOAD[:173])
        with server() as (url, requests):
            self.transfer(url)
            self.assertEqual(requests, ['bytes=173-'])
            self.assertEqual(self.archive.read_bytes(), PAYLOAD)
            self.assertEqual(C.cached_archive(self.folder, 'SRR1', url, len(PAYLOAD)), self.archive)
            self.assertFalse(self.archive.with_suffix('.part').exists())

    def test_server_ignoring_range_restarts_without_duplicate_bytes(self):
        self.archive.with_suffix('.part').write_bytes(PAYLOAD[:173])
        with server(ranges=False) as (url, _):
            self.transfer(url)
        self.assertEqual(self.archive.read_bytes(), PAYLOAD)

    def test_truncation_keeps_only_resumable_partial(self):
        with server(truncate=True) as (url, _):
            with self.assertRaises((RuntimeError, OSError)):
                self.transfer(url)
        self.assertFalse(self.archive.exists())
        self.assertFalse(self.archive.with_suffix('.json').exists())
        self.assertGreater(self.archive.with_suffix('.part').stat().st_size, 0)

    def test_reserve_and_budget_prevent_any_request(self):
        with server() as (url, requests):
            with self.assertRaisesRegex(RuntimeError, 'reserve'):
                C.download(url, self.archive, len(PAYLOAD), rate=1, reserve=10**20, budget=10**20)
            self.archive.with_suffix('.part').write_bytes(PAYLOAD[:173])
            (self.folder / 'SRR2.part').write_bytes(b'x')
            with self.assertRaisesRegex(RuntimeError, 'capacity'):
                self.transfer(url)
            self.assertEqual(requests, [])

    def test_checksum_corruption_is_rejected(self):
        with server() as (url, _):
            self.transfer(url)
            self.archive.write_bytes(b'x' * len(PAYLOAD))
            self.assertIsNone(C.cached_archive(self.folder, 'SRR1', url, len(PAYLOAD)))

    def test_consumer_waits_for_atomic_completion(self):
        with server() as (url, _):
            self.transfer(url)
            result = []
            started = threading.Event()

            def consume():
                started.set()
                result.append(C.cached_archive(self.folder, 'SRR1', url, len(PAYLOAD)))

            with C.archive_lock(self.archive):
                thread = threading.Thread(target=consume)
                thread.start()
                started.wait(1)
                self.assertEqual(result, [])
            thread.join(2)
            self.assertEqual(result, [self.archive])

    def test_rate_and_pressure_waits_are_enforced(self):
        waits = []
        with server() as (url, _), mock.patch.object(C.time, 'sleep', side_effect=waits.append):
            C.download(url, self.archive, len(PAYLOAD), rate=1000, reserve=0,
                       budget=len(PAYLOAD), guard=mock.Mock(side_effect=['busy', '']))
        self.assertEqual(waits[0], 5)
        self.assertGreater(waits[-1], 59)

    def test_worker_uses_cache_and_cleans_it_after_successful_expansion(self):
        scratch = self.folder / 'scratch'
        scratch.mkdir()
        with server() as (url, _):
            self.transfer(url)
            commands = []

            def fake_run(command, **kwargs):
                commands.append(command)
                if 'fasterq-dump' in command[0]:
                    for mate in (1, 2):
                        (scratch / f'SRR1_{mate}.fastq').write_text('synthetic placeholder')

            with mock.patch.multiple(Q, DOWNLOAD_CACHE=str(self.folder), SCRATCH=str(scratch), SRA_ODP=url), \
                    mock.patch.object(Q, '_head_size', return_value=len(PAYLOAD)), \
                    mock.patch.object(Q, '_free_gb', return_value=1000), \
                    mock.patch.object(Q, '_cache_module', return_value=C), \
                    mock.patch.object(Q, 'fetch', side_effect=AssertionError('duplicate download')), \
                    mock.patch.object(Q, 'log'), mock.patch.object(Q, 'run', side_effect=fake_run):
                reads = Q._reads({'run_accession': 'SRR1'})
            self.assertEqual(len(reads), 2)
            self.assertEqual(commands[0][-1], str(self.archive))
            self.assertFalse(self.archive.exists())
            self.assertFalse(self.archive.with_suffix('.json').exists())


if __name__ == '__main__':
    unittest.main()
