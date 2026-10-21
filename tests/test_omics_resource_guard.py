"""Check scheduler protection with fake processes; never alter real worker affinity."""
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from scripts import omics_resource_guard as G


class Guard(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.proc = self.root / 'proc'
        self.proc.mkdir()

    def process(self, pid, command):
        process = self.proc / str(pid)
        process.mkdir()
        (process / 'cmdline').write_bytes(command)
        return process

    def test_tutorial_process_is_detected_but_shell_text_is_not(self):
        self.process(1, b'/bin/bash\0-c\0python run_native_capture.py\0')
        self.process(2, b'/usr/bin/python3\0omics_resource_guard.py\0')
        self.assertFalse(G.tutorial_active(self.proc))
        self.process(3, b'/some/env/bin/python\0/tools/tutorials/run_native_capture.py\0')
        self.assertTrue(G.tutorial_active(self.proc))

    def test_validation_worker_keeps_tutorial_protection(self):
        self.process(1, b'/env/bin/python\0/capture_clip_gate.py\0')
        self.assertTrue(G.tutorial_active(self.proc))

    def test_protection_is_immediate_and_boost_waits_for_quiet(self):
        mode, quiet = G.select_mode('', 100, None, 30)
        self.assertEqual((mode, quiet), ('protected', 100))
        self.assertEqual(G.select_mode('', 130, quiet, 30), ('fast', 100))
        self.assertEqual(G.select_mode('tutorial', 131, quiet, 30), ('protected', None))
        self.assertEqual(G.select_mode('', 132, None, 30), ('protected', 132))

    def test_all_worker_threads_change_but_unrelated_process_is_untouched(self):
        cgroup = self.root / 'cgroup'
        cgroup.mkdir()
        (cgroup / 'cgroup.procs').write_text('11\n22\n')
        for pid, tids in ((11, (11, 12)), (22, (22, 23)), (33, (33,))):
            p = self.process(pid, b'/usr/bin/python\0')
            (p / 'task').mkdir()
            for tid in tids:
                (p / 'task' / str(tid)).mkdir()
        with mock.patch.object(G.os, 'sched_getaffinity', return_value={0, 1}), \
                mock.patch.object(G.os, 'sched_setaffinity') as setter:
            self.assertEqual(G.apply_affinity(cgroup, {10, 11}, self.proc), 4)
            self.assertEqual({call.args[0] for call in setter.call_args_list}, {11, 12, 22, 23})

    def test_worker_absent_or_exiting_is_harmless(self):
        self.assertEqual(G.apply_affinity(self.root / 'absent', {10}, self.proc), 0)
        cgroup = self.root / 'cgroup'
        cgroup.mkdir()
        (cgroup / 'cgroup.procs').write_text('11\n')
        p = self.process(11, b'/usr/bin/python\0')
        (p / 'task' / '11').mkdir(parents=True)
        with mock.patch.object(G.os, 'sched_getaffinity', side_effect=ProcessLookupError):
            self.assertEqual(G.apply_affinity(cgroup, {10}, self.proc), 0)


if __name__ == '__main__':
    unittest.main()
