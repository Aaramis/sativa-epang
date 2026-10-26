#!/usr/bin/env python3
import json
import os
import stat
import sys
import tempfile
import unittest

lib_path = os.path.abspath('..')
sys.path.append(lib_path)

from epac import epang_l1o
from epac.epang_l1o import MANIFEST_VERSION, l1o_exit_code, place_l1o_tasks


class PlaceExitCodeTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)
        self.taskdir = os.path.join(self.dir.name, "tasks")
        os.makedirs(os.path.join(self.taskdir, "fold_000"))
        with open(os.path.join(self.taskdir, "manifest.json"), "w") as handle:
            json.dump({"manifest_version": MANIFEST_VERSION, "model": "GTR", "epang_args": [],
                       "folds": [{"id": 0, "dir": "fold_000"}]}, handle)
        self.saved = epang_l1o.EPANG
        self.addCleanup(setattr, epang_l1o, "EPANG", self.saved)

    def fake_epang(self, body):
        path = os.path.join(self.dir.name, "epa-ng")
        with open(path, "w") as handle:
            handle.write("#!/bin/sh\n" + body + "\n")
        os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)
        epang_l1o.EPANG = path

    def place(self):
        failures = []
        placed = place_l1o_tasks(self.taskdir, failures=failures)
        return placed, failures

    def test_placed_fold_is_not_a_failure(self):
        self.fake_epang('touch "$(dirname "$0")/tasks/fold_000/epa_result.jplace"')
        placed, failures = self.place()
        self.assertEqual((placed, failures), (1, []))
        self.assertEqual(l1o_exit_code(placed, 1, failures), 0)

    def test_killed_epang_exits_137(self):
        self.fake_epang("kill -9 $$")
        placed, failures = self.place()
        self.assertEqual((placed, failures), (0, [(0, -9)]))
        self.assertEqual(l1o_exit_code(placed, 1, failures), 137)

    def test_other_failure_exits_1(self):
        self.fake_epang("exit 3")
        placed, failures = self.place()
        self.assertEqual((placed, failures), (0, [(0, 3)]))
        self.assertEqual(l1o_exit_code(placed, 1, failures), 1)

    def test_a_kill_among_other_failures_exits_137(self):
        self.assertEqual(l1o_exit_code(0, 3, [(0, 3), (1, -9), (2, 1)]), 137)


if __name__ == '__main__':
    unittest.main()
