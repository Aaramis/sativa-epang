#!/usr/bin/env python3
import os
import sys
import unittest

lib_path = os.path.abspath('..')
sys.path.append(lib_path)

from epac.epang_l1o import parse_fold_selection


class FoldSelectionTests(unittest.TestCase):
    IDS = list(range(25))

    def test_everything_by_default(self):
        self.assertEqual(parse_fold_selection(None, self.IDS), self.IDS)
        self.assertEqual(parse_fold_selection("", self.IDS), self.IDS)

    def test_named_folds(self):
        self.assertEqual(parse_fold_selection("7", self.IDS), [7])
        self.assertEqual(parse_fold_selection("0-4", self.IDS), [0, 1, 2, 3, 4])
        self.assertEqual(parse_fold_selection("0,7,24", self.IDS), [0, 7, 24])
        self.assertEqual(parse_fold_selection("0-2,10,24", self.IDS), [0, 1, 2, 10, 24])

    def test_shards_cover_everything_once(self):
        shards = [parse_fold_selection("%d/5" % i, self.IDS) for i in range(5)]
        self.assertEqual(sorted(sum(shards, [])), self.IDS)
        self.assertEqual(shards[3], [3, 8, 13, 18, 23])

    def test_bad_selections_are_refused(self):
        for spec in ("99", "0-30", "5-2", "abc", "-1", "10/3", "3/0"):
            self.assertRaises(ValueError, parse_fold_selection, spec, self.IDS)


if __name__ == '__main__':
    unittest.main()
