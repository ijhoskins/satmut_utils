#!/usr/bin/env/python
""" Tests for reading NGS reads by ReadEnumerator
"""

import unittest
import read_enumerators


class TestReadGenerator(unittest.TestCase):
    """ Tests for Read Generator
    """

    def setUpClass(cls):

        with read_enumerators.ReadEnumerator(r1=) as rg:
            cls.rg = rg


    def tearDownClass(cls):
        pass

    def test_filter_and_count_reads(self):
        pass

    def test_read_filter(self):
        pass

    def test_num_mismatches(self):
        pass

    def test_num_indels(self):
        pass

    def test_slop_target_features(self):
        pass

    def test_mapped_counts(self):
        pass

    def test_unmapped_counts(self):
        pass

    def test_mapped_proportion(self):
        pass

    def test_unique_filtered_counts(self):
        pass

    def test_unique_unfiltered_counts(self):
        pass

    def test_unfiltered_proportion(self):
        pass

    def test_get_qnames_from_bedtool(self):
        pass

    def test_get_on_target(self):
        pass

    def test_get_primer_counts(self):
        pass

    def test_get_target_counts(self):
        pass

    def test_check_files(self):
        pass

