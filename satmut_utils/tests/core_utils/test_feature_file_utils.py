#!/usr/bin/env/python
""" Tests for core_utils.feature_file_utils."""

import collections
import pybedtools
import tempfile
import unittest

from analysis.seq_utils import Strand, COORD_FORMAT_STRAND
import core_utils.feature_file_utils as ffu
from tests.analysis.test_seq_utils import TestSeqUtilsSetup, TEST_BED


class TestBedIntersect(TestSeqUtilsSetup, unittest.TestCase):
    """Tests for core_utils.feature_file_utils.intersect_features."""

    def test_intersect_features(self):
        """Test intersection of BED, GFF/GTF files."""

        with tempfile.NamedTemporaryFile(suffix=".test.intersect.bed", delete=False) as intersect_bed:
            output_fn = intersect_bed.name
            ffu.intersect_features(self.test_bed_a, self.test_bed_b, output_fn)

        # Note: I tried flushing the file in context and then re-reading but that does not seem to work;
        # Ensure the file's __exit__ method is ran for proper closing
        with open(output_fn, "r") as output_fh:
            observed = output_fh.read()
            self.assertEqual(observed, self.test_bed_b_str)

    def test_intersect_features_kwarg(self):
        """Test intersection of BED, GFF/GTF files with curry of kwarg."""

        with tempfile.NamedTemporaryFile(suffix=".test.intersect.bed", delete=False) as intersect_bed:
            output_fn = intersect_bed.name
            ffu.intersect_features(self.test_bed_a, self.test_bed_b, output_fn, v=True)

        with open(output_fn, "r") as output_fh:
            observed = output_fh.read()
            expected = "\n".join(TEST_BED.splitlines()[1:]) + "\n"
            self.assertEqual(observed, expected)

    def test_intersect_features_bedtool(self):
        """Test intersection of BED, GFF/GTF files."""

        intersect_bedtool = ffu.intersect_features(self.test_bed_a, self.test_bed_b, as_bedtool=True)
        self.assertTrue(pybedtools.BedTool == type(intersect_bedtool))


class TestFeatureFileUtils(TestSeqUtilsSetup, unittest.TestCase):
    """Tests for core_utils.feature_file_utils general functions."""

    def test_sort_feature_file(self):
        """Test that a feature file is lexicographically sorted."""

        unsorted_lines = """chr19	59066354	59066491	gi|571026644|ref|NM_014453.3|:1-137	0	-
chr19	59065411	59065603	gi|571026644|ref|NM_014453.3|:138-329	0	-
"""
        with tempfile.NamedTemporaryFile("w", suffix=".test.unsorted.bed", delete=False) as unsorted:
            unsorted.write(unsorted_lines)
            unsorted_fn = unsorted.name

        with tempfile.NamedTemporaryFile("w+", suffix=".test.sorted.bed", delete=False) as sorted:
            sorted_fn = sorted.name
            ffu.sort_feature_file(unsorted_fn, sorted_fn)

        with open(sorted_fn, "r") as output_fh:
            observed = output_fh.read()

            expected = """chr19	59065411	59065603	gi|571026644|ref|NM_014453.3|:138-329	0	-
chr19	59066354	59066491	gi|571026644|ref|NM_014453.3|:1-137	0	-
"""
            self.assertEqual(observed, expected)

    def test_store_coords(self):
        """Test that we can store coordinates from a feature file."""

        # COORD_TUPLE = collections.namedtuple("coord_tuple", "contig, start, stop, name, strand, score, allowable_coords")

        observed = ffu.store_coords(self.test_bed_b)
        expected = collections.OrderedDict({"gi|571026644|ref|NM_014453.3|:1-137": ffu.COORD_TUPLE(
            "chr19", 59066354, 59066491, "gi|571026644|ref|NM_014453.3|:1-137", Strand("-"), 0.0, frozenset(range(59066488, 59066494 + 1)))
                                           }.items())

        self.assertEqual(observed, expected)

    def test_store_coords_primer_allowable(self):
        """Test that we can store coordinates from a feature file, including the feature length in the allowable starts."""

        observed = ffu.store_coords(self.test_bed_b, primer_allowable=True)
        expected = collections.OrderedDict({"gi|571026644|ref|NM_014453.3|:1-137": ffu.COORD_TUPLE(
            "chr19", 59066354, 59066491, "gi|571026644|ref|NM_014453.3|:1-137", Strand("-"), 0.0, frozenset(range(59066354, 59066494 + 1)))
                                           }.items())

        self.assertEqual(observed, expected)

    def test_store_coords_use_name(self):
        """Test that we can store coordinates from a feature file, keying by coordinate instead of name."""

        observed = ffu.store_coords(self.test_bed_b, use_name=False)
        expected = collections.OrderedDict({COORD_FORMAT_STRAND.format("chr19", 59066354, 59066491, "-"): ffu.COORD_TUPLE(
            "chr19", 59066354, 59066491, "gi|571026644|ref|NM_014453.3|:1-137", Strand("-"), 0.0, frozenset(range(59066488, 59066494 + 1)))
                                           }.items())

        self.assertEqual(observed, expected)