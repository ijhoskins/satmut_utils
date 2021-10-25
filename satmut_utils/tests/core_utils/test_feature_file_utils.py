#!/usr/bin/env/python
""" Tests for core_utils.feature_file_utils."""

import collections
import pybedtools
import tempfile
import unittest

from analysis.seq_utils import Strand, COORD_FORMAT_STRAND
import core_utils.feature_file_utils as ffu
import core_utils.file_utils as fu
from tests.analysis.test_seq_utils import TEST_FASTA, TEST_BED


class TestBedIntersect(unittest.TestCase):
    """Tests for core_utils.feature_file_utils.intersect_features."""

    APPRIS_REF = "GRCh38.chr21.fa.gz"

    @classmethod
    def setUpClass(cls):
        """Setup for TestBedIntersect."""

        cls.tempdir = tempfile.mkdtemp()
        cls.seq = "ATCGATTACG"

        # We'll need another BED file to check intersection, etc.
        cls.test_bed_b_str = "".join(TEST_BED.split("\n")[0]) + fu.FILE_NEWLINE

        with tempfile.NamedTemporaryFile("w", suffix=".test.fa", delete=False, dir=cls.tempdir) as test_fasta, \
                tempfile.NamedTemporaryFile("w", suffix=".test.a.bed", delete=False, dir=cls.tempdir) as test_bed_a, \
                tempfile.NamedTemporaryFile("w", suffix=".test.b.bed", delete=False, dir=cls.tempdir) as test_bed_b:
            test_fasta.write(TEST_FASTA)
            test_bed_a.write(TEST_BED)
            test_bed_b.write(cls.test_bed_b_str)

            cls.test_fasta = test_fasta.name
            cls.test_bed_a = test_bed_a.name
            cls.test_bed_b = test_bed_b.name

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestBedIntersect."""

        fu.safe_remove((cls.tempdir,), force_remove=True)

    def test_intersect_features(self):
        """Test intersection of BED, GFF/GTF files."""

        with tempfile.NamedTemporaryFile(suffix=".test.intersect.bed", delete=False, dir=self.tempdir) as intersect_bed:
            output_fn = intersect_bed.name
            ffu.intersect_features(self.test_bed_a, self.test_bed_b, output_fn)

        with open(output_fn, "r") as output_fh:
            observed = output_fh.read()
            self.assertEqual(observed, self.test_bed_b_str)

    def test_intersect_features_kwarg(self):
        """Test intersection of BED, GFF/GTF files with curry of kwarg."""

        with tempfile.NamedTemporaryFile(suffix=".test.intersect.bed", delete=False, dir=self.tempdir) as intersect_bed:
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


class TestFeatureFileUtils(unittest.TestCase):
    """Tests for core_utils.feature_file_utils functions."""

    @classmethod
    def setUpClass(cls):
        """Setup for TestFeatureFileUtils."""

        cls.tempdir = tempfile.mkdtemp()
        cls.seq = "ATCGATTACG"

        # We'll need another BED file to check intersection, etc.
        cls.test_bed_b_str = "".join(TEST_BED.split("\n")[0]) + fu.FILE_NEWLINE

        with tempfile.NamedTemporaryFile("w", suffix=".test.b.bed", delete=False, dir=cls.tempdir) as test_bed_b:
            test_bed_b.write(cls.test_bed_b_str)
            cls.test_bed_b = test_bed_b.name

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestFeatureFileUtils."""

        fu.safe_remove((cls.tempdir,), force_remove=True)

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

        observed = ffu.store_coords(self.test_bed_b)
        expected = collections.OrderedDict({"gi|571026644|ref|NM_014453.3|:1-137": ffu.COORD_TUPLE(
            "chr19", 59066354, 59066491, "gi|571026644|ref|NM_014453.3|:1-137",
            Strand("-"), 0.0, frozenset(range(59066488, 59066494 + 1)))
                                           }.items())

        self.assertEqual(observed, expected)

    def test_store_coords_primer_allowable(self):
        """Test that we can store coordinates from a feature file, including the feature length in the allowable starts."""

        observed = ffu.store_coords(self.test_bed_b, primer_allowable=True)
        expected = collections.OrderedDict({"gi|571026644|ref|NM_014453.3|:1-137": ffu.COORD_TUPLE(
            "chr19", 59066354, 59066491, "gi|571026644|ref|NM_014453.3|:1-137",
            Strand("-"), 0.0, frozenset(range(59066354, 59066494 + 1)))
                                           }.items())

        self.assertEqual(observed, expected)

    def test_store_coords_use_name(self):
        """Test that we can store coordinates from a feature file, keying by coordinate instead of name."""

        observed = ffu.store_coords(self.test_bed_b, use_name=False)
        expected = collections.OrderedDict({COORD_FORMAT_STRAND.format("chr19", 59066354, 59066491, "-"): ffu.COORD_TUPLE(
            "chr19", 59066354, 59066491, "gi|571026644|ref|NM_014453.3|:1-137",
            Strand("-"), 0.0, frozenset(range(59066488, 59066494 + 1)))
                                           }.items())

        self.assertEqual(observed, expected)