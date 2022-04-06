#!/usr/bin/env python3
"""Tests for core_utils.feature_file_utils."""

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
        with tempfile.NamedTemporaryFile("w", suffix=".test.unsorted.bed", delete=False, dir=self.tempdir) as unsorted_fh:
            unsorted_fh.write(unsorted_lines)
            unsorted_fn = unsorted_fh.name

        with tempfile.NamedTemporaryFile("w", suffix=".test.sorted.bed", delete=False, dir=self.tempdir) as sorted_fh:
            sorted_fn = sorted_fh.name
            ffu.sort_feature_file(unsorted_fn, sorted_fn)

        with open(sorted_fn, "r") as output_fh:
            observed = output_fh.read()

            expected = """chr19	59065411	59065603	gi|571026644|ref|NM_014453.3|:138-329	0	-
chr19	59066354	59066491	gi|571026644|ref|NM_014453.3|:1-137	0	-
"""
            self.assertEqual(observed, expected)

    def test_store_coords(self):
        """Test that we can store coordinates from a feature file."""

        observed = ffu.store_coords(self.test_bed_b, use_name=True)
        expected = collections.OrderedDict({"gi|571026644|ref|NM_014453.3|:1-137": ffu.COORD_TUPLE(
            "chr19", 59066354, 59066491, "gi|571026644|ref|NM_014453.3|:1-137",
            Strand("-"), 0.0, frozenset(range(59066355, 59066492)))}.items())

        self.assertEqual(observed, expected)

    def test_store_coords_use_name(self):
        """Test that we can store coordinates from a feature file, keying by coordinate instead of name."""

        observed = ffu.store_coords(self.test_bed_b, use_name=False)
        expected = collections.OrderedDict({COORD_FORMAT_STRAND.format("chr19", 59066354, 59066491, "-"):
            ffu.COORD_TUPLE(
                "chr19", 59066354, 59066491, "gi|571026644|ref|NM_014453.3|:1-137",
                Strand("-"), 0.0, frozenset(range(59066355, 59066492)))}.items())

        self.assertEqual(observed, expected)


class TestSlopFeatures(unittest.TestCase):
    """Tests for core_utils.feature_file_utils.slop_features."""

    @classmethod
    def setUpClass(cls):
        """Setup for TestFeatureFileUtils."""

        cls.tempdir = tempfile.mkdtemp()
        cls.test_bed_b_str = "".join(TEST_BED.split("\n")[0]) + fu.FILE_NEWLINE

        with tempfile.NamedTemporaryFile("w", suffix=".test.b.bed", delete=False, dir=cls.tempdir) as test_bed_b:
            test_bed_b.write(cls.test_bed_b_str)
            cls.test_bed_b = test_bed_b.name

        with tempfile.NamedTemporaryFile("w", suffix=".genome_file.txt", delete=False, dir=cls.tempdir) as genome_file:
            genome_file_str = fu.FILE_DELIM.join(("chr19", str(59128983),)) + fu.FILE_NEWLINE
            genome_file.write(genome_file_str)
            cls.genome_file = genome_file.name

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestFeatureFileUtils."""

        fu.safe_remove((cls.tempdir,), force_remove=True)

    def test_slop_features(self):
        """Test that features can be slopped in both directions."""

        with tempfile.NamedTemporaryFile(suffix=".test.slop.bed", delete=False, dir=self.tempdir) as slop_bed:
            output_fn = slop_bed.name
            ffu.slop_features(self.test_bed_b, genome_file=self.genome_file, bp_left=10, bp_right=5,
                              by_strand=True, output=output_fn)

        with open(output_fn, "r") as output_fh:
            observed = output_fh.read()

            slopped_test_bed_b_fields = self.test_bed_b_str.split(fu.FILE_DELIM)
            slopped_test_bed_b_fields[1] = str(int(slopped_test_bed_b_fields[1]) - 5)
            slopped_test_bed_b_fields[2] = str(int(slopped_test_bed_b_fields[2]) + 10)
            expected = fu.FILE_DELIM.join(slopped_test_bed_b_fields)

            self.assertEqual(observed, expected)

    def test_slop_features_strand_unaware(self):
        """Test that features can be slopped in both directions."""

        with tempfile.NamedTemporaryFile(suffix=".test.slop.bed", delete=False, dir=self.tempdir) as slop_bed:
            output_fn = slop_bed.name
            ffu.slop_features(self.test_bed_b, genome_file=self.genome_file, bp_left=10, bp_right=5,
                              by_strand=False, output=output_fn)

        with open(output_fn, "r") as output_fh:
            observed = output_fh.read()

            slopped_test_bed_b_fields = self.test_bed_b_str.split(fu.FILE_DELIM)
            slopped_test_bed_b_fields[1] = str(int(slopped_test_bed_b_fields[1]) - 10)
            slopped_test_bed_b_fields[2] = str(int(slopped_test_bed_b_fields[2]) + 5)
            expected = fu.FILE_DELIM.join(slopped_test_bed_b_fields)

            self.assertEqual(observed, expected)
