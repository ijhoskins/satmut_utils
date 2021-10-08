#!/usr/bin/env/python
"""Tests for analysis.aligners."""

import pysam
import tempfile
import unittest

import analysis.aligners as al
import analysis.seq_utils as su
import core_utils.file_utils as fu
from definitions import *

tempfile.tempdir = os.getenv("SCRATCH", "/tmp")


class TestAligners(unittest.TestCase):
    """Tests for read aligners."""

    CBS_PEZY3_REF = "CBS_pEZY3.fa"

    @classmethod
    def setUpClass(cls):
        """Set up for TestAligners."""

        cls.tempdir = tempfile.mkdtemp()
        cls.test_dir = os.path.dirname(__file__)
        cls.test_data_dir = os.path.abspath(os.path.join(cls.test_dir, "..", "test_data"))
        cls.ref = os.path.join(cls.test_data_dir, cls.CBS_PEZY3_REF)

        cls.wt_seq = "GGAGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAATCCTGGGAATGGTGACGCTTGGGAACATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA"
        cls.snp_seq = su.introduce_error(cls.wt_seq, 0.05)
        cls.indel_seq = su.introduce_indels(cls.wt_seq, 0.05)
        fasta_lines = [">WT", cls.wt_seq, ">SNPs", cls.snp_seq, ">InDels", cls.indel_seq]
        fasta_string = fu.FILE_NEWLINE.join(fasta_lines) + fu.FILE_NEWLINE

        cls.test_fasta = tempfile.NamedTemporaryFile(suffix=".test.fa", delete=False).name
        with open(cls.test_fasta, "w") as test_fa_fh:
            test_fa_fh.write(fasta_string)

        cls.test_bam = tempfile.NamedTemporaryFile(suffix=".test.bam", delete=False).name

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestAligners."""

        fu.safe_remove((cls.test_fasta, cls.test_bam,))


class TestBowtie2(TestAligners):
    """Tests for Bowtie2."""

    def test_args(self):
        """Test proper passing of args."""

        # Pass the -f flag to align a FASTA
        bowtie_config = al.BowtieConfig(self.ref, True, 1, "f")
        _ = al.Bowtie2(config=bowtie_config, f1=self.test_fasta, output_bam=self.test_bam)

        # To count reads we have to index
        su.index_bam(self.test_bam)

        with pysam.AlignmentFile("rb", self.test_bam) as test_af:
            mapped_count = test_af.mapped

        self.assertEqual(mapped_count, 3)

    def test_kwargs(self):
        """Test proper passing of kwargs."""

        # Pass the --skip flag to align the last two reads
        bowtie_config = al.BowtieConfig(self.ref, local=True, nthreads=1, skip=1)
        _ = al.Bowtie2(config=bowtie_config, f1=self.test_fasta, output_bam=self.test_bam)

        observed_qnames = {align_seg.query_name for align_seg in pysam.AlignmentFile(self.test_bam, "rb").fetch()}

        self.assertTrue("WT" not in observed_qnames)

    def test_alignment(self):
        """Test proper alignment of various types of reads."""

        with pysam.AlignmentFile(self.test_bam, "rb") as af:
            for align_seg in af.fetch():
                if align_seg.query_name == "WT":
                    # This read should match the reference exactly
                    matches_ref = align_seg.cigarstring == "105M"
                if align_seg.query_name == "SNPs":
                    has_snps = align_seg.get_tag(su.SAM_EDIT_DIST_TAG) >= 0
                if align_seg.query_name == "InDels":
                    has_ins = "I" in align_seg.cigarstring
                    has_dels = "D" in align_seg.cigarstring

        self.assertTrue(matches_ref and has_snps and (has_ins or has_dels))
