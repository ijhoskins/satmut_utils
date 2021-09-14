#!/usr/bin/env/python
"""Tests for analysis.aligners."""

import pysam
import tempfile
import unittest

import analysis.aligners as al
import analysis.seq_utils as su
import core_utils.file_utils as fu
from definitions import *

tempfile.tempdir = "/tmp"


class TestAligners(unittest.TestCase):
    """Tests for read aligners."""

    @classmethod
    def setUpClass(cls):
        """Set up for TestAligners."""

        cls.wt_seq = "GACGGACACCCGTTTCCGCTTCCGGGTCACGCTTCTCTTTCTGGGATCCCCGACTTGCCCACCAACTAAGGCCCCTCGGTCCTGTCGCCGCCGCCGCCGTTTCCGGATTAAACGACGTGACGTAACATGCCCCGCCCGCACCCGGAACGT"
        cls.snp_seq = su.introduce_error(cls.wt_seq, 0.05)
        cls.indel_seq = su.introduce_indels(cls.wt_seq, 0.05)
        fasta_lines = [">WT", cls.wt_seq, ">SNPs", cls.snp_seq, ">InDels", cls.indel_seq]
        fasta_string = fu.FILE_NEWLINE.join(fasta_lines) + fu.FILE_NEWLINE

        cls.test_fasta = tempfile.NamedTemporaryFile("w+", suffix=".test.fa", delete=False).name
        with open(cls.test_fasta) as test_fa_fh:
            test_fa_fh.write(fasta_string)

        cls.test_bam = tempfile.NamedTemporaryFile("w+b", suffix=".test.bam", delete=False).name

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestAligners."""

        fu.safe_remove((cls.test_fasta, cls.test_bam,))


class TestBowtie2(TestAligners):
    """Tests for Bowtie2."""

    def test_args(self):
        """Test proper passing of args."""

        # Pass the -f flag to align a FASTA
        bowtie_config = al.BowtieConfig(GRCH38_FASTA, True, 1, "f")
        _ = al.Bowtie2(config=bowtie_config, f1=self.test_fasta, output_bam=self.test_bam)

        # To count reads we have to index
        su.index_bam(self.test_bam)

        with pysam.AlignmentFile("rb", self.test_bam) as test_af:
            mapped_count = test_af.mapped

        self.assertEqual(mapped_count, 3)

    def test_kwargs(self):
        """Test proper passing of kwargs."""

        # Pass the --skip flag to align the last two reads
        bowtie_config = al.BowtieConfig(GRCH38_FASTA, local=True, nthreads=1, skip=1)
        _ = al.Bowtie2(config=bowtie_config, f1=self.test_fasta, output_bam=self.test_bam)

        valid_qnames = {align_seg.query_name for align_seg in
                        pysam.AlignmentFile("rb", self.test_bam).fetch(until_eof=True)}

        self.assertTrue("WT" not in valid_qnames)

    def test_alignment(self):
        """Test proper alignment of various types of reads."""

        has_soft_clip = False
        has_snps = False
        has_ins = False
        has_dels = False

        with pysam.AlignmentFile(self.test_bam.name, "rb") as af:
            for align_seg in af.fetch(until_eof=True):

                if align_seg.is_unmapped:
                    continue

                # These could potentially change with different alignment params
                if align_seg.query_name == "WT":
                    # A terminal error in this read causes soft clipping
                    has_soft_clip = "S" in align_seg.cigarstring
                if align_seg.query_name == "SNPs":
                    has_snps = align_seg.get_tag(su.SAM_EDIT_DIST_TAG) >= 0
                if align_seg.query_name == "InDels":
                    has_ins = "I" in align_seg.cigarstring
                    has_dels = "D" in align_seg.cigarstring

        self.assertTrue(has_soft_clip and has_snps and has_ins and has_dels)
