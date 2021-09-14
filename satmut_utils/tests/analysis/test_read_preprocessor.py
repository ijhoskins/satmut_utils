#!/usr/bin/env/python
"""Tests for analysis.read_preprocessor"""

import os
import pysam
import re
import tempfile
import unittest

import analysis.read_preprocessor as rp
import analysis.seq_utils as su
import core_utils.file_utils as fu
from definitions import *

tempfile.tempdir = "/tmp"

TEST_DATA_DIR = "/Users/ianhoskins/bioinfo/tests/test_data"

# First pair contains no adapter readthrough and is first pair in CBS1_65

# Second pair contains typical adapter readthrough and was selected from CBS1_63 by grep first instance of adapter match

# Third pair contains an unexpected presence of the 5' adapter in R1, which appears to be caused by PCR1 target primer
# carry-through into PCR2 and mispriming off the adapter; was selected from CBS1_63 by match first instance of
# 5' adapter in R1; the R2 for this contains two instances of the RC of this 5' P7 adapter

TEST_R1_TILESEQ_FASTQ = """@A00738:253:HF5HHDSX2:2:1101:27516:1078 1:N:0:GAATTCGT+ACGTCCTG
GAGGAGGCGTTCACCTTTGCCCGCATGCTGATCGCGCAAGAGGGGCTGCTGTGCGGTGGCAGTGCTGGCAGCACGGTGGCGGTGGCCGTGAAGGCTGCGCAGGAGCTGCAGGAGGGCCAGCGCTGCGTGGTCATTCTGCCCGACTCAGTG
+
FFFFFFFFFFFFFFF:FF:FFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFFFFFF,FFFFFFFFFFF:FFFFFFFFFFF:FFFFFFFF::FFFFFFFFFFFF:FFFFFFFFFFFFFFFFF:FFFFFFFF
@A00738:253:HF5HHDSX2:2:1101:27290:1063 1:N:0:CGCTCATT+CTTCGCCT
CCAATTCTCACATCCTAGACCAACACGACGCTCTTCCGATCTCCAATTCTCACATCCTAGACCATCACGGGCATTGCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGCTCATTATCGCGGATGCCGTCGTGGGGGGGGGGGGGG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FF,:F,:F:FF:F:F::F,,F:F::FFFF,,F:FFFF:,:,FF:FFFFFFFFFF:FF:F:F,F,FFFFFFFFFFFFF:FF,,,,:::,F,,,,,:F,,,,,,F,,:F,,,,::F
@A00738:253:HF5HHDSX2:2:1101:14281:2237 1:N:0:CGCTCATT+CTTCGCCT
CCAATTCTCACATCCTAGACCATACACGACGCTCTTCCGATCTCCAATTCTCACATCCTATACCATCACGGGCATTGCAGATCGGAAGAGCACACGGCTGAACTCCAGTCACCGCTCATTATCGGGGGTGCCGGGGGGGGGGGGGGGGGG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFF,:::F,:::FFF:,,F,F:,,FF,,,,:F::FFFF:F,,,:F,F:::,,F:,,FFF::FFF,,FFF,FFF::,FFF:,:F:,,,,:,,,:,F,,:F,,F,,,,,:F,F:F:,,,F:
"""

TEST_R2_TILESEQ_FASTQ = """@A00738:253:HF5HHDSX2:2:1101:27516:1078 2:N:0:GAATTCGT+ACGTCCTG
GCACTGAGTCGGGCAGAATGACCACGCAGCGCTGGCCCTCCTGCAGCTCCTGCGCAGCCTTCACGGCCACCGCCACCGTGCTGCCAGCACTGCCACCGCACAGCAGCCCCTCTTGCGCGATCAGCATGCGGGCAAAGGTGAACGCCTCCT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00738:253:HF5HHDSX2:2:1101:27290:1063 2:N:0:CGCTCATT+CTTCGCCT
GCAATGCCCGTGATGGTCTAGGATGTGAGAATTGGAGATCGGAAGAGCGTCGTGTTGGTCTAGGATGTGAGAATTGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCTTCGCCTGTGTAGATCTCGGTGGTCGCCGTATCATTAAA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF,FFFFF,FFFFFFFFFFFFFFFFFFF:F:FFF:FFFFFF
@A00738:253:HF5HHDSX2:2:1101:14281:2237 2:N:0:CGCTCATT+CTTCGCCT
GCAATGCCCGTGATGGTCTAGGATGTGAGAATTGGAGATCGGAAGAGCGTCGTGTATGGTCTAGGATGTGAGAATTGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCTTCGCCTGTGTAGATCTCGGTGGTCGCCGTATCATTAA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FF:F::FFFF:FFFFFFFFFFFF:FFF:FFFFF:FFFFF:FF,F:F,FFFFF:FF:FFFFFFFF,F:FFFFF:F,FFF
"""

# These contain an 8 bp UMI with downstream attB sequence and then target to test trim of PEZY3_ATTB1_P7 and its RC
# first pair of CBS1_67
# UMI in this dataset was 5' of flanking vector-target sequence TGCTGGAATTGATTAATA
TEST_R1_UMI_TILESEQ_FASTQ = """@A00738:253:HF5HHDSX2:2:1101:13964:1063 1:N:0:TAATGCGC+AGGCTATA
GCTCGGCGTGCTGGAATTGATTAATACAAGTTTGTACAAAAAAGTTGGCATGCCTTCTGAGACCCCCCAGGCAGAAGTGGGGCCCACAGGCTGCCCCCACCGCTCAGGGCCACACTCGGCGAAGGGGAGCCTGGAGAAGGGGTCCCCAGA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
"""


TEST_R2_UMI_TILESEQ_FASTQ = """@A00738:253:HF5HHDSX2:2:1101:13964:1063 2:N:0:TAATGCGC+AGGCTATA
CACAGGGGCTCCTTGGCTTCCTTATCCTCTGGGGACCCCTTCTCCAGGCTCCCCTTCGCCGAGTGTGGCCCTGAGCGGTGGGGGCAGCCTGTGGGCCCCACTTCTGCCTGGGGGGTCTCAGAAGGCATGCCAACTTTTTTGTACAAACTT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFFFFFFFF:FFFF
"""

# Two Gs in the primer were changed to C
TEST_R2_UMI_TILESEQ_PRIMER_ERROR_FASTQ = """@A00738:253:HF5HHDSX2:2:1101:13964:1063 2:N:0:TAATGCGC+AGGCTATA
CACAGCCGCTCCTTGGCTTCCTTATCCTCTGGGGACCCCTTCTCCAGGCTCCCCTTCGCCGAGTGTGGCCCTGAGCGGTGGGGGCAGCCTGTGGGCCCCACTTCTGCCTGGGGGGTCTCAGAAGGCATGCCAACTTTTTTGTACAAACTT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFFFFFFFF:FFFF
"""

# The primer was removed
TEST_R2_UMI_TILESEQ_NOPRIMER_FASTQ = """@A00738:253:HF5HHDSX2:2:1101:13964:1063 2:N:0:TAATGCGC+AGGCTATA
TCCTTATCCTCTGGGGACCCCTTCTCCAGGCTCCCCTTCGCCGAGTGTGGCCCTGAGCGGTGGGGGCAGCCTGTGGGCCCCACTTCTGCCTGGGGGGTCTCAGAAGGCATGCCAACTTTTTTGTACAAACTT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFFFFFFFF:FFFF
"""

TILESEQ_CBS_1R_FASTA = """>CBS_pEZY3:1083-1101(-)
CACAGGGGCTCCTTGGCT
"""


# First pair in CBS1_35 was selected, which has a very short insert; in R1 common region observed in addition to GSP2 RC
# in R2, RC of common region observed
TEST_R1_UMI_AMP_FASTQ = """@MG01HS02:1483:HG7MTBCX3:1:1107:1230:2106 1:N:0:TAAGGCGA+NAGGCTTA
TATGGGCGAACCGCCAGGAGTCACTTGGCCAAGAGCTCACACTTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTGAAAAAANAAAAAGTAACAGTCAAATAGAAAGGTTAAATTT
+
DCDDAHCHIIIICFHI<C@GHHFFGHHIIIIHHHHIFHHIH?EGHIHIHHIHIH?C@HEFHIHHHIHGCDFHHHE@GHIHHHH?HCEH/<DG?GCEG1CE<<<<<FD1DH11DC1CE#<<11101<11<1111001<000<00<00<00<<
"""

TEST_R2_UMI_AMP_FASTQ = """@MG01HS02:1483:HG7MTBCX3:1:1107:1230:2106 2:N:0:TAAGGCGA+NAGGCTTA
NNNNNGTGAGCTCTTGGCCAAGTGACTCCTGGCGGTTCGCCCATAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTTAAGCCTCGTTTAGATGATCGTGGCCTCAGCGTCAAATAAAAAAAAAAAAACAAAACCGAAGAAAAAACCAAA
+
#####<<<C@?C@FG1<CGEHIIIIIIIHHIIHH?EE1FHHHHHHE@GHCHEEED?@HHHID/DDF?HIGHHIIHEC<F@F<FF1<C11111111111<1<</<111<1/<<D0<1<1CGHII//<</<CG0000=/<<<DDH/F/=//=/
"""

# Contains duplicate reads with no error, mismatches, and deletions
DUP_TEST_BAM = """@HD	VN:1.0	SO:coordinate
@SQ	SN:CBS_pEZY3	LN:7108
@RG	ID:CBS1_35_comb_R
@PG	ID:bowtie2	PN:bowtie2	VN:2.4.2	CL:"/home/ihoskins/miniconda3/envs/CBS_variants/bin/bowtie2-align-s --wrapper basic-0 -p 5 --maxins 1000 --no-discordant --fr --rf --mp 4 --rdg 6,4 --rfg 6,4 --local --rg-id CBS1_35_comb_R -x /home/ihoskins/reference_files/CBS_pEZY3.fa -1 /scratch/users/ihoskins/CBS_02APR2020/pEZY3_default_dedup/CBS1_35_comb_R1.umi.trimmed.fq -2 /scratch/users/ihoskins/CBS_02APR2020/pEZY3_default_dedup/CBS1_35_comb_R2.umi.trimmed.fq"
@PG	ID:samtools	PN:samtools	PP:bowtie2	VN:1.10 (pysam)	CL:samtools sort -o /scratch/users/ihoskins/tmpq77y8hyj.bam -O BAM -@ 5 -n /scratch/users/ihoskins/CBS_02APR2020/pEZY3_default_dedup/CBS1_35_comb_R.group.bam
@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.10 (pysam)	CL:samtools view -b -f 66 -o /scratch/users/ihoskins/tmp23xggiwu.bam -O BAM -@ 5 /scratch/users/ihoskins/tmpq77y8hyj.bam
@PG	ID:samtools.2	PN:samtools	PP:samtools.1	VN:1.10 (pysam)	CL:samtools sort -o /scratch/users/ihoskins/CBS_02APR2020/pEZY3_default_dedup/CBS1_35_comb_R.preprocess.bam -O BAM -@ 5 -t UG /scratch/users/ihoskins/tmpbwcoeu1f.tag.bam
MG01HS02:1483:HG7MTBCX3:1:2108:10337:91772_GCATATCG	83	CBS_pEZY3	1062	44	6M1D63M	=	1062	-70	AGGGGTCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCAC	EHGC<<EHFD1GEEHHIHH@EIIHEHHHHE?E=IHHHEG<0/HDHHHHHIHHIIHHHIIIHIHFHIHHI	AS:i:128	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:6^C63	YS:i:128	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:GCATATCG	UG:Z:1000009_R1
MG01HS02:1483:HG7MTBCX3:1:2108:10337:91772_GCATATCG	163	CBS_pEZY3	1062	44	6M1D63M	=	1062	70	AGGGGTCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCAC	ABAADEGHHHHI?C1FHIE@FHIIIIEH<GEHIIIC<FHHEHIIIIIIHIIHIIHHEEHCCHIIECGHI	AS:i:128	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:6^C63	YS:i:128	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:1000009_R2
MG01HS02:1483:HG7MTBCX3:1:2108:10337:91772_GCATATCG	83	CBS_pEZY3	1062	44	6M1D63M	=	1062	-70	AGGGGTCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCAC   EHGC<<EHFD1GEEHHIHH@EIIHEHHHHE?E=IHHHEG<0/HDHHHHHIHHIIHHHIIIHIHFHIHHI	AS:i:128	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:6^C63	YS:i:128	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:GCATATCG	UG:Z:1000009_R1
MG01HS02:1483:HG7MTBCX3:1:2213:19735:62464_GCGTATCG	83	CBS_pEZY3	1062	44	69M	=	1062	-70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   HE=DHDCFC<C11D<1ECHHC<<1CG@FEC/C/<1HEHHFDE<HHHHGHGHEHHCCDDIHGCCGGCGFCE	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:GCATATCG   UG:Z:1000009_R1
MG01HS02:1483:HG7MTBCX3:2:2210:18370:5712_GCGTATCG	83	CBS_pEZY3	1062	44	69M	=	1062	-70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   HIIIIIIHHHEHHIHIIHHIIHHHIHG@IHGHHDHF?GEEDHCCCC</GHHCEHHIIHDCIIHHHHIHGI	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:GCATATCG   UG:Z:1000009_R1
MG01HS02:1506:HGKGCBCX3:2:1102:18758:39821_GCATATCG	83	CBS_pEZY3	1062	44	69M	=	1062	-70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   IIIIIIIHEIHIIIIIIHIIIIIIIIHGIHIIIGIIIIIIIIIIDIIHIIIIIIIIIIIIIIIIIIIIII	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:GCATATCG   UG:Z:1000009_R1
MG01HS02:1506:HGKGCBCX3:2:1105:17502:94042_GCATATCG	83	CBS_pEZY3	1062	44	69M	=	1062	-70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   GCIIIIIIHHFHIHGIIIIIIIHHHCH?IHIIHHEIIIIIIIHIIHHIHIIHHIHHIHIHIIIIIHIHIH	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:GCATATCG   UG:Z:1000009_R1
MG01HS02:1483:HG7MTBCX3:1:2108:10337:91772_GCATATCG	163	CBS_pEZY3	1062	44	6M1D63M	=	1062	70	AGGGGTCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA    CABAADEGHHHHI?C1FHIE@FHIIIIEH<GEHIIIC<FHHEHIIIIIIHIIHIIHHEEHCCHIIECGHI	AS:i:128	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:6^C63	YS:i:128	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:1000009_R2
MG01HS02:1483:HG7MTBCX3:1:2213:19735:62464_GCGTATCG	163	CBS_pEZY3	1062	44	69M	=	1062	70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   D@@00<0FGEE@F1<D@H1FHHE?1DGCC@1F1CGHHGHEHHH?1FHH?EHCHHHCEHHHDHHICD1F?G	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:1000009_R2
MG01HS02:1483:HG7MTBCX3:2:2210:18370:5712_GCGTATCG	163	CBS_pEZY3	1062	44	69M	=	1062	70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   D<D@DEHEHHHGHHIIHHIIIIHFHHHICHEHHIHHEHHHHFFEHFHIIIEHHIIIIIHIHIHIHIHEHH	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:1000009_R2
MG01HS02:1506:HGKGCBCX3:2:1102:18758:39821_GCATATCG	163	CBS_pEZY3	1062	44	69M	=	1062	70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   DDDDDIIIIIHHIIIIHIIIIIICHIIIIIIIIIIIIIIIIIHIIIGIIIIIIIIIIIIIIIIIHCHHHI	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:1000009_R2
MG01HS02:1506:HGKGCBCX3:2:1105:17502:94042_GCATATCG	163	CBS_pEZY3	1062	44	69M	=	1062	70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   DDDDDHGHIHIIHHGHIIGIIIIHHHIHIHHHIIIIIIIIHIGHHHIIHIIIIIIFHIIIIHGHHHIIII	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:1000009_R2
MG01HS02:1483:HG7MTBCX3:1:2108:6635:96297_AAATGGGG	83	CBS_pEZY3	1062	44	69M	=	1062	-70	AGAGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   <C@EF</@ECF<@C<<<<FD1F@HEHHHHHHCGIEHHF1F?EDE?EC/HCHHHHDEHHCEH@E@HEHEHH	AS:i:135	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:2G67	YS:i:135	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:AAATGGGG	UG:Z:1000004_R1
MG01HS02:1483:HG7MTBCX3:1:2108:12739:36922_AAATGGGG	83	CBS_pEZY3	1062	44	69M	=	1062	-70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   IIIGIIHHIIIIIIIIHHIHIIIIIIIIIHIIIIIIIIIIIIIIIIIIIIIIIIHHHIIHIIIGIIIIII	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:AAATGGGG   UG:Z:1000004_R1
MG01HS02:1506:HGKGCBCX3:1:2106:15083:48801_AAATGGGG	83	CBS_pEZY3	1062	44	69M	=	1062	-70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   IIIHHEHIGIHGHIIIIIIIIIIIIHIIIIIIIFIHIIIIIIHIIHHHEIIIIIIIIIIIIIIHIHIHHH	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:AAATGGGG   UG:Z:1000004_R1
MG01HS02:1483:HG7MTBCX3:1:2108:6635:96297_AAATGGGG	163	CBS_pEZY3	1062	44	69M	=	1062	70	AGAGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   BDD@@FE1D@HC1C<DF?GF?HGEHHCEF@@<EHE@CEH?GC1F@GHHIIGHDEGEH?HD</CHE<@CCH	AS:i:135	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:2G67	YS:i:135	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:1000004_R2
MG01HS02:1483:HG7MTBCX3:1:2108:12739:36922_AAATGGGG	163	CBS_pEZY3	1062	44	69M	=	1062	70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   DADDDIIIIIIIIIGHIHHIHHIGGIIIIHIGIIIIIIHHHHIHIIDCHIIHHIIIIIIHIGIHF?FEHI	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:1000004_R2
MG01HS02:1506:HGKGCBCX3:1:2106:15083:48801_AAATGGGG	163	CBS_pEZY3	1062	44	69M	=	1062	70	AGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTGGATCCGGCCCGATGCTCCGAGCAGGTGCA   DDADCHIIIIIHHHHHIIIHFHHHGHIGHGIIIIIIIFHHHEHHFEHIIHHHIHHIIIIIIHHIHCGHII	AS:i:140	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YS:i:140	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:1000004_R2
"""


class TestFastqPreprocessor(unittest.TestCase):
    """Tests for FastqPreprocessor."""

    @classmethod
    def setUpClass(cls):
        """Set up for TestFastqPreprocessor."""

        cls.tempdir = tempfile.mkdtemp()
        cls.tileseq_r1_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.R1.fastq", delete=False).name
        cls.tileseq_r2_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.R2.fastq", delete=False).name
        cls.tileseq_r1_umi_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.umi.R1.fastq", delete=False).name
        cls.tileseq_r2_umi_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.umi.R2.fastq", delete=False).name
        cls.amp_r1_fastq = tempfile.NamedTemporaryFile(suffix=".amp.umi.R1.fastq", delete=False).name
        cls.amp_r2_fastq = tempfile.NamedTemporaryFile(suffix=".amp.umi.R2.fastq", delete=False).name

        with open(cls.tileseq_r1_fastq, "w") as tileseq_r1_fh, \
                open(cls.tileseq_r2_fastq, "w") as tileseq_r2_fh, \
                open(cls.tileseq_r1_umi_fastq, "w") as tileseq_r1_umi_fh, \
                open(cls.tileseq_r2_umi_fastq, "w") as tileseq_r2_umi_fh, \
                open(cls.amp_r1_fastq, "w") as amp_r1_fh, \
                open(cls.amp_r2_fastq, "w") as amp_r2_fh:

            tileseq_r1_fh.write(TEST_R1_TILESEQ_FASTQ)
            tileseq_r2_fh.write(TEST_R2_TILESEQ_FASTQ)
            tileseq_r1_umi_fh.write(TEST_R1_UMI_TILESEQ_FASTQ)
            tileseq_r2_umi_fh.write(TEST_R2_UMI_TILESEQ_FASTQ)
            amp_r1_fh.write(TEST_R1_UMI_AMP_FASTQ)
            amp_r2_fh.write(TEST_R2_UMI_AMP_FASTQ)

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestFastqPreprocessor."""

        fu.safe_remove((cls.tileseq_r1_fastq, cls.tileseq_r2_fastq, cls.tileseq_r1_umi_fastq, cls.tileseq_r2_umi_fastq,
                        cls.amp_r1_fastq, cls.amp_r2_fastq))

    def test_run_cutadapt_tileseq_r1_no_adapter(self):
        """Tests lack of adapter trimming of Tile-seq R1 without adapter readthrough."""

        expected = "GAGGAGGCGTTCACCTTTGCCCGCATGCTGATCGCGCAAGAGGGGCTGCTGTGCGGTGGCAGTGCTGGCAGCACGGTGGCGGTGGCCGTGAAGGCTGCGCAGGAGCTGCAGGAGGGCCAGCGCTGCGTGGTCATTCTGCCCGACTCAGTG"

        # For CDS-terminal PCR tiles we also include PEZY3 flanking sequences
        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   outdir=self.tempdir, validate=False)

        with open(fqp.trimmed_f1, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):
                if i == 1:
                    self.assertEqual(expected, line)

    def test_run_cutadapt_tileseq_r1_normal_pair(self):
        """Tests adapter trimming of Tile-seq R1 with adapter readthrough."""

        expected = "CCAATTCTCACATCCTAGACCAACACGACGCTCTTCCGATCTCCAATTCTCACATCCTAGACCATCACGGGCATTGC"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   outdir=self.tempdir, validate=False)

        with open(fqp.trimmed_f1, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):
                if i == 5:
                    self.assertEqual(expected, line)

    def test_run_cutadapt_tileseq_r1_abnormal_pair(self):
        r"""Tests adapter trimming of Tile-seq R1 with target priming off the 5' adapter, leading to abnormal
        5' adapter seq presence."""

        expected = "CCAATTCTCACATCCTATACCATCACGGGCATTGC"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   outdir=self.tempdir, validate=False)

        with open(fqp.trimmed_f1, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):
                if i == 9:
                    self.assertEqual(expected, line)

    def test_run_cutadapt_tileseq_r2_no_adapter(self):
        """Tests lack of adapter trimming of Tile-seq R2 without adapter readthrough."""

        expected = "GCACTGAGTCGGGCAGAATGACCACGCAGCGCTGGCCCTCCTGCAGCTCCTGCGCAGCCTTCACGGCCACCGCCACCGTGCTGCCAGCACTGCCACCGCACAGCAGCCCCTCTTGCGCGATCAGCATGCGGGCAAAGGTGAACGCCTCCT"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   outdir=self.tempdir, validate=False)

        with open(fqp.trimmed_f2, "r") as r2_out_fh:
            for i, line in enumerate(r2_out_fh):
                if i == 1:
                    self.assertEqual(expected, line)

    def test_run_cutadapt_tileseq_r2_normal_pair(self):
        """Tests lack of adapter trimming of Tile-seq R1 without adapter readthrough."""

        expected = "GCAATGCCCGTGATGGTCTAGGATGTGAGAATTGG"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   outdir=self.tempdir, validate=False)

        with open(fqp.trimmed_f2, "r") as r2_out_fh:
            for i, line in enumerate(r2_out_fh):
                if i == 5:
                    self.assertEqual(expected, line)

    def test_run_cutadapt_tileseq_r2_abnormal_pair(self):
        r"""Tests adapter trimming of Tile-seq R2 with target priming off the R1 5' adapter, leading to abnormal
        R1 5' adapter seq presence."""

        expected = "GCAATGCCCGTGATGGTCTAGGATGTGAGAATTGG"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   outdir=self.tempdir, validate=False)

        with open(fqp.trimmed_f2, "r") as r2_out_fh:
            for i, line in enumerate(r2_out_fh):
                if i == 9:
                    self.assertEqual(expected, line)

    def test_run_cutadapt_tileseq_r1_umi(self):
        """Tests adapter trimming of Tile-seq UMI-containing R1 with adapter readthrough."""

        expected = "ATGCCTTCTGAGACCCCCCAGGCAGAAGTGGGGCCCACAGGCTGCCCCCACCGCTCAGGGCCACACTCGGCGAAGGGGAGCCTGGAGAAGGGGTCCCCAGA"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_umi_fastq, f2=self.tileseq_r2_umi_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   outdir=self.tempdir, validate=False)

        with open(fqp.trimmed_f1, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):
                if i == 1:
                    self.assertEqual(expected, line)

    def test_run_cutadapt_tileseq_r2_umi(self):
        """Tests adapter trimming of Tile-seq R2 with adapter readthrough paired with UMI-containing R1."""

        expected = "CACAGGGGCTCCTTGGCTTCCTTATCCTCTGGGGACCCCTTCTCCAGGCTCCCCTTCGCCGAGTGTGGCCCTGAGCGGTGGGGGCAGCCTGTGGGCCCCACTTCTGCCTGGGGGGTCTCAGAAGGCAT"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_umi_fastq, f2=self.tileseq_r2_umi_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   outdir=self.tempdir, validate=False)

        with open(fqp.trimmed_f2, "r") as r2_out_fh:
            for i, line in enumerate(r2_out_fh):
                if i == 1:
                    self.assertEqual(expected, line)

    def test_run_cutadapt_amp_r1_normal_pair(self):
        """Tests adapter trimming of AMP R1 containing 5' and 3' adapters."""

        # The anchored common region should be trimmed at the 5' end, and the GSP2 tail RC trimmed at the 3' end
        expected = "CACTTGGCCAAGAGCTCACACTTC"

        fqp = rp.FastqPreprocessor(f1=self.amp_r1_fastq, f2=self.amp_r2_fastq,
                                   r1_fiveprime_adapters=AMP_CR, r1_threeprime_adapters=AMP_GSP2_TAIL_RC,
                                   outdir=self.tempdir, validate=False)

        with open(fqp.trimmed_f1, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):
                if i == 1:
                    self.assertEqual(expected, line)

    def test_run_cutadapt_amp_r2_normal_pair(self):
        """Tests adapter trimming of AMP R2 containing 3' adapter."""

        # The anchored common region RC should be trimmed at the 3' end
        expected = "NNNNNGTGAGCTCTTGGCCAAGTG"

        fqp = rp.FastqPreprocessor(f1=self.amp_r1_fastq, f2=self.amp_r2_fastq,
                                   r1_fiveprime_adapters=AMP_CR, r1_threeprime_adapters=AMP_GSP2_TAIL_RC,
                                   outdir=self.tempdir, validate=False)

        with open(fqp.trimmed_f2, "r") as r2_out_fh:
            for i, line in enumerate(r2_out_fh):
                if i == 1:
                    self.assertEqual(expected, line)

    def test_validate(self):
        """Test that fastqc output files are generated."""

        _ = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   outdir=self.tempdir, validate=True)

        matches = [f for f in os.listdir(self.tempdir) if re.search("fastqc.html", f)]
        # We should have a file for R1 and R2
        self.assertEqual(len(matches), 2)


class TestUmiExtractor(unittest.TestCase):
    """Tests for UmiExtractor."""

    @classmethod
    def setUpClass(cls):
        """Set up for TestUmiExtractor."""

        cls.tempdir = tempfile.mkdtemp()
        cls.tileseq_r1_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.umi.R1.fastq", delete=False).name
        cls.tileseq_r2_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.umi.R2.fastq", delete=False).name
        cls.amp_r1_fastq = tempfile.NamedTemporaryFile(suffix=".amp.umi.R1.fastq", delete=False).name
        cls.amp_r2_fastq = tempfile.NamedTemporaryFile(suffix=".amp.umi.R2.fastq", delete=False).name
        cls.tileseq_umi_1r_fasta = tempfile.NamedTemporaryFile(suffix=".tileseq.umi.1R.fasta", delete=False).name

        with open(cls.tileseq_r1_fastq, "w") as tileseq_r1_fh, \
                open(cls.tileseq_r2_fastq, "w") as tileseq_r2_fh, \
                open(cls.amp_r1_fastq, "w") as amp_r1_fh, \
                open(cls.amp_r2_fastq, "w") as amp_r2_fh, \
                open(cls.tileseq_umi_1r_fasta, "w") as tileseq_umi_1r_fh:

            tileseq_r1_fh.write(TEST_R1_UMI_TILESEQ_FASTQ)
            tileseq_r2_fh.write(TEST_R2_UMI_TILESEQ_FASTQ)
            amp_r1_fh.write(TEST_R1_UMI_AMP_FASTQ)
            amp_r2_fh.write(TEST_R2_UMI_AMP_FASTQ)
            tileseq_umi_1r_fh.write(TEST_R2_UMI_TILESEQ_FASTQ)

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestUmiExtractor."""

        fu.safe_remove((cls.tileseq_r1_fastq, cls.tileseq_r2_fastq, cls.amp_r1_fastq, cls.amp_r2_fastq))

    def test_umitools_extract_tileseq(self):
        """Tests proper extraction of UMIs from Tile-seq data."""

        umi_extractor = rp.UMIExtractor(r1_fastq=self.tileseq_r1_fastq, r2_fastq=self.tileseq_r2_fastq,
                                        umi_regex=TILESEQ_UMI_REGEX, outdir=self.tempdir)

        # Here we expect two things:
        # 1) the 8 bp UMI GCTCGGCG has been moved to the read name using _ as the delimiter
        # 2) The UMI and adjacent anchor sequence are discarded from the sequence per TILESEQ_UMI_REGEX
        expected_1 = "@A00738:253:HF5HHDSX2:2:1101:13964:1063_GCTCGGCG 1:N:0:TAATGCGC+AGGCTATA"
        expected_2 = "CAAGTTTGTACAAAAAAGTTGGCATGCCTTCTGAGACCCCCCAGGCAGAAGTGGGGCCCACAGGCTGCCCCCACCGCTCAGGGCCACACTCGGCGAAGGGGAGCCTGGAGAAGGGGTCCCCAGA"

        with open(umi_extractor.r1_out_fastq, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):

                # Test that the UMI has moved to the qname
                if i == 0:
                    test_1 = line == expected_1

                # Test that the UMI sequence and adjacent adapter sequence were trimmed
                if i == 1:
                    test_2 = line ==  expected_2
                    break

        self.assertTrue(all((test_1, test_2)))

    def test_umitools_extract_amp(self):
        """Tests proper extraction of UMIs from RACE-like (e.g. AMP) data."""

        umi_extractor = rp.UMIExtractor(r1_fastq=self.amp_r1_fastq, r2_fastq=self.amp_r2_fastq,
                                        umi_regex=AMP_UMI_REGEX, outdir=self.tempdir)

        # Here we expect two things:
        # 1) the 8 bp UMI TATGGGCG has been moved to the read name using _ as the delimiter
        # 2) The UMI and adjacent anchor sequence are discarded from the sequence per AMP_UMI_REGEX
        expected_1 = "@MG01HS02:1483:HG7MTBCX3:1:1107:1230:2106_TATGGGCG 1:N:0:TAAGGCGA+NAGGCTTA"
        expected_2 = "CACTTGGCCAAGAGCTCACACTTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTGAAAAAANAAAAAGTAACAGTCAAATAGAAAGGTTAAATTT"

        with open(umi_extractor.r1_out_fastq, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):

                # Test that the UMI has moved to the qname
                if i == 0:
                    test_1 = line == expected_1

                # Test that the UMI sequence and adjacent adapter sequence were trimmed
                if i == 1:
                    test_2 = line == expected_2
                    break

        self.assertTrue(all((test_1, test_2)))

    def test_get_orig_primer_exists(self):
        """Test that we can find the primer in the R2 with perfect match."""

        umi_extractor = rp.UMIExtractor(r1_fastq=self.tileseq_r1_fastq, r2_fastq=self.tileseq_r2_fastq,
                                        umi_regex=TILESEQ_UMI_REGEX, primer_fasta=self.tileseq_umi_1r_fasta,
                                        outdir=self.tempdir)

        expected = TILESEQ_CBS_1R_FASTA[0:16]
        observed = umi_extractor.get_orig_r2_primer(r2_seq=TEST_R2_UMI_TILESEQ_FASTQ)

        self.assertEqual(expected, observed)

    def test_get_orig_primer_error(self):
        """Test that we can find the primer in the R2 with edit distance of 2."""

        umi_extractor = rp.UMIExtractor(r1_fastq=self.tileseq_r1_fastq, r2_fastq=self.tileseq_r2_fastq,
                                        umi_regex=TILESEQ_UMI_REGEX, primer_fasta=self.tileseq_umi_1r_fasta,
                                        outdir=self.tempdir)

        expected = TILESEQ_CBS_1R_FASTA[0:16]
        observed = umi_extractor.get_orig_r2_primer(r2_seq=TEST_R2_UMI_TILESEQ_PRIMER_ERROR_FASTQ)

        self.assertEqual(expected, observed)

    def test_get_orig_primer_nonexistent(self):
        """Test that we do not find a nonexistent primer in the R2."""

        umi_extractor = rp.UMIExtractor(r1_fastq=self.tileseq_r1_fastq, r2_fastq=self.tileseq_r2_fastq,
                                        umi_regex=TILESEQ_UMI_REGEX, primer_fasta=self.tileseq_umi_1r_fasta,
                                        outdir=self.tempdir)

        expected = "".join([rp.UMIExtractor.UNKNOWN_PRIMER_CHAR] * rp.UMIExtractor.PRIMER_SEQ_LEN)
        observed = umi_extractor.get_orig_r2_primer(r2_seq=TEST_R2_UMI_TILESEQ_NOPRIMER_FASTQ)

        self.assertEqual(expected, observed)

    def test_append_primer_name(self):
        """Tests that primers are appended to R1 and R2 names."""

        expected_1 = "@A00738:253:HF5HHDSX2:2:1101:13964:1063.CACAGGGGCTCCTTGG 1:N:0:TAATGCGC+AGGCTATA"
        expected_2 = "@A00738:253:HF5HHDSX2:2:1101:13964:1063.CACAGGGGCTCCTTGG 2:N:0:TAATGCGC+AGGCTATA"

        umi_extractor = rp.UMIExtractor(r1_fastq=self.tileseq_r1_fastq, r2_fastq=self.tileseq_r2_fastq,
                                        umi_regex=TILESEQ_UMI_REGEX, primer_fasta=self.tileseq_umi_1r_fasta,
                                        outdir=self.tempdir)

        with open(umi_extractor.r1_out_fastq, "r") as r1_out_fh, \
                open(umi_extractor.r2_out_fastq, "r") as r2_out_fh:

            for i, (r1_line, r2_line) in enumerate(zip(r1_out_fh, r2_out_fh)):

                # Test that the UMI has moved to the qname
                if i == 0:
                    test_1 = r1_line == expected_1
                    test_2 = r2_line == expected_1

                # Test that the UMI sequence and adjacent adapter sequence were trimmed
                if i == 1:
                    test_3 = r1_line == expected_1
                    test_4 = r2_line == expected_1
                    break

        self.assertTrue(all((test_1, test_2, test_3, test_4)))

# Skip testing ReadGrouper and ReadDeduplicator as these simply call umi_tools


class TestDeduplicators(unittest.TestCase):
    """Tests for read deduplicators."""

    @classmethod
    def setUpClass(cls):
        """Set up for TestDeduplicators."""

        cls.test_dir = os.path.dirname(__file__)
        cls.test_data_dir = os.path.join(cls.test_dir, "test_data")
        cls.group_bam = os.path.join(TEST_DATA_DIR, "CBS1_32_R.group.bam")
        cls.preprocess_bam = os.path.join(TEST_DATA_DIR, "CBS1_32_R.preprocess.bam")
        cls.dedup_bam = os.path.join(TEST_DATA_DIR, "CBS1_32_R.dedup.bam")

        # Write the test BAM which contains duplicated reads
        cls.ref = os.path.join(TEST_DATA_DIR, "CBS_pEZY3.fa")
        if not os.path.exists(cls.ref + ".fai"):
            pysam.faidx(cls.ref)

        with tempfile.NamedTemporaryFile(suffix=".test.sam", mode="w", delete=False) as test_sam:
            test_sam.write(DUP_TEST_BAM)
            cls.test_sam = test_sam.name

        cls.test_bam = su.sam_view(am=cls.test_sam)

        # Select a test read with a mismatch
        with pysam.AlignmentFile(cls.test_bam, "rb", check_sq=False) as test_af:
            for i, read in enumerate(test_af.fetch(until_eof=True)):
                if i == 0:
                    cls.test_align_seg_del = read
                if i == 13:
                    cls.test_align_seg_mismtach = read
            test_af.reset()

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestDeduplicators."""

        fu.safe_remove((cls.test_sam, cls.test_bam))


class TestConsensusDeduplicator(TestDeduplicators):
    """Tests for ConsensusDeduplicator."""

    @classmethod
    def setUpClass(cls):
        """Constructor for TestConsensusDeduplicator."""

        super(TestConsensusDeduplicator, cls).setUpClass()
        cls.cd = rp.ConsensusDeduplicator(in_bam=cls.test_bam, ref=cls.ref, contig_del_thresh=3)
        cls.consensus_quals = [0, 0, 0, 40, 40, None, 40, 40, 40, None, None, None, None, 40, 0, 0]
        cls.consensus_seq = ["A", "A", "A", "T", "T", "", "G", "G", "G", "", "", "", "", "T", "A", "A"]

    def test_extract_umi_network(self):
        """Tests the ability to extract the UMI group from BAM tag."""

        res = self.cd._extract_umi_network(align_seg=self.test_align_seg_mismtach)
        self.assertEqual(res, "1000004")

    def test_construct_cigar(self):
        """Tests construction of the CIGAR string with matches and a deletion."""

        # Returns a pysam cigartuple
        res = rp.ConsensusDeduplicator._construct_cigar(self.test_align_seg_del.query_alignment_qualities)
        expect = [(0, 6), (2, 1), (0, 63)]
        self.assertEqual(res, expect)

    def test_construct_align_seg(self):
        """Tests generation of a consensus read object from scratch."""

        res = self.cd._construct_align_seg(
            read_umi_network="0_R1", curr_mate_strand=rp.MATE_STRAND_POS_TUPLE(
                su.ReadMate("R1"), su.Strand("+"), pos=None, ref="dummy"), start_pos=99,
            consensus_seq=self.consensus_seq, consensus_quals=self.consensus_quals, n_duplicates=1)

        self.assertEquals(type(res), pysam.AlignedSegment)

    def test_get_consensus_read_attrs(self):
        """Tests for generation of a proper SAM flag."""

        res = self.cd._get_consensus_read_attrs(
            self.test_align_seg_del,
            rp.MATE_STRAND_POS_TUPLE(mate=su.ReadMate.R1, strand=su.Strand.MINUS, pos=None, ref=self.ref))

        # PAIRED,PROPER_PAIR,REVERSE,READ1
        expect = ("1000009", su.SAM_FLAG_PAIRED + su.SAM_FLAG_PROPER_PAIR + su.SAM_FLAG_REVERSE + su.SAM_FLAG_R1)

        self.assertEqual(res, expect)

    def test_get_missing_base_indices(self):
        """Tests that the proper indices are returned for missing bases."""

        res = self.cd._get_missing_base_indices(self.consensus_quals)
        expect = set(list(range(9,13)))
        self.assertEqual(res, expect)

    def test_set_missing_bases(self):
        """Tests that the proper bases are set to N."""

        res = self.cd._set_missing_bases(self.consensus_seq, self.consensus_quals)
        expect = (["A", "A", "A", "T", "T", "", "G", "G", "G", "N", "N", "N", "N", "T", "A", "A"],
                  [0,0,0,40,40,None,40,40,40,su.DEFAULT_MAX_BQ,su.DEFAULT_MAX_BQ,su.DEFAULT_MAX_BQ,su.DEFAULT_MAX_BQ,40,0,0])
        self.assertEqual(res, expect)
