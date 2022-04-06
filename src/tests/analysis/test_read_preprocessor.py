#!/usr/bin/env python3
"""Tests for analysis.read_preprocessor"""

import collections
import numpy as np
import pysam
import tempfile
import unittest

import analysis.read_preprocessor as rp
import analysis.seq_utils as su
import core_utils.file_utils as fu
from satmut_utils.definitions import *

tempfile.tempdir = DEFAULT_TEMPDIR

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

GROUP_TEST_SAM = """@HD	VN:1.0	SO:queryname
@SQ	SN:CBS_pEZY3	LN:7108
@RG	ID:CBS1_35_comb_R
MG01HS02:1483:HG7MTBCX3:1:1101:1237:33394_TTACCTGA	99	CBS_pEZY3	1971	44	131M	=	2074	254	CCTTTGCCCGCATGCTGATCGCGCAAGAGGGGCTGCTGTGCGGTGGCAGTGCTGGCAGCACGGTGGCGGTGGCCGTGAAGGCTGCGCAGGAGCTGCANGAGGGCCAGCGCTGCGTGGTCATTCTGCCCGAC	IHIIIHIIIIIIIIIIIIIIIIIIIIIIIHIIHIIHIIIHIIIIIIHHHHIIIHIIIHIIICHHHHIIIIIIIIIHIIIIIDGHIIIIIHIIIIIII#<DGHIIIHIIIIIIIIIIIIIIGIIIIIIIIHH	AS:i:259	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:97G33	YS:i:277	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:i:6956334	BX:Z:TTACCTGA
MG01HS02:1483:HG7MTBCX3:1:1101:1237:33394_TTACCTGA	147	CBS_pEZY3	2074	44	146M5S	=	1971	-254	CCAGCGCTGCGTGGTCATTCTGCCCGACTCAGTGCGGAACTACATGACCAAGTTCCTGAGCGACAGGTGGATGCTGCAGAAGGGCTTTCTGAAGGAGGAGGACCTCACGGAGAAGNAGCCCNGGTGNTGGCACCTCCGTGTNCANGNNNNN	HHIHDH?HCCFGHHHIIIHHHHIHIIIHHGGIIHHIIIIHIIGIIIIIHGIIHHCHIHHHHIIHIHIIIIIHHHHHHHIIIIIIHIIIIIIIIIIIIIIIIHHGHHEHIIIHE<<#IHG<<#FGD<#IHIHHECDIHHG<<#<<#<#####	AS:i:277	XN:i:0	XM:i:5	XO:i:0	XG:i:0	NM:i:5	MD:Z:115A5T4G14T2G1	YS:i:259	YT:Z:CP	RG:Z:CBS1_35_comb_R
MG01HS02:1483:HG7MTBCX3:1:1101:1225:82159_TTTTTCTA	83	CBS_pEZY3	1302	44	130M	=	1286	-156	AGTTCTTCAACGCGGGCGGGAGCGTGAAGGACCNCANCAGCCTGCGGATGATTGAGGATGCTGAGCGCGACGGGACGCTGAAGCCCGGGGACACGATTATCGAGCCGACATCCGGGAACACCGGGATCGG	@HHHGIHHHIIIIIIIIIIHIIIIHIIIHHD<<#<<#EIIIHHHHGIIHHHG@HIIIIIHIIIIIIIIIIHIHHIHIHIIIIIIGHIHIIIHHHHEHIHHIDIIIIIHDIIGIIIIIHIIIIIIHIHIHE	AS:i:254	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:33G2T93	YS:i:270	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:i:3275192	BX:Z:TTTTTCTA
MG01HS02:1483:HG7MTBCX3:1:1101:1225:82159_TTTTTCTA	163	CBS_pEZY3	1286	44	10S141M	=	1302	156	NNNNNNNGNNCNCTTGGCCAAGTGNGAGTNCTTCANCGCGGGCGGGAGCGTGAAGGACCGCATCAGCCTGCGGATGATTGAGGATGCTGAGCGCGACGGGACGCTGAAGCCCGGGGACACGATTATCGAGCCGACATCCGGGAACACCGGG	#######<##<#<<CEHIIIIIHG#<<EH#<<GHI#<<DHHIIIIHHHIICHIIIHIIHHHHHHIIIIHIIIIHHH1DHHHIHCH?HHEHHHICHHHHHIICHHHHHIIHIIIIIDEHHIIIIIIIFHHIIHHHIIGHIIIGHHHHH/DEH	AS:i:270	XN:i:0	XM:i:4	XO:i:0	XG:i:0	NM:i:4	MD:Z:1T12T4T5A115	YS:i:254	YT:Z:CP	RG:Z:CBS1_35_comb_R
"""

PREPROC_TEST_SAM = """@HD	VN:1.0	SO:unknown
@SQ	SN:CBS_pEZY3	LN:7108
@RG	ID:CBS1_35_comb_R
MG01HS02:1483:HG7MTBCX3:2:1213:4430:45262_ACTTTTGA	83	CBS_pEZY3	2290	44	105M	=	2290	-105	GGTGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAAGCCTGGGAATGGTGACGCTTGGGATCATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA	111<1<11CCGCC<<@1CC</C</0CCG<F<1C<<C/C<<H@G<C<<111<1<1CEC<000<<1HHF<CG<1D1<</0/<<<1/11FFD<</D@D1HFHF?GHIG	AS:i:198	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:2A44T23A33	YS:i:167	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:ACTTTTGA	UG:Z:10000001_R1
MG01HS02:1506:HGKGCBCX3:2:2205:7402:22634_ACTTTTGA	83	CBS_pEZY3	2290	44	105M	=	2290	-105	GGAGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAATCCTGGGAATGGTGACGCTTGGGAACATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA	IHIIIIHHIHIIIIIIIIIIIIHHHIGIIIHHIIHGIHHIIIHGEHDGHIIIIIIIIIIIIIHHIIHIIIIHIIIHHIHHIIIIHIIIHIIIIIIIIIIIIIIIH	AS:i:210	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:105	YS:i:210	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:ACTTTTGA	UG:Z:10000001_R1
MG01HS02:1483:HG7MTBCX3:2:1213:4430:45262_ACTTTTGA	163	CBS_pEZY3	2290	44	105M	=	2290	105	GGAGAAGGGCTTACAAGAGGCGCCCGAGGAGGATGCGGCGGGGGTAATCCTGGGAATGGTGAAGCTTGGGAACATGCTATCGTCCCTGCTTGCCGGGAAAGTGCA	0<00011<0111<1111<1D0/<C/<////0<0<<1/</<FEHH/0<CH11<1DH11DC111<DF1<1D?1D<11<1<<11<<10011111<1</C/C?11=1DD	AS:i:167	XN:i:0	XM:i:10	XO:i:0	XG:i:0	NM:i:10	MD:Z:12C0G1C0C9T2T5A26C15C20G5	YS:i:198	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10000001_R2
MG01HS02:1506:HGKGCBCX3:2:2205:7402:22634_ACTTTTGA	163	CBS_pEZY3	2290	44	105M	=	2290	105	GGAGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAATCCTGGGAATGGTGACGCTTGGGAACATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA	DDDCDIGIGHIIIIIIIIIIIIIIIIIGIIIIIHEHICHHHIHGHHI@GEHHII1D@HHIIIHHHIIIHHHHHHEGHEHHHII=HGIHHHICEGHIHEDEEGGHH	AS:i:210	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:105	YS:i:210	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10000001_R2
MG01HS02:1483:HG7MTBCX3:1:1116:9256:33330_ACTTTGGT	83	CBS_pEZY3	2290	44	105M	=	2290	-105	GGAGAAGGGCTTCGACCAGGCGCCCGTTGTGGATGAGGCGGGGGTAATCCTGGGAATGGTGACGCTTGGGAACATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA	GCC0FC0G?EC</E<?E</HCDE<<CD1HHGIIHE0C/<IIHIHECD<EEEHG<@@GD0=GEHG@GF1<1<1HHHC<0DHHHEE?C1?EHC?DHD1G1H@@FC1G	AS:i:206	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:27G77	YS:i:206	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:ATTTTGCT	UG:Z:10000004_R1
MG01HS02:1483:HG7MTBCX3:2:1214:5820:81145_ATTTTGCT	83	CBS_pEZY3	2290	42	105M	=	2290	-105	GGAGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAATCCTGGGAATGGTGACGCTTGGGATCATGCTCTCCTCCCTGCTTGCCGGGAAGGTGCA	HCIIHEIHCEHCDG<<1<</C/CC<1F<FD11?@CD</C/C<1GF1<<@C<<11CD110ED<<1?F<C@<11FF@1G<00<01/C<111C</C0C1D1HFFCC1F	AS:i:202	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:71A9G23	YS:i:150	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:ATTTTGCT	UG:Z:10000004_R1
MG01HS02:1483:HG7MTBCX3:2:2211:4802:70306_ATTTTGCT	83	CBS_pEZY3	2290	44	105M	=	2290	-105	GGAGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAATCCTGGGAATGGTGACGCTTGGGAACATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA	HIIIIIIIHIIIIIIHGIHHIIHECIIIIHHIIIIIIIIIIIHHIIIIIHHIIIIHHHHIHHHIHHHEIIHFIIHFHHHIHHFHHIIHHHHHH=F@IIIIHIIFG	AS:i:210	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:105	YS:i:210	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:ATTTTGCT	UG:Z:10000004_R1
MG01HS02:1506:HGKGCBCX3:2:1207:8268:71580_ACTTTGCT	83	CBS_pEZY3	2290	44	105M	=	2290	-105	GGAGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAATCCTGGGAATGGTGACGCTTGGGAACATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA	IIIIIIIIIIIIIIIIIIIIIIIHHIIFIIIIIIIIIIIIIIHIIIIIIHHECHHGIHGIIHHEGIHHHIGIIIHH<IHHHIHIHCIIIIIIIIIIIIIIIIIIH	AS:i:210	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:105	YS:i:210	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:ATTTTGCT	UG:Z:10000004_R1
MG01HS02:1483:HG7MTBCX3:1:1116:9256:33330_ACTTTGGT	163	CBS_pEZY3	2290	44	105M	=	2290	105	GGAGAAGGGCTTCGACCAGGCGCCCGTGATGGATGAGGCGGGGGTAATCCTGGGAATGGTGACGCTTGGGAACATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA	0<<0<DFC0111<CG<C/1C1<C/<?D=111<11<1<D<</EDE/<GC1D1<F?C1<CE@F?CCE<EH=11<<FH??<1D1<CEHHHCC<C1DHD<ECC00/<<C	AS:i:206	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:28G76	YS:i:206	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10000004_R2
MG01HS02:1483:HG7MTBCX3:2:1214:5820:81145_ATTTTGCT	163	CBS_pEZY3	2290	42	105M	=	2290	105	GGAGACCGGGTTCGACCAGGCCCCCGTGGTGGGTGAGCGGGGGGCAATCCGGGGGATGGAAACGCTTGGAAAAATGCTCTCGTCCCTGCTTGCCGGGGAGGTGCA	<<0<011/</0<</10/<E1D0<<</CEE11<//<<11///<///<<1110</D//=0<11<1<<00<01111111111<1<10<01<11111<DE//<//.1<D	AS:i:150	XN:i:0	XM:i:15	XO:i:0	XG:i:0	NM:i:15	MD:Z:5A0G2C11G10A4G0C5T5T3A4T0G8G2C24A7	YS:i:202	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10000004_R2
MG01HS02:1483:HG7MTBCX3:2:2211:4802:70306_ATTTTGCT	163	CBS_pEZY3	2290	44	105M	=	2290	105	GGAGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAATCCTGGGAATGGTGACGCTTGGGAACATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA	@D?B@CGI=ECCEHHIIIIFHIIGDHDHCGHHE<DCGHEHHHCDCGH?@H?D@H1<G1<1<CC@GGEHEGHGHFHEHHIIGDHIIIICHHIIHHHDGHHIIHIIF	AS:i:210	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:105	YS:i:210	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10000004_R2
MG01HS02:1506:HGKGCBCX3:2:1207:8268:71580_ACTTTGCT	163	CBS_pEZY3	2290	44	105M	=	2290	105	GGAGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAATCCTGGGAATGGTGACGCTTGGGAACATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA	DDDDCG??=HHIIHIIIIIIICCHDHF=DHHI1CDGHIGHIIIIHHHFHIGIGFHHFHI1FHC0EEGHCIHIGGIHI?FGCCGEGHE1GHHCHHHHIDHDFFHHE	AS:i:210	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:105	YS:i:210	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10000004_R2
MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC	99	CBS_pEZY3	2397	44	64M1D66M	=	2408	163	CGTCAGACCAAGTTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACGCGCTGGCAGGCTCTCGCACATCCTGGAGATGGACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTAC	HHHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHIIIHHHIIIIIIIHHIIIIHIIIGIHIIIIIIIIHHIIIIIIIIIIIIIIIIIIHIIIIIIIHIIHIHIIIEHHHIEHIIHIIC	AS:i:245	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:59A4^G66	YS:i:286	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:CGTTGATC	UG:Z:10015877_R1
MG01HS02:1483:HG7MTBCX3:1:2103:16630:87701_CGTTGATC	99	CBS_pEZY3	2397	44	64M1D66M	=	2408	163	CGTCAGACCAAGTTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACACGCTGGCAGGCTCTCGCACATCCTGGAGATGGACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTAC	HIHHHIIIIGIHHIIHIIIIHHIIIIIIIIIIIHIHGIIIIIIIIHIIIIIIIIIIIHIIIIIIHHHIIIIIHIIIIIIICHHGHGHIIIIGHIIIIIIIIIHGHHIIIGIIIIHHIIIIIHIHHHHIHD	AS:i:250	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:64^G66	YS:i:292	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:CGTTGATC	UG:Z:10015877_R1
MG01HS02:1483:HG7MTBCX3:2:1208:5073:30436_CGTTGATC	99	CBS_pEZY3	2397	44	64M1D66M	=	2409	163	CGTCAGACCAAGTTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACACGCTGGCAGGCTCTCGCACATCCTGGAGATGTACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTAC	HIIIIIIIIIIIIIIIIIIIIIIHHIIHIHIHHIIIIIIIIIIIIIIIIIHHIIIIIIIIIIIIIIIIHIIIIIIIIIIIGCGIHIIIGHIIIIIIIIIIIIIIIIIIIIIIHIIGHHIIHHHHHIHCH@	AS:i:244	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:64^G27G38	YS:i:285	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:CGTTGATC	UG:Z:10015877_R1
MG01HS02:1506:HGKGCBCX3:2:2205:19957:35646_CGTTGATC	99	CBS_pEZY3	2397	44	64M1D66M	=	2409	163	CGTCAGACCAAGTTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACACGCTGGCAGGCTCTCGCACATCCTGGAGATGGACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTAC	GHHIIHHIIIIIIIIIIIGIIHHIIIIIIIIIIIIHHIHIIIIHIIIIIIHIIIIIIIHHIIGDHHHIIIIIIIIIIIIIIIIGHIIEHIIIIIHHIIIIIIIIIIHIIIIIGIIHIIIIIIHIIIIIIC	AS:i:250	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:64^G66	YS:i:290	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:CGTTGATC	UG:Z:10015877_R1
MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC	147	CBS_pEZY3	2408	44	53M1D98M	=	2397	-163	GTTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACGCGCTGGCAGGCTCTCGCACATCCTGGAGATGGACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTACCACAGCACCGGGAAGTCCAGTCAGCGGCAGAT	CIHHHEHGIIIHHFFHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHIHIIHIIIIIIIIIIHIIIHIIIIIIHCIIIIIIIIHIIIIHIIIIIIIHHHIHIGIIIIGIIIIIIIIIIIHHDHGGIIIIIIGGIHIIHIIIIIDCDDD	AS:i:286	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:48A4^G98	YS:i:245	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10015877_R2
MG01HS02:1483:HG7MTBCX3:1:2103:16630:87701_CGTTGATC	147	CBS_pEZY3	2408	44	53M1D98M	=	2397	-163	GTTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACACGCTGGCAGGCTCTCGCACATCCTGGAGATGGACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTACCACAGCACCGGGAAGTCCAGTCAGCGGCAGAT	CHHHHHIHHHIHIHCHHIIIIIIIIIIIIIIIIIIIIIHIIIHHGIIHHEIIIIIIIIIIIIIIHIIIIIIIHHIIIIIIIIIHIHIGIIIIIIIIIIHIIIIIIHHHIHEGHGIIHHHIIIIIIHEHIIIIIHIHIGIIIIHIIIADBDD	AS:i:292	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:53^G98	YS:i:250	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10015877_R2
MG01HS02:1483:HG7MTBCX3:2:1208:5073:30436_CGTTGATC	147	CBS_pEZY3	2409	44	52M1D98M	=	2397	-163	TTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACACGCTGGCAGGCTCTCGCACATCCTGGAGATGTACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTACCACAGCACCGGGAAGTCCAGTCAGCGGCAGAT	HFCHHIHHIIHIHHIHHGIIIIIIIIIIIIHGIIGIHHHIHEHHFHHHHHIIHIIIIIIIIIIIIIIIHHEGIIHIIHIHIIIIIIIIIIHIIIIIIIHDFIIIIIIIHFHHEHIIIHHHHHIIIIIIIHIIIHHIIIIIGIIIIDDCDD	AS:i:285	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:52^G27G70	YS:i:244	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10015877_R2
MG01HS02:1506:HGKGCBCX3:2:2205:19957:35646_CGTTGATC	147	CBS_pEZY3	2409	44	52M1D98M	=	2397	-163	TTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACACGCTGGCAGGCTCTCGCACATCCTGGAGATGGACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTACCACAGCACCGGGAAGTCCAGTCAGCGGCAGAT	IHHIGIIHIIIIIHHHHHGIIHIIIIIHIIIGIHHIIIIFIIIIHHHGCIIIIIHIIIH?DHHIHHHCIIGIIIIIIIIIHIIIIIIIHEGIIIIIIIHHIIIIIIIIIIIHIIIIHHHEHGIIIIHIHEIIHGHIHIIHIIIIIDDDDD	AS:i:290	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:52^G98	YS:i:250	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10015877_R2
"""

TEST_PRIMERS = """CBS_pEZY3	2289	2309	CBSpEZY3_13F_GSP2	0	+
CBS_pEZY3	2318	2338	CBSpEZY3_12R_GSP2	0	-
CBS_pEZY3	2401	2428	CBSpEZY3_14F_GSP2	0	+
CBS_pEZY3	2430	2452	CBSpEZY3_13R_GSP2	0	-
CBS_pEZY3	2514	2536	CBSpEZY3_15F_GSP2	0	+
CBS_pEZY3	2538	2559	CBSpEZY3_14R_GSP2	0	-
"""

MASKING_TEST_PRIMERS = TEST_PRIMERS + """CBS_pEZY3	2380	2394	CBSpEZY3_mockR_GSP2	0	-
CBS_pEZY3	2289	2395	CBSpEZY3_containF_GSP2	0	+
CBS_pEZY3	2289	2394	CBSpEZY3_flushF_GSP2	0	+
CBS_pEZY3	2288	2394	CBSpEZY3_containR_GSP2	0	-
CBS_pEZY3	2289	2394	CBSpEZY3_flushR_GSP2	0	-
"""


class TestQnameVerification(unittest.TestCase):
    """Read name format verification for processing."""

    @classmethod
    def setUpClass(cls):
        """Set up for TestQnameVerification."""

        # Test both a valid FASTQ and BAM
        cls.tileseq_r1_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.R1.fastq", delete=False).name

        with open(cls.tileseq_r1_fastq, "w") as tileseq_r1_fh:
            tileseq_r1_fh.write(TEST_R1_TILESEQ_FASTQ)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".group.sam") as group_sam:
            group_sam.write(GROUP_TEST_SAM)
            fu.flush_files((group_sam,))
            cls.grouped_bam = su.sam_view(group_sam.name, None, "BAM", 0, "b")

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestQnameVerification."""

        fu.safe_remove((cls.tileseq_r1_fastq, cls.grouped_bam,))

    def test_verify_qname_format_Illumina(self):
        """Test that a standard Illumina qname format is verified."""

        expected = (True, ILLUMINA_FORMAT_INDEX)
        observed = rp.QnameVerification().verify_qname_format("MG01HS02:1483:HG7MTBCX3:1:1107:1230:2106")
        self.assertEqual(expected, observed)

    def test_verify_qname_format_int(self):
        """Test that an int qname format is verified."""

        expected = (True, INT_FORMAT_INDEX)
        observed = rp.QnameVerification().verify_qname_format("1")
        self.assertEqual(expected, observed)

    def test_verify_qname_format_invalid(self):
        """Test that an unrecognized qname format is not verified."""

        expected = (False, None)
        observed = rp.QnameVerification().verify_qname_format("strange_qname")
        self.assertEqual(expected, observed)

    def test_verify_fastq(self):
        """Test that qnames can be verified from a FASTQ."""

        expected = (True, ILLUMINA_FORMAT_INDEX)
        qv = rp.QnameVerification(fastq=self.tileseq_r1_fastq)
        observed = qv.compatible, qv.format_index
        self.assertEqual(expected, observed)

    def test_verify_bam(self):
        """Test that qnames can be verified from a BAM."""

        expected = (True, ILLUMINA_FORMAT_INDEX)
        qv = rp.QnameVerification(bam=self.grouped_bam)
        observed = qv.compatible, qv.format_index
        self.assertEqual(expected, observed)


class TestFastqPreprocessor(unittest.TestCase):
    """Tests for FastqPreprocessor."""

    @classmethod
    def setUpClass(cls):
        """Set up for TestFastqPreprocessor."""

        cls.tempdir = tempfile.mkdtemp()
        cls.tileseq_r1_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.R1.fastq", delete=False, dir=cls.tempdir).name
        cls.tileseq_r2_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.R2.fastq", delete=False, dir=cls.tempdir).name
        cls.tileseq_r1_umi_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.umi.R1.fastq", delete=False, dir=cls.tempdir).name
        cls.tileseq_r2_umi_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.umi.R2.fastq", delete=False, dir=cls.tempdir).name
        cls.amp_r1_fastq = tempfile.NamedTemporaryFile(suffix=".amp.umi.R1.fastq", delete=False, dir=cls.tempdir).name
        cls.amp_r2_fastq = tempfile.NamedTemporaryFile(suffix=".amp.umi.R2.fastq", delete=False, dir=cls.tempdir).name

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

        fu.safe_remove((cls.tempdir,), force_remove=True)

    def test_run_cutadapt_tileseq_r1_no_adapter(self):
        """Tests lack of adapter trimming of Tile-seq R1 without adapter readthrough."""

        expected = "GAGGAGGCGTTCACCTTTGCCCGCATGCTGATCGCGCAAGAGGGGCTGCTGTGCGGTGGCAGTGCTGGCAGCACGGTGGCGGTGGCCGTGAAGGCTGCGCAGGAGCTGCAGGAGGGCCAGCGCTGCGTGGTCATTCTGCCCGACTCAGTG"

        # For CDS-terminal PCR tiles we also include PEZY3 flanking sequences
        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   r2_fiveprime_adapters=",".join((NEB_ADAPTER_P5, PEZY3_ATTB2_P5)),
                                   r2_threeprime_adapters=",".join((NEB_ADAPTER_P7_RC, PEZY3_ATTB1_P7_RC)),
                                   outdir=self.tempdir)

        with open(fqp.trimmed_f1, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):
                if i == 1:
                    self.assertEqual(expected, line.strip(fu.FILE_NEWLINE))

    def test_run_cutadapt_tileseq_r1_normal_pair(self):
        """Tests adapter trimming of Tile-seq R1 with adapter readthrough."""

        expected = "CCAATTCTCACATCCTAGACCATCACGGGCATTGC"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   r2_fiveprime_adapters=",".join((NEB_ADAPTER_P5, PEZY3_ATTB2_P5)),
                                   r2_threeprime_adapters=",".join((NEB_ADAPTER_P7_RC, PEZY3_ATTB1_P7_RC)),
                                   outdir=self.tempdir)

        with open(fqp.trimmed_f1, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):
                if i == 5:
                    self.assertEqual(expected, line.strip(fu.FILE_NEWLINE))

    def test_run_cutadapt_tileseq_r1_abnormal_pair(self):
        r"""Tests adapter trimming of Tile-seq R1 with target priming off the 5' adapter, leading to abnormal
        5' adapter seq presence."""

        expected = "CCAATTCTCACATCCTATACCATCACGGGCATTGC"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   r2_fiveprime_adapters=",".join((NEB_ADAPTER_P5, PEZY3_ATTB2_P5)),
                                   r2_threeprime_adapters=",".join((NEB_ADAPTER_P7_RC, PEZY3_ATTB1_P7_RC)),
                                   outdir=self.tempdir)

        with open(fqp.trimmed_f1, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):
                if i == 9:
                    self.assertEqual(expected, line.strip(fu.FILE_NEWLINE))

    def test_run_cutadapt_tileseq_r2_no_adapter(self):
        """Tests lack of adapter trimming of Tile-seq R2 without adapter readthrough."""

        expected = "GCACTGAGTCGGGCAGAATGACCACGCAGCGCTGGCCCTCCTGCAGCTCCTGCGCAGCCTTCACGGCCACCGCCACCGTGCTGCCAGCACTGCCACCGCACAGCAGCCCCTCTTGCGCGATCAGCATGCGGGCAAAGGTGAACGCCTCCT"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   r2_fiveprime_adapters=",".join((NEB_ADAPTER_P5, PEZY3_ATTB2_P5)),
                                   r2_threeprime_adapters=",".join((NEB_ADAPTER_P7_RC, PEZY3_ATTB1_P7_RC)),
                                   outdir=self.tempdir)

        with open(fqp.trimmed_f2, "r") as r2_out_fh:
            for i, line in enumerate(r2_out_fh):
                if i == 1:
                    self.assertEqual(expected, line.strip(fu.FILE_NEWLINE))

    def test_run_cutadapt_tileseq_r2_normal_pair(self):
        """Tests lack of adapter trimming of Tile-seq R1 without adapter readthrough."""

        expected = "GCAATGCCCGTGATGGTCTAGGATGTGAGAATTGG"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   r2_fiveprime_adapters=",".join((NEB_ADAPTER_P5, PEZY3_ATTB2_P5)),
                                   r2_threeprime_adapters=",".join((NEB_ADAPTER_P7_RC, PEZY3_ATTB1_P7_RC)),
                                   outdir=self.tempdir)

        with open(fqp.trimmed_f2, "r") as r2_out_fh:
            for i, line in enumerate(r2_out_fh):
                if i == 5:
                    self.assertEqual(expected, line.strip(fu.FILE_NEWLINE))

    def test_run_cutadapt_tileseq_r2_abnormal_pair(self):
        r"""Tests adapter trimming of Tile-seq R2 with target priming off the R1 5' adapter, leading to abnormal
        R1 5' adapter seq presence."""

        expected = "GCAATGCCCGTGATGGTCTAGGATGTGAGAATTGG"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_fastq, f2=self.tileseq_r2_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   r2_fiveprime_adapters=",".join((NEB_ADAPTER_P5, PEZY3_ATTB2_P5)),
                                   r2_threeprime_adapters=",".join((NEB_ADAPTER_P7_RC, PEZY3_ATTB1_P7_RC)),
                                   outdir=self.tempdir)

        with open(fqp.trimmed_f2, "r") as r2_out_fh:
            for i, line in enumerate(r2_out_fh):
                if i == 9:
                    self.assertEqual(expected, line.strip(fu.FILE_NEWLINE))

    def test_run_cutadapt_tileseq_r1_umi(self):
        """Tests adapter trimming of Tile-seq UMI-containing R1 with adapter readthrough."""

        expected = "ATGCCTTCTGAGACCCCCCAGGCAGAAGTGGGGCCCACAGGCTGCCCCCACCGCTCAGGGCCACACTCGGCGAAGGGGAGCCTGGAGAAGGGGTCCCCAGA"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_umi_fastq, f2=self.tileseq_r2_umi_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   r2_fiveprime_adapters=",".join((NEB_ADAPTER_P5, PEZY3_ATTB2_P5)),
                                   r2_threeprime_adapters=",".join((NEB_ADAPTER_P7_RC, PEZY3_ATTB1_P7_RC)),
                                   outdir=self.tempdir)

        with open(fqp.trimmed_f1, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):
                if i == 1:
                    self.assertEqual(expected, line.strip(fu.FILE_NEWLINE))

    def test_run_cutadapt_tileseq_r2_umi(self):
        """Tests adapter trimming of Tile-seq R2 with adapter readthrough paired with UMI-containing R1."""

        expected = "CACAGGGGCTCCTTGGCTTCCTTATCCTCTGGGGACCCCTTCTCCAGGCTCCCCTTCGCCGAGTGTGGCCCTGAGCGGTGGGGGCAGCCTGTGGGCCCCACTTCTGCCTGGGGGGTCTCAGAAGGCAT"

        fqp = rp.FastqPreprocessor(f1=self.tileseq_r1_umi_fastq, f2=self.tileseq_r2_umi_fastq,
                                   r1_fiveprime_adapters=",".join((NEB_ADAPTER_P7, PEZY3_ATTB1_P7)),
                                   r1_threeprime_adapters=",".join((NEB_ADAPTER_P5_RC, PEZY3_ATTB2_P5_RC)),
                                   r2_fiveprime_adapters=",".join((NEB_ADAPTER_P5, PEZY3_ATTB2_P5)),
                                   r2_threeprime_adapters=",".join((NEB_ADAPTER_P7_RC, PEZY3_ATTB1_P7_RC)),
                                   outdir=self.tempdir)

        with open(fqp.trimmed_f2, "r") as r2_out_fh:
            for i, line in enumerate(r2_out_fh):
                if i == 1:
                    self.assertEqual(expected, line.strip(fu.FILE_NEWLINE))

    def test_run_cutadapt_amp_r1_normal_pair(self):
        """Tests adapter trimming of AMP R1 containing 5' and 3' adapters."""

        # The anchored common region should be trimmed at the 5' end, and the GSP2 tail RC trimmed at the 3' end
        expected = "CACTTGGCCAAGAGCTCACACTTC"

        fqp = rp.FastqPreprocessor(f1=self.amp_r1_fastq, f2=self.amp_r2_fastq,
                                   r1_fiveprime_adapters=AMP_CR, r1_threeprime_adapters=AMP_GSP2_TAIL_RC,
                                   r2_fiveprime_adapters=AMP_GSP2_TAIL, r2_threeprime_adapters=AMP_CR_RC,
                                   outdir=self.tempdir)

        with open(fqp.trimmed_f1, "r") as r1_out_fh:
            for i, line in enumerate(r1_out_fh):
                if i == 1:
                    self.assertEqual(expected, line.strip(fu.FILE_NEWLINE))

    def test_run_cutadapt_amp_r2_normal_pair(self):
        """Tests adapter trimming of AMP R2 containing 3' adapter."""

        # The anchored common region RC should be trimmed at the 3' end
        expected = "NNNNNGTGAGCTCTTGGCCAAGTG"

        fqp = rp.FastqPreprocessor(f1=self.amp_r1_fastq, f2=self.amp_r2_fastq,
                                   r1_fiveprime_adapters=AMP_CR, r1_threeprime_adapters=AMP_GSP2_TAIL_RC,
                                   r2_fiveprime_adapters=AMP_GSP2_TAIL, r2_threeprime_adapters=AMP_CR_RC,
                                   outdir=self.tempdir)

        with open(fqp.trimmed_f2, "r") as r2_out_fh:
            for i, line in enumerate(r2_out_fh):
                if i == 1:
                    self.assertEqual(expected, line.strip(fu.FILE_NEWLINE))


class TestUmiExtractor(unittest.TestCase):
    """Tests for UmiExtractor."""

    @classmethod
    def setUpClass(cls):
        """Set up for TestUmiExtractor."""

        cls.tempdir = tempfile.mkdtemp()
        cls.tileseq_r1_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.umi.R1.fastq", delete=False, dir=cls.tempdir).name
        cls.tileseq_r2_fastq = tempfile.NamedTemporaryFile(suffix=".tileseq.umi.R2.fastq", delete=False, dir=cls.tempdir).name
        cls.amp_r1_fastq = tempfile.NamedTemporaryFile(suffix=".amp.umi.R1.fastq", delete=False, dir=cls.tempdir).name
        cls.amp_r2_fastq = tempfile.NamedTemporaryFile(suffix=".amp.umi.R2.fastq", delete=False, dir=cls.tempdir).name
        cls.tileseq_primer_1r_fasta = tempfile.NamedTemporaryFile(suffix=".tileseq.primer.1R.fasta", delete=False, dir=cls.tempdir).name

        with open(cls.tileseq_r1_fastq, "w") as tileseq_r1_fh, \
                open(cls.tileseq_r2_fastq, "w") as tileseq_r2_fh, \
                open(cls.amp_r1_fastq, "w") as amp_r1_fh, \
                open(cls.amp_r2_fastq, "w") as amp_r2_fh, \
                open(cls.tileseq_primer_1r_fasta, "w") as tileseq_umi_1r_fh:

            tileseq_r1_fh.write(TEST_R1_UMI_TILESEQ_FASTQ)
            tileseq_r2_fh.write(TEST_R2_UMI_TILESEQ_FASTQ)
            amp_r1_fh.write(TEST_R1_UMI_AMP_FASTQ)
            amp_r2_fh.write(TEST_R2_UMI_AMP_FASTQ)
            tileseq_umi_1r_fh.write(TILESEQ_CBS_1R_FASTA)

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestUmiExtractor."""

        fu.safe_remove((cls.tempdir,), force_remove=True)

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
                    test_1 = line.strip(fu.FILE_NEWLINE) == expected_1

                # Test that the UMI sequence and adjacent adapter sequence were trimmed
                if i == 1:
                    test_2 = line.strip(fu.FILE_NEWLINE) == expected_2
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
                    test_1 = line.strip(fu.FILE_NEWLINE) == expected_1

                # Test that the UMI sequence and adjacent adapter sequence were trimmed
                if i == 1:
                    test_2 = line.strip(fu.FILE_NEWLINE) == expected_2
                    break

        self.assertTrue(all((test_1, test_2)))

    def test_get_orig_primer_exists(self):
        """Test that we can find the primer in the R2 with perfect match."""

        umi_extractor = rp.UMIExtractor(r1_fastq=self.tileseq_r1_fastq, r2_fastq=self.tileseq_r2_fastq,
                                        umi_regex=TILESEQ_UMI_REGEX, primer_fasta=self.tileseq_primer_1r_fasta,
                                        outdir=self.tempdir)

        expected = "CAGGGGCTCCTTGGCT"
        observed = umi_extractor.get_orig_r2_primer(r2_seq=TEST_R2_UMI_TILESEQ_FASTQ.splitlines()[1])

        self.assertEqual(expected, observed)

    def test_get_orig_primer_error(self):
        """Test that we can find the primer in the R2 with edit distance of 2."""

        umi_extractor = rp.UMIExtractor(r1_fastq=self.tileseq_r1_fastq, r2_fastq=self.tileseq_r2_fastq,
                                        umi_regex=TILESEQ_UMI_REGEX, primer_fasta=self.tileseq_primer_1r_fasta,
                                        outdir=self.tempdir)

        expected = "CAGGGGCTCCTTGGCT"
        observed = umi_extractor.get_orig_r2_primer(r2_seq=TEST_R2_UMI_TILESEQ_PRIMER_ERROR_FASTQ.splitlines()[1])

        self.assertEqual(expected, observed)

    def test_get_orig_primer_nonexistent(self):
        """Test that we do not find a nonexistent primer in the R2."""

        umi_extractor = rp.UMIExtractor(r1_fastq=self.tileseq_r1_fastq, r2_fastq=self.tileseq_r2_fastq,
                                        umi_regex=TILESEQ_UMI_REGEX, primer_fasta=self.tileseq_primer_1r_fasta,
                                        outdir=self.tempdir)

        expected = "".join([rp.UMIExtractor.UNKNOWN_PRIMER_CHAR] * rp.UMIExtractor.PRIMER_SEQ_LEN)
        observed = umi_extractor.get_orig_r2_primer(r2_seq=TEST_R2_UMI_TILESEQ_NOPRIMER_FASTQ)

        self.assertEqual(expected, observed)

    def test_append_primer_name(self):
        """Tests that primers are appended to R1 and R2 names."""

        expected_1 = "@A00738:253:HF5HHDSX2:2:1101:13964:1063.CAGGGGCTCCTTGGCT_GCTCGGCG 1:N:0:TAATGCGC+AGGCTATA"
        expected_2 = "@A00738:253:HF5HHDSX2:2:1101:13964:1063.CAGGGGCTCCTTGGCT_GCTCGGCG 2:N:0:TAATGCGC+AGGCTATA"

        umi_extractor = rp.UMIExtractor(r1_fastq=self.tileseq_r1_fastq, r2_fastq=self.tileseq_r2_fastq,
                                        umi_regex=TILESEQ_UMI_REGEX, primer_fasta=self.tileseq_primer_1r_fasta,
                                        outdir=self.tempdir)

        with open(umi_extractor.r1_out_fastq, "r") as r1_out_fh, \
                open(umi_extractor.r2_out_fastq, "r") as r2_out_fh:

            for i, (r1_line, r2_line) in enumerate(zip(r1_out_fh, r2_out_fh)):

                # Test that the UMI has moved to the qname
                if i == 0:
                    test_1 = r1_line.strip(fu.FILE_NEWLINE) == expected_1
                    test_2 = r2_line.strip(fu.FILE_NEWLINE) == expected_2

                break

        self.assertTrue(all((test_1, test_2,)))

# Skip testing ReadGrouper and ReadDeduplicator as these are wrappers of umi_tools


class TestConsensusDeduplicatorPreprocessor(unittest.TestCase):
    """Tests for ConsensusDeduplicatorPreprocessor."""

    @classmethod
    def setUpClass(cls):
        """Set up for TestConsensusDeduplicatorPreprocessor."""

        cls.tempdir = tempfile.mkdtemp()

        with tempfile.NamedTemporaryFile(mode="w", suffix=".group.sam") as group_sam:
            group_sam.write(GROUP_TEST_SAM)
            fu.flush_files((group_sam,))
            cls.grouped_bam = su.sam_view(group_sam.name, None, "BAM", 0, "b")

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestConsensusDeduplicator."""

        fu.safe_remove((cls.tempdir, cls.grouped_bam,), force_remove=True)

    def test_update_tags(self):
        """Tests that we properly update the tags in a grouped BAM with only R1 marked with group tag."""

        expected = ["6956334_R1", "6956334_R2", "3275192_R1", "3275192_R2"]

        updated_bam = rp.ConsensusDeduplicatorPreprocessor.update_tags(qname_sorted=self.grouped_bam)

        with pysam.AlignmentFile(updated_bam, "rb") as updated_af:
            observed = [align_seg.get_tag(rp.UMITOOLS_UG_TAG) for align_seg in updated_af.fetch(until_eof=True)]

        fu.safe_remove((updated_bam,))
        self.assertEqual(expected, observed)

    def test_update_tags_from_grouped_bam(self):
        """Tests that we properly qname-sorted the BAM and update tags."""

        expected = ["3275192_R1", "3275192_R2", "6956334_R1", "6956334_R2"]

        sorted_updated_bam = rp.ConsensusDeduplicatorPreprocessor.update_tags_from_grouped_bam(in_bam=self.grouped_bam)

        with pysam.AlignmentFile(sorted_updated_bam, "rb") as updated_af:
            observed = [align_seg.get_tag(rp.UMITOOLS_UG_TAG) for align_seg in updated_af.fetch(until_eof=True)]

        fu.safe_remove((sorted_updated_bam,))
        self.assertEqual(expected, observed)

    def test_workflow(self):
        """Tests that we output a BAM with final sorting on group tag."""

        expected = ["3275192_R1", "3275192_R2", "6956334_R1", "6956334_R2"]

        cdp = rp.ConsensusDeduplicatorPreprocessor(group_bam=self.grouped_bam, outdir=self.tempdir)

        with pysam.AlignmentFile(cdp.preprocess_bam, "rb") as preprocess_af:
            observed = [align_seg.get_tag(rp.UMITOOLS_UG_TAG) for align_seg in preprocess_af.fetch(until_eof=True)]

        fu.safe_remove((cdp.preprocess_bam,))
        self.assertEqual(expected, observed)


class TestConsensusDeduplicator(unittest.TestCase):
    """Tests for ConsensusDeduplicator."""

    @classmethod
    def setUpClass(cls):
        """Constructor for TestConsensusDeduplicator."""

        cls.tempdir = tempfile.mkdtemp()
        cls.test_dir = os.path.dirname(__file__)
        cls.test_data_dir = os.path.join(os.path.split(cls.test_dir)[0], "test_data")
        cls.ref = os.path.join(cls.test_data_dir, "CBS_pEZY3.fa")

        with tempfile.NamedTemporaryFile(mode="w", suffix=".preprocess.sam") as preproc_sam:
            preproc_sam.write(PREPROC_TEST_SAM)
            fu.flush_files((preproc_sam,))
            cls.preproc_bam = su.sam_view(preproc_sam.name, None, "BAM", 0, "b")

        cls.cd = rp.ConsensusDeduplicator(in_bam=cls.preproc_bam, ref=cls.ref, outdir=cls.tempdir, contig_del_thresh=3)

        # Create simpler mock data for certain methods
        # Select reads for various method testing
        with pysam.AlignmentFile(cls.preproc_bam, "rb", check_sq=False) as test_af:
            cls.test_header = test_af.header

            for i, read in enumerate(test_af.fetch(until_eof=True)):
                if i == 0:
                    cls.test_align_seg_mismatches = read
                if i == 1:
                    cls.test_align_seg_clean = read
                if i == 12:
                    cls.test_align_seg_del = read
                if i == 13:
                    cls.test_align_seg_del_dup = read
                if i == 14:
                    cls.test_align_seg_del_dup_minority_call = read
                if i == 17:
                    cls.test_align_seg_del_r2_a = read

            test_af.reset()

        # For testing R2 merging into contigs and del_threshold
        cls.consensus_quals = [0, 0, 0, 40, 40, None, 40, 40, 40, None, None, None, None, 40, 0, 0]
        cls.consensus_seq = ["A", "A", "A", "T", "T", "", "G", "G", "G", "", "", "", "", "T", "A", "A"]

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestConsensusDeduplicator."""

        fu.safe_remove((cls.tempdir, cls.preproc_bam, cls.cd.out_bam), force_remove=True)

    def test_extract_umi_network(self):
        """Tests the ability to extract the UMI group from BAM tag."""

        res = self.cd._extract_umi_network(align_seg=self.test_align_seg_mismatches)
        self.assertEqual(res, "10000001")

    def test_init_positions(self):
        """Tests that we properly initialize positions for a dedup contig."""

        # 0-based start compared to SAM 1-based start
        expected = [rp.MATE_STRAND_POS_TUPLE(
            mate=su.ReadMate("R1"), strand=su.Strand("-"), pos=i, ref="CBS_pEZY3") for i in range(2289, 2289+105)]

        consensus_dict = collections.OrderedDict()

        test_mate_strand = rp.MATE_STRAND_POS_TUPLE(
            mate=su.ReadMate(self.test_align_seg_mismatches.is_read1),
            strand=su.Strand(self.test_align_seg_mismatches.is_reverse),
            pos=None, ref=self.test_align_seg_mismatches.reference_name)

        self.cd._init_positions(align_seg=self.test_align_seg_mismatches,
                                mate_strand=test_mate_strand, consensus_dict=consensus_dict)

        # Just check the keys as checking empty arrays gives a ValueError from unittest
        observed = list(consensus_dict.keys())

        self.assertEqual(expected, observed)

    def test_construct_cigar(self):
        """Tests construction of the CIGAR string with matches and a deletion."""

        # Generate a consensus_quals list with None at the del position
        quals = list(self.test_align_seg_del.query_alignment_qualities)
        quals.insert(64, None)

        # Returns a pysam cigartuple
        observed = rp.ConsensusDeduplicator._construct_cigar(quals)
        expected = [(0, 64), (2, 1), (0, 66)]
        self.assertEqual(expected, observed)

    def test_construct_align_seg(self):
        """Tests generation of a consensus read object from scratch."""

        res = self.cd._construct_align_seg(
            read_umi_network="0_R1", curr_mate_strand=rp.MATE_STRAND_POS_TUPLE(
                su.ReadMate("R1"), su.Strand("+"), pos=None, ref="dummy"), start_pos=99,
            consensus_seq=self.consensus_seq, consensus_quals=self.consensus_quals, n_duplicates=1)

        self.assertEquals(type(res), pysam.AlignedSegment)

    def test_get_missing_base_indices(self):
        """Tests that the proper indices are returned for missing bases."""

        observed = self.cd._get_missing_base_indices(self.consensus_quals)
        expected = set(range(9, 13))
        self.assertEqual(expected, observed)

    def test_set_missing_bases(self):
        """Tests that the proper bases are set to N."""

        observed = self.cd._set_missing_bases(self.consensus_seq, self.consensus_quals)
        expected = (["A", "A", "A", "T", "T", "", "G", "G", "G", "N", "N", "N", "N", "T", "A", "A"],
                  [0,0,0,40,40,None,40,40,40,su.DEFAULT_MAX_BQ,su.DEFAULT_MAX_BQ,su.DEFAULT_MAX_BQ,su.DEFAULT_MAX_BQ,40,0,0])
        self.assertEqual(expected, observed)

    def test_update_consensus_dict(self):
        """Tests that we add read information to the consensus dictionary."""

        # Test that we add bases to the consensus dict and the pos dict contains the start position
        # 64 match, 1 del, 66 match; sum is 131 ref positions in the consensus dict
        dict_values = [np.zeros(shape=10, dtype=np.int32)] * 130

        for i, base in enumerate(self.test_align_seg_del.query_sequence):
            idx = self.cd._index_from_base(base)
            dict_values[i][idx] += 1

        # Finally insert the del base into the consensus dict
        dict_values.insert(64, np.zeros(shape=10, dtype=np.int32))

        expected_1 = 1
        expected_2 = 0
        expected_3 = [2396]

        consensus_dict = collections.OrderedDict()
        pos_list = []

        # Updates the dict and set by side effect
        self.cd._update_consensus_dict(self.test_align_seg_del, consensus_dict, pos_list)

        # Sum the base counts (omit qualities which are [5:]
        observed_1 = sum(consensus_dict[rp.MATE_STRAND_POS_TUPLE(
            mate=su.ReadMate("R1"), strand=su.Strand("+"), pos=2396, ref="CBS_pEZY3")][0:5])

        observed_2 = sum(consensus_dict[rp.MATE_STRAND_POS_TUPLE(
            mate=su.ReadMate("R1"), strand=su.Strand("+"), pos=2460, ref="CBS_pEZY3")][0:5])

        test_res = [expected_1 == observed_1, expected_2 == observed_2, expected_3 == pos_list]
        self.assertTrue(all(test_res))

    def test_get_consensus_del(self):
        """Tests that a consensus base at a position is properly generated."""

        expected = ("", None)

        consensus_dict = collections.OrderedDict()
        pos_list = []
        self.cd._update_consensus_dict(self.test_align_seg_del, consensus_dict, pos_list)
        self.cd._update_consensus_dict(self.test_align_seg_del_dup, consensus_dict, pos_list)

        # We have a deletion in both reads at this index
        consensus_key, consensus_array = list(consensus_dict.items())[64]
        observed = self.cd._get_consensus(consensus_key, consensus_array)

        self.assertEqual(expected, observed)

    def test_get_consensus_one_match_ref(self):
        r"""Tests that a consensus base at a position is called as the reference base with two supporting bases and \
         only one matching the reference base."""

        # The read with the base matching the reference has BQ 40
        expected = ("A", 40)

        consensus_dict = collections.OrderedDict()
        pos_list = []
        self.cd._update_consensus_dict(self.test_align_seg_del, consensus_dict, pos_list)
        self.cd._update_consensus_dict(self.test_align_seg_del_dup, consensus_dict, pos_list)

        # Here we have a mismatch
        consensus_key, consensus_array = list(consensus_dict.items())[59]
        observed = self.cd._get_consensus(consensus_key, consensus_array)

        self.assertEqual(expected, observed)

    def test_get_consensus_majority(self):
        r"""Tests that a consensus base at a position is called as the reference base with two supporting bases and \
         only one matching the reference base."""

        # The two BQs are 39 and 40, with int() taking the floor of the mean
        expected = ("G", 40)

        consensus_dict = collections.OrderedDict()
        pos_list = []
        self.cd._update_consensus_dict(self.test_align_seg_del, consensus_dict, pos_list)
        self.cd._update_consensus_dict(self.test_align_seg_del_dup, consensus_dict, pos_list)
        self.cd._update_consensus_dict(self.test_align_seg_del_dup_minority_call, consensus_dict, pos_list)

        # Here we have a mismatch
        consensus_key, consensus_array = list(consensus_dict.items())[91]
        observed = self.cd._get_consensus(consensus_key, consensus_array)

        self.assertEqual(expected, observed)

    def test_get_consensus_read_attrs(self):
        """Tests that we get the consensus read SAM flags and new qname."""

        expected = ("10000001", 83)

        mate_strand = rp.MATE_STRAND_POS_TUPLE(mate=su.ReadMate("R1"), strand=su.Strand("-"), pos=None, ref="CBS_pEZY3")
        observed = self.cd._get_consensus_read_attrs("10000001_R1", mate_strand)

        self.assertEqual(expected, observed)

    def test_write_consensus_one_match_ref(self):
        """Tests that a consensus read is generated where one duplicate has mismatches and the other matches the ref."""

        # Test that the consensus read has no mismatches, that it starts at the expected coordinate, and that
        # the number duplicates (ND) tag is 2
        expected_1 = self.test_align_seg_clean.query_alignment_sequence
        expected_2 = self.test_align_seg_clean.query_alignment_start
        expected_3 = 2

        consensus_dict = collections.OrderedDict()
        pos_list = []
        self.cd._update_consensus_dict(self.test_align_seg_mismatches, consensus_dict, pos_list)
        self.cd._update_consensus_dict(self.test_align_seg_clean, consensus_dict, pos_list)

        # Write the consensus read
        with tempfile.NamedTemporaryFile(suffix=".consensus.bam", delete=False, dir=self.tempdir) as consensus_bam, \
                pysam.AlignmentFile(consensus_bam, mode="wb", header=self.test_header) as consensus_af:

            self.cd._write_consensus(consensus_af, consensus_dict, pos_list, "10000001_R1")
            consensus_bam_name = consensus_bam.name

        with pysam.AlignmentFile(consensus_bam_name, "rb") as res_af:
            for align_seg in res_af.fetch(until_eof=True):
                observed_1 = align_seg.query_alignment_sequence
                observed_2 = align_seg.query_alignment_start
                observed_3 = align_seg.get_tag(self.cd.N_DUPLICATES_TAG)

        self.assertTrue(all((expected_1 == observed_1, expected_2 == observed_2, expected_3 == observed_3)))

    def test_write_consensus_majority_with_del(self):
        """Tests that a consensus read is generated by majority vote while retaining a deletion."""

        # Test that the consensus read has no mismatches but retains the del, that it starts at the expected coordinate,
        # and that the number duplicates (ND) tag is 3
        expected_1 = self.test_align_seg_del_dup.query_alignment_sequence
        expected_2 = self.test_align_seg_del_dup.query_alignment_start
        expected_3 = 3

        consensus_dict = collections.OrderedDict()
        pos_list = []
        self.cd._update_consensus_dict(self.test_align_seg_del, consensus_dict, pos_list)
        self.cd._update_consensus_dict(self.test_align_seg_del_dup, consensus_dict, pos_list)
        self.cd._update_consensus_dict(self.test_align_seg_del_dup_minority_call, consensus_dict, pos_list)

        # Write the consensus read
        with tempfile.NamedTemporaryFile(suffix=".consensus.bam", delete=False, dir=self.tempdir) as consensus_bam, \
                pysam.AlignmentFile(consensus_bam, mode="wb", header=self.test_header) as consensus_af:

            self.cd._write_consensus(consensus_af, consensus_dict, pos_list, "10015877_R1")
            consensus_bam_name = consensus_bam.name

        with pysam.AlignmentFile(consensus_bam_name, "rb") as res_af:
            for align_seg in res_af.fetch(until_eof=True):
                observed_1 = align_seg.query_alignment_sequence
                observed_2 = align_seg.query_alignment_start
                observed_3 = align_seg.get_tag(self.cd.N_DUPLICATES_TAG)

        fu.safe_remove((consensus_bam_name,))

        test_res = [expected_1 == observed_1, expected_2 == observed_2, expected_3 == observed_3]

        self.assertTrue(all(test_res))

    def test_generate_consensus_reads(self):
        """Tests that three fragments are properly deduplicated in succession."""

        # All reads starting at 2290 should match the reference after deduplication
        expected_1 = [self.test_align_seg_clean.query_alignment_sequence] * 4

        # All reads starting at 2397 should match the reference with the exception of the del after dedup
        expected_1.append(self.test_align_seg_del_dup.query_alignment_sequence)

        # R2 duplicates start at 2408 and 2409; we expect the start for the consensus read to be 2408
        expected_1.append(self.test_align_seg_del_r2_a.query_alignment_sequence)

        expected_2 = ["10000001", "10000001", "10000004", "10000004", "10015877", "10015877"]

        consensus_bam = self.cd._generate_consensus_reads()

        observed_1 = []
        observed_2 = []
        with pysam.AlignmentFile(consensus_bam, "rb") as consensus_af:
            for align_seg in consensus_af.fetch(until_eof=True):
                observed_1.append(align_seg.query_alignment_sequence)
                observed_2.append(align_seg.query_name)

        fu.safe_remove((consensus_bam,))

        test_res = [expected_1 == observed_1, expected_2 == observed_2]
        self.assertTrue(all(test_res))


class TestReadMasker(unittest.TestCase):
    """Tests for ReadMasker."""

    @classmethod
    def setUpClass(cls):
        """Constructor for TestReadMasker."""

        cls.tempdir = tempfile.mkdtemp()
        cls.test_dir = os.path.dirname(__file__)
        cls.test_data_dir = os.path.join(os.path.split(cls.test_dir)[0], "test_data")

        with tempfile.NamedTemporaryFile(mode="w", suffix=".preprocess.sam") as preproc_sam, \
                tempfile.NamedTemporaryFile(mode="w", suffix=".primers.bed", delete=False, dir=cls.tempdir) as primer_bed, \
                tempfile.NamedTemporaryFile(mode="w", suffix=".primers2.bed", delete=False, dir=cls.tempdir) as primer_bed2:

            preproc_sam.write(PREPROC_TEST_SAM)
            fu.flush_files((preproc_sam,))
            cls.preproc_bam = su.sam_view(preproc_sam.name, None, "BAM", 0, "b")

            primer_bed.write(MASKING_TEST_PRIMERS)
            cls.primer_bed = primer_bed.name

            primer_bed2.write(TEST_PRIMERS)
            cls.primer_bed2 = primer_bed2.name

        # Select reads for various method testing
        with pysam.AlignmentFile(cls.preproc_bam, "rb") as test_af:
            for i, read in enumerate(test_af.fetch(until_eof=True)):
                if i == 0:
                    cls.test_align_seg_r1_reverse = read
                if i == 2:
                    cls.test_align_seg_r2_positive = read
                if i > 2:
                    break

            test_af.reset()

        cls.rm_tileseq = rp.ReadMasker(in_bam=cls.preproc_bam, feature_file=cls.primer_bed, outdir=cls.tempdir)

        cls.rm_race = rp.ReadMasker(
            in_bam=cls.preproc_bam, feature_file=cls.primer_bed, race_like=True, outdir=cls.tempdir)

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestReadMasker."""

        fu.safe_remove((cls.tempdir, cls.preproc_bam,), force_remove=True)

    def test_get_mask_base_indices_tileseq_r1(self):
        """Test that we return the read indices to mask for a Tile-seq R1 starting and ending at a primer start."""

        expected = set(list(range(0, 20)) + list(range(91, 105)))

        associated_primers = {"CBS_pEZY3:2289-2309:+", "CBS_pEZY3:2380-2394:-"}
        observed = self.rm_tileseq._get_mask_base_indices(
            align_seg=self.test_align_seg_r1_reverse, associated_primers=associated_primers)

        self.assertEqual(0, len(expected - observed))

    def test_get_mask_base_indices_tileseq_r2(self):
        """Test that we return the read indices to mask for a Tile-seq R2 starting and ending at a primer start."""

        expected = set(list(range(0, 20)) + list(range(91, 105)))

        associated_primers = {"CBS_pEZY3:2289-2309:+", "CBS_pEZY3:2380-2394:-"}
        observed = self.rm_tileseq._get_mask_base_indices(
            align_seg=self.test_align_seg_r2_positive, associated_primers=associated_primers)

        self.assertEqual(0, len(expected - observed))

    def test_get_mask_base_indices_race_r1_end(self):
        """Test that we return the read indices to mask for a RACE-like R1 ending at a primer."""

        expected = set(range(0, 20))

        associated_primers = {"CBS_pEZY3:2289-2309:+", "CBS_pEZY3:2380-2394:-"}
        observed = self.rm_race._get_mask_base_indices(
            align_seg=self.test_align_seg_r1_reverse, associated_primers=associated_primers)

        self.assertEqual(0, len(expected - observed))

    def test_get_mask_base_indices_race_r2_start(self):
        """Test that we return the read indices to mask for a RACE-like R2 starting at a primer."""

        expected = set(range(0, 20))

        associated_primers = {"CBS_pEZY3:2289-2309:+", "CBS_pEZY3:2380-2394:-"}
        observed = self.rm_race._get_mask_base_indices(
            align_seg=self.test_align_seg_r2_positive, associated_primers=associated_primers)

        self.assertEqual(0, len(expected - observed))

    def test_get_mask_base_indices_race_r1_start(self):
        """Test that we return no indices for a RACE-like R1 starting at a primer start."""

        associated_primers = {"CBS_pEZY3:2380-2394:-"}
        observed = self.rm_race._get_mask_base_indices(
            align_seg=self.test_align_seg_r1_reverse, associated_primers=associated_primers)

        self.assertEqual(0, len(observed))

    def test_get_mask_base_indices_race_r2_end(self):
        """Test that we return no indices for a RACE-like R2 ending at a primer start."""

        associated_primers = {"CBS_pEZY3:2380-2394:-"}
        observed = self.rm_race._get_mask_base_indices(
            align_seg=self.test_align_seg_r2_positive, associated_primers=associated_primers)

        self.assertEqual(0, len(observed))

    def test_get_mask_base_indices_r1_contained(self):
        """Tests that a R1 completely contained in a primer is fully masked."""

        associated_primers = {"CBS_pEZY3:2289-2395:+"}

        observed = self.rm_tileseq._get_mask_base_indices(
            align_seg=self.test_align_seg_r1_reverse, associated_primers=associated_primers)

        self.assertEqual(len(observed), self.test_align_seg_r1_reverse.query_length)

    def test_get_mask_base_indices_r2_contained(self):
        """Tests that a R2 completely contained in a primer is fully masked."""

        associated_primers = {"CBS_pEZY3:2288-2394:-"}

        observed = self.rm_tileseq._get_mask_base_indices(
            align_seg=self.test_align_seg_r2_positive, associated_primers=associated_primers)

        self.assertEqual(len(observed), self.test_align_seg_r2_positive.query_length)

    def test_get_mask_base_indices_r1_flush(self):
        """Tests that a R1 flush with a primer at both ends is fully masked."""

        associated_primers = {"CBS_pEZY3:2289-2394:+"}

        observed = self.rm_tileseq._get_mask_base_indices(
            align_seg=self.test_align_seg_r1_reverse, associated_primers=associated_primers)

        self.assertEqual(len(observed), self.test_align_seg_r1_reverse.query_length)

    def test_get_mask_base_indices_r2_flush(self):
        """Tests that a R2 flush with a primer at both ends is fully masked."""

        associated_primers = {"CBS_pEZY3:2289-2394:-"}

        observed = self.rm_tileseq._get_mask_base_indices(
            align_seg=self.test_align_seg_r2_positive, associated_primers=associated_primers)

        self.assertEqual(len(observed), self.test_align_seg_r2_positive.query_length)

    def test_workflow(self):
        """Tests that the appropriate reads are masked."""

        # 10000001_R1 and 10015877_R1 duplicates should not have been masked with the true set of primers
        expected = {"10000004_R1", "10000004_R2", "10015877_R2"}

        # The proprocess BAM qname format is not compatible with masking
        # Remove the UMIs from the qnames so that it is compatible and does not raise a RuntimeError in _mask_reads()
        with tempfile.NamedTemporaryFile(mode="wb", suffix=".test.bam", delete=False, dir=self.tempdir) as test_bam, \
                pysam.AlignmentFile(self.preproc_bam, mode="rb") as in_af, \
                pysam.AlignmentFile(test_bam, mode="wb", header=in_af.header) as out_af:

            for align_seg in in_af.fetch(until_eof=True):
                align_seg.query_name = align_seg.query_name.split("_")[0]
                out_af.write(align_seg)

            test_bam_fn = test_bam.name

        rm_normal = rp.ReadMasker(in_bam=test_bam_fn, feature_file=self.primer_bed2, race_like=True, outdir=self.tempdir)
        rm_normal.workflow()

        observed = set()
        with pysam.AlignmentFile(rm_normal.out_bam, "rb") as test_af:
            for align_seg in test_af.fetch(until_eof=True):
                if 0 in set(align_seg.query_qualities):
                    observed.add(align_seg.get_tag(rp.UMITOOLS_UG_TAG))

        self.assertEqual(0, len(expected - observed))

# Skip testing VariantCallerPreprocessor as we simply make samtools calls which are tested in test_seq_utils
