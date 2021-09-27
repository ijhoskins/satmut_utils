#!/usr/bin/env/python
"""Collection of project definitions."""

import os

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
PROJECT_LAB = "Can_Cenik_lab"
PROJECT_AUTHOR = "Ian_Hoskins"

APPRIS_CONTIG_IDS = "appris_human_v1_actual_regions_contigs.txt"
APPRIS_TRX_IDS = "appris_human_v1_actual_regions_trx_ids.txt"
APPRIS_GENE_IDS = "appris_human_v1_actual_regions_gene_ids.txt"
APPRIS_TRX_FASTA = "appris_human_24_01_2019_selected.fa"
APPRIS_TRX_GFF = "gencode.v29.annotation.gtf"
GRCH38_FASTA = "GRCh38.fa"

# CTCCTTTC
NEXTERA_R1_R2_ADAPTER_TRIM_1 = "CTGTCTCTTATACACATCT"
NEXTERA_R1_R2_ADAPTER_TRIM_2 = "AGATGTGTATAAGAGACAG"


NEB_ADAPTER_P5 = "AGACGTGTGCTCTTCCGATCT"  # 5' tail of tile-seq R primers
NEB_ADAPTER_P5_RC = "AGATCGGAAGAGCACACGTCT"  # to be trimmed from 3' end of R1 in tile-seq with NEBNext adapters
NEB_ADAPTER_P7 = "TACACGACGCTCTTCCGATCT"  # 5' tail of tile-seq F primers
NEB_ADAPTER_P7_RC = "AGATCGGAAGAGCGTCGTGTA"   # to be trimmed from 3' end of R2 in tile-seq with NEBNext adapters

PEZY3_ATTB2_P5 = "ACCACTTTGTACAAGAAAGTTGG"  # exists only on the terminal R primer (for last tile)
PEZY3_ATTB2_P5_RC = "CCAACTTTCTTGTACAAAGTGGT"
PEZY3_ATTB1_P7 = "CAAGTTTGTACAAAAAAGTTGGC"  # exists only on the terminal F primer (for tile 1)
PEZY3_ATTB1_P7_RC = "GCCAACTTTTTTGTACAAACTTG"

PEZY3_ATTB1_EXT_P7 = "TGCTGGAATTGATTAATA"  # 5' Extension for tile_1_F primer

TILESEQ_UMI_LEN = 8
TILESEQ_NM_ALLOWANCE = 2

TILESEQ_UMI_REGEX = "(?P<umi_1>[ATCG]{%i})" % TILESEQ_UMI_LEN + \
                    "(?P<discard_1>%s){e<=%i}" % (PEZY3_ATTB1_EXT_P7, TILESEQ_NM_ALLOWANCE)

# CR=12 bp adapter duplex "common region" plus T from A tail ligation (13 bp)
AMP_CR = "AACCGCCAGGAGT"
AMP_CR_RC = "ACTCCTGGCGGTT"
AMP_GSP2_TAIL = "GTTCAGACGTGTGCTCTTCCGATCT"
AMP_GSP2_TAIL_RC = "AGATCGGAAGAGCACACGTCTGAAC"

AMP_UMI_LEN = 8
AMP_UMI_BUFFER = 6
AMP_CR_NM_ALLOWANCE = 2
AMP_UMI_NM_ALLOWANCE = 1

# The final UMI will be composed of the 8 bp MBC, and the first bases of the start site determined by the buffer
# This uses the python regex library syntax with named capture groups and fuzzy matching
AMP_UMI_REGEX = "(?P<umi_1>[ATCG]{%i})" % AMP_UMI_LEN + \
                "(?P<discard_1>%s){e<=%i}" % (AMP_CR, AMP_CR_NM_ALLOWANCE)
