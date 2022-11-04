#!/usr/bin/env python3
"""Objects for converting reads to other saturation mutagenesis input or output read requirements."""

import collections
import logging
import os
import pysam
import pybedtools
import random
import statistics
import tempfile
import warnings

from analysis.seq_utils import extract_seq, sort_bam, bam_to_fastq, sam_view, reverse_complement, \
    DEFAULT_MIN_BQ, DEFAULT_MAX_BQ, DNA_BASES, SAM_FLAG_SUPPL, SAM_FLAG_SECONDARY, SAM_FLAG_UNMAP, \
    SAM_FLAG_MUNMAP, SAM_CIGAR_INS, SAM_CIGAR_DEL
import core_utils.file_utils as fu
from core_utils.string_utils import make_random_str
from core_utils.vcf_utils import VCF_CAO_ID, VCF_CONTIG_INDEX, VCF_POS_INDEX, VCF_REF_INDEX, VCF_ALT_INDEX, VCF_DP_ID
from satmut_utils.definitions import *
from scripts.run_bowtie2_aligner import workflow as align_workflow

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

tempfile.tempdir = os.getenv("SCRATCH", "/tmp")
logger = logging.getLogger(__name__)


class FilterForPairs(object):
    """Filters BAM files to contain only primary, paired alignments."""

    DEFAULT_OUTDIR = "."
    DEFAULT_NTHREADS = 0
    DEFAULT_SUFFIX = "pairs.bam"

    def __init__(self, in_bam, outdir=DEFAULT_OUTDIR, nthreads=DEFAULT_NTHREADS):
        """Constructor for FilterForPairs.

        :param str in_bam: input alignments for a single tile
        :param str outdir: output directory. Default current directory.
        :param int nthreads: Number of threads to use for BAM operations. Default 0 (autodetect).
        """

        self.in_bam = in_bam
        self.outdir = outdir
        self.nthreads = nthreads
        self.out_bam = os.path.join(self.outdir, fu.replace_extension(self.in_bam, self.DEFAULT_SUFFIX))
        self.workflow()

    @staticmethod
    def _get_nonpaired_qnames(in_bam):
        """Determines the non-paired qnames in the input BAM.

        :param str in_bam: input BAM
        :return set: non-paired qnames to filter from the BAM to generate paired reads only.
        """

        with pysam.AlignmentFile(in_bam, mode="rb") as in_af:
            qname_list = [align_seg.query_name for align_seg in in_af.fetch(until_eof=True)]
            in_af.reset()
            qname_counter = collections.Counter(qname_list)
            nonpaired_qnames = {qname for qname, count in qname_counter.items() if count != 2}
            return nonpaired_qnames

    def workflow(self):
        """Runs the pair-filtering workflow."""

        logger.info("Filtering out unmapped, supplementary, secondary, and any non-paired reads.")

        # First discard supplementary, secondary, unmapped alignments
        preprocessed_bam = sam_view(
            self.in_bam, None, "BAM", self.nthreads,
            F=SAM_FLAG_SUPPL + SAM_FLAG_SECONDARY + SAM_FLAG_UNMAP + SAM_FLAG_MUNMAP)

        # This may still not be enough preprocessing if the user intersected the reads and not updated the SAM flags
        # For instance, if one read of a pair does not intersect, it will not be in the input but still contain a
        # paired bit in the flag. Thus determine the non-paired qnames and filter them manually
        nonpaired_qnames = self._get_nonpaired_qnames(preprocessed_bam)

        with pysam.AlignmentFile(preprocessed_bam, mode="rb") as in_af, \
                pysam.AlignmentFile(self.out_bam, mode="wb", header=in_af.header) as out_af:

            for align_seg in in_af.fetch(until_eof=True):

                if align_seg.query_name in nonpaired_qnames:
                    continue

                out_af.write(align_seg)

        fu.safe_remove((preprocessed_bam,))
        return self.out_bam


class DesignConverter(object):
    """Converts reads to dms_tools or Enrich2 required input format."""

    DEFAULT_UMI_LEN = 8
    DEFAULT_NO_UMI = False
    DEFAULT_OUTDIR = "."
    DEFAULT_NTHREADS = 0
    POS_RANGE_DELIM = ","
    OUTPUT_PREFIX = "design_convert"
    OUTBAM_SUFFIX = "design_convert.bam"

    def __init__(self, in_bam, ref, pos_range, umi_len=DEFAULT_UMI_LEN, no_umis=DEFAULT_NO_UMI, outdir=DEFAULT_OUTDIR,
                 nthreads=DEFAULT_NTHREADS):
        """Constructor for DesignConverter.

        :param str in_bam: input alignments for a single tile
        :param str ref: reference FASTA
        :param str pos_range: 1-based positions flush with codons spanning the target
        :param int umi_len: length of randomer UMI/barcode. Default 8.
        :param bool no_umis: do not add UMIs, only make reads flush with the pos_range. Default False.
        :param str outdir: output directory. Default current directory.
        :param int nthreads: Number of threads to use for BAM operations. Default 0 (autodetect).
        """

        self.in_bam = in_bam
        self.ref = ref
        self.pos_range = tuple(map(int, pos_range.split(self.POS_RANGE_DELIM)))
        self.umi_len = umi_len
        self.no_umis = no_umis
        self.outdir = outdir
        self.nthreads = nthreads
        self.out_prefix = os.path.join(self.outdir, fu.replace_extension(os.path.basename(in_bam), self.OUTPUT_PREFIX))

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        self.temp_flush_bam = tempfile.NamedTemporaryFile(suffix=".flush.bam", mode="wb", delete=False).name
        self.temp_sort_bam = tempfile.NamedTemporaryFile(suffix=".qsort.bam", mode="wb", delete=False).name

    @staticmethod
    def _trim_read_start(seq, quals, idx):
        """Trims read sequence and qualities at the start.

        :param list seq: sequence
        :param list quals: base qualities
        :param int idx: 0-based start index
        :return tuple: (list, list) sequence and quals trimmed at the start
        """

        seq_res = seq[idx:]
        quals_res = quals[idx:]
        return seq_res, quals_res

    @staticmethod
    def _trim_read_end(seq, quals, idx):
        """Trims read sequence and qualities at the end.

        :param list seq: sequence
        :param list quals: base qualities
        :param int idx: 0-based end index
        :return tuple: (list, list) sequence and quals trimmed at the end
        """

        seq_res = seq[:idx + 1]
        quals_res = quals[:idx + 1]
        return seq_res, quals_res

    def _add_read_start(self, contig, ref_pos, seq, quals, nbases):
        """Adds read sequence and qualities at the start.

        :param str contig: reference name
        :param int ref_pos: reference position of the read start, 1-based
        :param list seq: sequence
        :param list quals: base qualities
        :param int nbases: number bases to add
        :return tuple: (list, list) sequence and quals with buffered sequence
        """

        add_seq = extract_seq(contig=contig, start=ref_pos - nbases, stop=ref_pos - 1, ref=self.ref)
        add_quals = [random.randint(DEFAULT_MIN_BQ, DEFAULT_MAX_BQ) for _ in range(nbases)]
        seq_res = add_seq + seq
        quals_res = add_quals + list(quals)
        return seq_res, quals_res

    def _add_read_end(self, contig, ref_pos, seq, quals, nbases):
        """Adds read sequence and qualities at the end.

        :param str contig: reference name
        :param int ref_pos: reference position of the read end, 1-based
        :param list seq: sequence
        :param list quals: base qualities
        :param int nbases: number bases to add
        :return tuple: (list, list) sequence and quals with buffered sequence
        """

        add_seq = extract_seq(contig=contig, start=ref_pos + 1, stop=ref_pos + nbases, ref=self.ref)
        add_quals = [random.randint(DEFAULT_MIN_BQ, DEFAULT_MAX_BQ) for _ in range(nbases)]
        seq_res = seq + add_seq
        quals_res = list(quals) + add_quals
        return seq_res, quals_res

    @staticmethod
    def _update_read(align_seg, seq, quals):
        """Updates read sequence and qualities by side effect.

        :param pysam.AlignedSegment align_seg: read object
        :param list seq: sequence to reassign
        :param list quals: quals to reassign
        """

        # Sequence must always be updated before qualities
        align_seg.query_sequence = "".join(seq)
        align_seg.query_qualities = quals

    def _convert_reads(self, paired_bam):
        """Makes reads flush with the reference and adds UMIs to the 5' end of the reads.

        :param str paired_bam: BAM with only paired reads
        """

        with pysam.AlignmentFile(paired_bam, mode="rb") as in_af, \
                pysam.AlignmentFile(self.temp_flush_bam, mode="wb", header=in_af.header) as out_af:

            used_umis = set()
            filt_reads = set()

            for align_seg in in_af.fetch(until_eof=True):

                # Filter any reads that have InDels, as all reads must be the same length as the WT reference
                cigar_set = set(align_seg.cigarstring)
                if SAM_CIGAR_INS in cigar_set or SAM_CIGAR_DEL in cigar_set:
                    filt_reads.add(align_seg.query_name)
                    continue

                # For any mate of a read that has InDel, also filter it because we must have only pairs in the output
                if align_seg.query_name in filt_reads:
                    continue

                ref_positions = align_seg.get_reference_positions()
                min_ref_pos = min(ref_positions) + 1
                max_ref_pos = max(ref_positions) + 1
                start_diff = self.pos_range[0] - min_ref_pos
                end_diff = self.pos_range[1] - max_ref_pos
                seq = align_seg.query_alignment_sequence
                quals = list(align_seg.query_alignment_qualities)

                if end_diff < 0:
                    seq, quals = self._trim_read_end(
                        seq=seq, quals=quals, idx=align_seg.query_alignment_length - abs(end_diff))

                if start_diff > 0:
                    seq, quals = self._trim_read_start(seq=seq, quals=quals, idx=start_diff)

                if start_diff < 0:
                    seq, quals = self._add_read_start(
                        contig=align_seg.reference_name, ref_pos=min_ref_pos,
                        seq=seq, quals=quals, nbases=abs(start_diff))

                if end_diff > 0:
                    seq, quals = self._add_read_end(
                        contig=align_seg.reference_name, ref_pos=max_ref_pos, seq=seq, quals=quals, nbases=end_diff)

                # Finally append UMIs/barcodes to each read
                if not self.no_umis:
                    seq, quals = self._append_umi(align_seg, used_umis, seq, quals)

                self._update_read(align_seg, seq, quals)

                # Now that reads have been modified, write paired reads to FASTQ and re-align
                align_seg.cigarstring = None
                out_af.write(align_seg)

    def _append_umi(self, align_seg, used_umis, seq, quals):
        """Appends a UMI to the 5' end of the read.

        :param pysam.AlignedSegment align_seg: read object
        :param set used_umis: set of UMIs already assigned
        :param list seq: codon-flush sequence
        :param list quals: codon-flush  base qualities
        :return tuple: sequence and quals lists
        """

        umi = make_random_str(str_len=self.umi_len, letters=DNA_BASES)
        while umi in used_umis:
            umi = make_random_str(str_len=self.umi_len, letters=DNA_BASES)

        used_umis.add(umi)
        umi_quals = [random.randint(DEFAULT_MIN_BQ, DEFAULT_MAX_BQ) for _ in range(self.umi_len)]

        if not align_seg.is_reverse:
            seq = umi + seq
            quals = umi_quals + quals
        else:
            umi = reverse_complement(umi)
            seq += umi
            quals += umi_quals

        return seq, quals

    def _write_fastqs(self):
        """Writes and gzips FASTQs.

        :return tuple: (str, str) paths of the R1 and R2 FASTQ files
        """

        # Discard any singletons
        sort_bam(bam=self.temp_flush_bam, output_am=self.temp_sort_bam,
                 output_format="BAM", by_qname=True, nthreads=self.nthreads)
        r1_fastq, r2_fastq = bam_to_fastq(self.temp_sort_bam, self.out_prefix, True, self.nthreads, s=os.devnull)
        zipped_r1_fastq = fu.gzip_file(r1_fastq, force=True)
        zipped_r2_fastq = fu.gzip_file(r2_fastq, force=True)
        fu.safe_remove((self.temp_sort_bam,))
        return zipped_r1_fastq, zipped_r2_fastq

    def workflow(self):
        """Runs the dms_tools2/Enrich2 design conversion workflow.

        :return tuple: (str, str) paths of the R1 and R2 FASTQ files
        """

        logger.info("Started design conversion workflow.")

        # We must first determine if there are any singleton reads or reads with secondary and supplementary alignments
        # because we have to write strictly paired FASTQs. Filter them out using SAM flags.
        # Note however that if in_bam consists of reads intersected with a target tile, we have to manually remove
        # singletons.
        ffp = FilterForPairs(in_bam=self.in_bam, outdir=self.outdir, nthreads=self.nthreads)

        # Now iterate over the reads, make them flush with the range positions, exclude reads/pairs with InDels,
        # then add UMIs/barcodes (for dms_tools2 conversion)
        logger.info("Converting reads: filtering pairs with InDels, and trimming/appending sequence to make reads "
                    "flush with the pos_range.")
        self._convert_reads(ffp.out_bam)

        # Data now in temp BAM, clean up the pairs.bam
        fu.safe_remove((ffp.out_bam,))

        logger.info("Writing FASTQs for converted reads.")
        zipped_r1_fastq, zipped_r2_fastq = self._write_fastqs()

        logger.info("Globally re-aligning converted reads.")
        nthreads = self.nthreads if self.nthreads != 0 else 1
        align_workflow(
            f1=zipped_r1_fastq, f2=zipped_r2_fastq, ref=self.ref, outdir=self.outdir, local=False, nthreads=nthreads)

        logger.info("Completed design conversion workflow.")
        return zipped_r1_fastq, zipped_r2_fastq


class DmsTools2ToSatmutUtils(object):
    """Class for converting reads produced by barcoded subamplicon sequencing to satmut_utils compatible input."""

    DEFAULT_OUTDIR = "."
    DEFAULT_UMI_LEN = 8
    DEFAULT_EXT = "satmut.fastq"

    def __init__(self, r1_fastq, r2_fastq, umi_len=DEFAULT_UMI_LEN, outdir=DEFAULT_OUTDIR):
        """Constructor for DmsTools2ToSatmutUtils.

        :param str r1_fastq: R1 FASTQ filepath.
        :param str r2_fastq: R2 FASTQ filepath.
        :param int umi_length: length of UMIs at start of each read. Default 8.
        :param str outdir: Optional output directory. Default current directory.

        Notes
        -----
        This class is meant to accept dms_tools2 reads that have had adapters trimmed. R1 and R2 should still contain \
        UMIs/barcodes. The R2 barcode will be transferred to the start of R1 to enable satmut_utils analysis with the \
        --umi_regex flag.
        """

        self.r1_fastq = r1_fastq
        self.r2_fastq = r2_fastq
        self.umi_length = umi_len
        self.outdir = outdir

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        self.output_r1_fastq = os.path.join(
            outdir, fu.replace_extension(os.path.basename(self.r1_fastq), self.DEFAULT_EXT))

        self.output_r2_fastq = os.path.join(
            outdir, fu.replace_extension(os.path.basename(self.r2_fastq), self.DEFAULT_EXT))

    def workflow(self):
        """Runs the conversion workflow- move R2 UMIs to start of R1.

        :return tuple: R1 and R2 gzipped FASTQs
        """

        logger.info("Started dms_tools2 to satmut_utils FASTQ conversion workflow.")

        with pysam.FastxFile(self.r1_fastq) as r1_input_fastq,\
                pysam.FastxFile(self.r2_fastq) as r2_input_fastq,\
                open(self.output_r1_fastq, "w") as r1_output_fastq,\
                open(self.output_r2_fastq, "w") as r2_output_fastq:

            for r1, r2 in zip(r1_input_fastq, r2_input_fastq):

                r2_umi = r2.sequence[:self.umi_length]
                r2_umi_quals = r2.quality[:self.umi_length]

                # Add the R2 UMI and qualities to R1
                r1.sequence = r2_umi + r1.sequence
                r1.quality = r2_umi_quals + r1.quality

                # Remove the UMI from R2
                r2.sequence = r2.sequence[self.umi_length:]
                r2.quality = r2.quality[self.umi_length:]

                r1_output_fastq.write(str(r1) + fu.FILE_NEWLINE)
                r2_output_fastq.write(str(r2) + fu.FILE_NEWLINE)

        # Gzip the FASTQs
        logger.info("Compressing converted FASTQ files.")
        zipped_r1_fastq = fu.gzip_file(self.output_r1_fastq, force=True)
        zipped_r2_fastq = fu.gzip_file(self.output_r2_fastq, force=True)

        logger.info("Completed dms_tools2 to satmut_utils FASTQ conversion workflow.")
        return zipped_r1_fastq, zipped_r2_fastq


class SatmutUtilsToDiMSum(object):
    """Class for converting VCF variants into full coding sequence strings."""

    DEFAULT_OUTDIR = "."
    DEFAULT_EXT = "DiMSum_counts.txt"
    DEFAULT_HEADER = ("VAR_ID", "nt_seq", "counts")
    VAR_ID_DELIM = ":"

    def __init__(self, vcf_summary, reference, cds_bed, outdir=DEFAULT_OUTDIR):
        """
        :param str vcf_summary: satmut_utils VCF summary.txt file to process
        :param str reference: reference FASTA
        :param str cds_bed: Single-line BED file with CDS annotation relative to the reference FASTA
        :param str outdir: Optional output directory. Default current directory.
        """

        self.vcf_summary = vcf_summary
        self.reference = reference
        self.cds_bed = cds_bed
        self.outdir = outdir
        self.output_file = os.path.join(
            outdir, fu.replace_extension(os.path.basename(vcf_summary), self.DEFAULT_EXT))
        self.contig_name, self.cds_seq, self.cds_start = self._get_cds_list()
        self.cds_len = len(self.cds_seq)

    def _get_cds_list(self):
        """Extracts the contig name and CDS sequence given the CDS BED and reference FASTA.

        :return tuple: (str, list, int) contig name, CDS sequence, 0-based CDS start relative to reference
        """

        # Extract the CDS sequence
        cds_bedtool = pybedtools.BedTool(self.cds_bed)
        contig = cds_bedtool[0].chrom
        cds_start = cds_bedtool[0].start
        extract_fa = cds_bedtool.sequence(fi=self.reference)

        with pysam.FastaFile(extract_fa.seqfn) as cds_fasta:
            cds_seq_list = list(cds_fasta.fetch(reference=cds_fasta.references[0]).upper())

        # Also return start position of the CDS so we can get the offset from the whole reference
        # start coordinate is always 0-based
        return contig, cds_seq_list, cds_start

    def vcf_to_string(self, pos, ref, alt):
        """Constructs a coding sequence string from a VCF variant call.

        :param int pos: 1-based variant position
        :param str ref: reference nucleotides
        :param str alt: alternate nucleotides
        :return str | None: coding sequence string
        """

        cds_pos = pos - self.cds_start - 1

        if cds_pos < 0 or pos > (self.cds_len + self.cds_start):
            warnings.warn("Variant %i:%s:%s is not coding. Skipping." % (pos, ref, alt))
            return None

        ref_len = len(ref)
        exp_ref = "".join(self.cds_seq[cds_pos:cds_pos + ref_len])

        if ref != exp_ref:
            warnings.warn("Variant REF field %s at POS %i does not match the reference sequence %s" % (ref, pos, exp_ref))

        cds_seq = "".join(self.cds_seq[:cds_pos] + list(alt) + self.cds_seq[cds_pos + ref_len:])

        return cds_seq

    @staticmethod
    def _extract_wt_counts(pos, dp, alt, wt_dict):
        """Collects WT counts from each position.

        :param int pos: variant position
        :param int dp: total depth of coverage at position
        :param int alt: variant count at position
        :param collections.OrderedDict wt_dict: dictionary holding DP and sum of variant counts
        """

        # Only collect DP for each position uniquely; satmut_utils output ensures DP for multiple variants at the
        # same position is the same
        if pos not in wt_dict:
            wt_dict[pos] = [0, 0]
            wt_dict[pos][0] += dp

        # Add the variant counts (possibly multiple per position)
        wt_dict[pos][0] += alt

    @staticmethod
    def _get_total_wt_count(wt_dict):
        """Aggregates per-position WT counts to generate a total WT count for the whole target region.

        :param collections.OrderedDict wt_dict: dictionary holding DP and sum of variant counts
        :return int: total WT count
        """

        wt_counts = []
        for dp, var_count in wt_dict.values():
            wt_counts.append(dp - var_count)

        total_wt_count = int(statistics.median(wt_counts))

        return total_wt_count

    def workflow(self):
        """Makes a three-column file containing variant ID, variant string, and counts."""

        with open(self.vcf_summary, "r") as vcf_summ_file, open(self.output_file, "w") as out_file:

            var_ids = set()
            wt_dict = collections.OrderedDict()

            for i, line in enumerate(vcf_summ_file):

                line_split = line.rstrip().split(fu.FILE_DELIM)

                if i == 0:

                    try:
                        dp_index = line_split.index(VCF_DP_ID)
                        cao_index = line_split.index(VCF_CAO_ID)
                    except ValueError:
                        raise RuntimeError("Input file does not have %s or %s fields." % (VCF_DP_ID, VCF_CAO_ID))

                    out_file.write(fu.FILE_DELIM.join(self.DEFAULT_HEADER) + fu.FILE_NEWLINE)
                    continue

                var_id = self.VAR_ID_DELIM.join(
                    line_split[0:2] + [line_split[VCF_REF_INDEX]] + [line_split[VCF_ALT_INDEX]])

                # Use because satmut_utils vcf.summary.txt files contain multiple records for each base in MNPs,
                # which contain the same VAR_ID and count/frequency
                if var_id in var_ids:
                    continue

                var_ids.add(var_id)
                var_pos = int(line_split[VCF_POS_INDEX])
                dp_count = line_split[dp_index]
                alt_count = line_split[cao_index]

                self._extract_wt_counts(pos=var_pos, dp=int(dp_count), alt=int(alt_count), wt_dict=wt_dict)

                nt_string = self.vcf_to_string(
                    pos=var_pos, ref=line_split[VCF_REF_INDEX], alt=line_split[VCF_ALT_INDEX])

                if nt_string is None:
                    continue

                out_file.write(fu.FILE_DELIM.join((var_id, nt_string, alt_count,)) + fu.FILE_NEWLINE)

            # Add the WT counts by taking the median of the DP-CAO across each targeted position
            # DiMSum requires a single count for the WT coverage
            wt_id = self.VAR_ID_DELIM.join((self.contig_name, str(self.cds_start),))

            # WARNING: assumes there are no mismatches to the reference to accurately report WT sequence
            wt_string = "".join(self.cds_seq)
            total_wt_count = self._get_total_wt_count(wt_dict=wt_dict)

            out_file.write(fu.FILE_DELIM.join((wt_id, wt_string, str(total_wt_count),)) + fu.FILE_NEWLINE)
