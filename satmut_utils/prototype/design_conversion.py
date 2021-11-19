#!/usr/bin/env/python
"""Objects for converting reads to other saturation mutagenesis designs."""

import collections
import logging
import pysam
import random
import tempfile

from analysis.seq_utils import extract_seq, bam_to_fastq, sam_view, reverse_complement, DEFAULT_MIN_BQ, DEFAULT_MAX_BQ, \
    DNA_BASES, SAM_FLAG_SUPPL, SAM_FLAG_SECONDARY, SAM_FLAG_UNMAP, SAM_FLAG_MUNMAP
import core_utils.file_utils as fu
from core_utils.string_utils import make_random_str
from definitions import *

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

tempfile.tempdir = os.getenv("SCRATCH", "/tmp")
_logger = logging.getLogger(__name__)


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

        # First discard supplementary, secondary, unmapped alignments
        preprocessed_bam = sam_view(
            self.in_bam, None, "BAM", self.nthreads, F=SAM_FLAG_SUPPL + SAM_FLAG_SECONDARY + SAM_FLAG_UNMAP + SAM_FLAG_MUNMAP)

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


class ConvertToDmstools(object):
    """Converts reads to dms_tools input format."""

    DEFAULT_UMI_LEN = 8
    DEFAULT_OUTDIR = "."
    DEFAULT_NTHREADS = 0
    OUTPUT_PREFIX = "dms_tools"
    OUTBAM_SUFFIX = "dms_tools.bam"

    def __init__(self, in_bam, ref, pos_range, umi_len, outdir=DEFAULT_OUTDIR, nthreads=DEFAULT_NTHREADS):
        """Constructor for ConvertToDmstools.

        :param str in_bam: input alignments for a single tile
        :param str ref: reference FASTA
        :param str pos_range: 1-based positions flush with codons spanning the target
        :param int umi_len: length of randomer UMI. Default 8.
        :param str outdir: output directory. Default current directory.
        """

        self.in_bam = in_bam
        self.ref = ref
        self.pos_range = tuple(map(int, ",".split(pos_range)))
        self.umi_len = umi_len
        self.outdir = outdir
        self.nthreads = nthreads
        self.out_prefix = os.path.join(self.outdir, fu.replace_extension(os.path.basename(in_bam), self.OUTPUT_PREFIX))

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        self.out_bam = os.path.join(self.outdir, fu.replace_extension(
            os.path.basename(self.in_bam), self.OUTBAM_SUFFIX))

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
        quals_res = add_quals + quals
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
        quals_res = quals + add_quals
        return seq_res, quals_res

    @staticmethod
    def _update_read(align_seg, seq, quals):
        """Updates read sequence and qualities by side effect.

        :param pysam.AlignedSegment align_seg: read object
        :param list seq: sequence to reassign
        :param list quals: quals to reassign
        """

        # Sequence must always be updated before qualities
        align_seg.query_alignment_sequence = "".join(seq)
        align_seg.query_alignment_qualities = quals

    def _convert_reads(self, paired_bam):
        """Makes reads flush with the reference and adds UMIs to the 5' end of the reads.

        :param str paired_bam: BAM with only paired reads
        """

        with pysam.AlignmentFile(paired_bam, mode="rb") as in_af, \
                pysam.AlignmentFile(self.out_bam, mode="wb", header=in_af.header) as out_af:

            used_umis = set()

            for align_seg in in_af.fetch(until_eof=True):

                ref_positions = align_seg.get_reference_positions()
                min_ref_pos = min(ref_positions) + 1
                max_ref_pos = max(ref_positions) + 1
                start_diff = self.pos_range[0] - min_ref_pos
                end_diff = self.pos_range[1] - max_ref_pos

                if start_diff > 0:
                    seq, quals = self._trim_read_start(
                        seq=align_seg.query_alignment_sequence, quals=align_seg.query_alignment_qualities,
                        idx=start_diff)
                    self._update_read(align_seg, seq, quals)

                if end_diff < 0:
                    seq, quals = self._trim_read_end(
                        seq=align_seg.query_alignment_sequence, quals=align_seg.query_alignment_qualities,
                        idx=align_seg.query_alignment_length - abs(end_diff))
                    self._update_read(align_seg, seq, quals)

                if start_diff < 0:
                    seq, quals = self._add_read_start(
                        contig=align_seg.reference_name, ref_pos=min_ref_pos,
                        seq=align_seg.query_alignment_sequence, quals=align_seg.query_alignment_qualities,
                        nbases=abs(start_diff))
                    self._update_read(align_seg, seq, quals)

                if end_diff > 0:
                    seq, quals = self._add_read_end(
                        contig=align_seg.reference_name, ref_pos=max_ref_pos,
                        seq=align_seg.query_alignment_sequence, quals=align_seg.query_alignment_qualities,
                        nbases=end_diff)
                    self._update_read(align_seg, seq, quals)

                # Finally append UMIs/barcodes to each read
                self._append_umi(align_seg, used_umis)

                # Now that reads have been modified, write paired reads to FASTQ and re-align
                align_seg.cigarstring = None
                out_af.write(align_seg)

    def _append_umi(self, align_seg, used_umis):
        """Appends a UMI to the 5' end of the read.

        :param pysam.AlignedSegment align_seg: read object
        :param set used_umis: set of UMIs already assigned
        """

        umi = make_random_str(str_len=self.umi_len, letters=DNA_BASES)
        while umi not in used_umis:
            umi = make_random_str(str_len=self.umi_len, letters=DNA_BASES)

        used_umis.add(umi)

        umi_quals = [random.randint(DEFAULT_MIN_BQ, DEFAULT_MAX_BQ) for _ in range(self.umi_len)]

        if not align_seg.is_reverse:
            align_seg.query_alignment_sequence = umi + align_seg.query_alignment_sequence
            align_seg.query_alignment_qualities = umi_quals + list(align_seg.query_alignment_qualities)
        else:
            umi = reverse_complement(umi)
            align_seg.query_alignment_sequence += umi
            align_seg.query_alignment_qualities = list(align_seg.query_alignment_qualities) + umi_quals

    def _write_fastqs(self):
        """Writes and gzips FASTQs.

        :return tuple: (str, str) paths of the R1 and R2 FASTQ files
        """

        r1_fastq, r2_fastq = bam_to_fastq(bam=self.out_bam, out_prefix=self.out_prefix)
        zipped_r1_fastq = fu.gzip_file(r1_fastq, force=True)
        zipped_r2_fastq = fu.gzip_file(r2_fastq, force=True)
        return zipped_r1_fastq, zipped_r2_fastq

    def workflow(self):
        """Runs the dms_tools conversion workflow.

        :return tuple: (str, str) paths of the R1 and R2 FASTQ files
        """

        _logger.info("Started dms_tools conversion workflow.")

        # We must first determine if there are any singleton reads or reads with secondary and supplementary alignments
        # because we have to write paired FASTQs.
        # Do not rely on SAM flags because the alignments may have been previously intersected, which may lead
        # to loss of mates (and SAM flags are not typically updated following manipulation).
        ffp = FilterForPairs(in_bam=self.in_bam, outdir=self.outdir, nthreads=self.nthreads)

        # Now iterate over the reads, make them flush with the range positions, add UMIs, and exclude non-paired reads
        self._convert_reads(ffp.out_bam)

        # Write the FASTQs
        zipped_r1_fastq, zipped_r2_fastq = self._write_fastqs()

        _logger.info("Completed dms_tools conversion workflow.")
        return zipped_r1_fastq, zipped_r2_fastq
