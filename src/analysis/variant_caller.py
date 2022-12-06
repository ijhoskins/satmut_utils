#!/usr/bin/env python3
"""Variant caller for SNPs and MNPs."""

import collections
import copy
import itertools
import logging
import numpy as np
import os
import pybedtools
import pysam
import statistics
import tempfile

import analysis.coordinate_mapper as cm
import analysis.read_preprocessor as rp
from analysis.references import APPRIS_CONTIG_DELIM, APPRIS_TRX_INDEX
import analysis.seq_utils as su

from core_utils.feature_file_utils import intersect_features
import core_utils.file_utils as fu
import core_utils.vcf_utils as vu

from satmut_utils.definitions import PROJECT_ROOT, PROJECT_LAB, PROJECT_AUTHOR, DEFAULT_MUT_SIG, VALID_MUT_SIGS, DEFAULT_TEMPDIR

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

COORDINATE_KEY = collections.namedtuple("COORDINATE_KEY", "contig, pos")
MM_TUPLE = collections.namedtuple("MM_TUPLE", "contig, pos, ref, alt, bq, read_pos")
CALL_TUPLE = collections.namedtuple("CALL_TUPLE", "contig, pos, ref, alt, refs, alts, positions")
PER_BP_STATS = collections.namedtuple("PER_BP_STATS", "r1_bqs, r2_bqs, r1_read_pos, r2_read_pos")
VARIANT_CALL_KEY_TUPLE = collections.namedtuple("VARIANT_CALL_KEY_TUPLE", "contig, pos, ref, alt, index")

VARIANT_CALL_SUMMARY_TUPLE = collections.namedtuple(
    "VARIANT_CALL_SUMMARY_TUPLE",
    (vu.VCF_POS_NT_ID, vu.VCF_REF_NT_ID, vu.VCF_ALT_NT_ID, vu.VCF_UP_REF_NT_ID, vu.VCF_DOWN_REF_NT_ID,
     vu.VCF_DP_ID, vu.VCF_CAO_ID, vu.VCF_NORM_CAO_ID, vu.VCF_CAF_ID,
     vu.VCF_R1_PLUS_AO_ID, vu.VCF_R1_MINUS_AO_ID, vu.VCF_R2_PLUS_AO_ID, vu.VCF_R2_MINUS_AO_ID,
     vu.VCF_R1_PLUS_MED_POS_ID, vu.VCF_R1_MINUS_MED_POS_ID, vu.VCF_R2_PLUS_MED_POS_ID, vu.VCF_R2_MINUS_MED_POS_ID,
     vu.VCF_R1_PLUS_MED_BQ_ID, vu.VCF_R1_MINUS_MED_BQ_ID, vu.VCF_R2_PLUS_MED_BQ_ID, vu.VCF_R2_MINUS_MED_BQ_ID,
     vu.VCF_R1_PLUS_MED_NM_ID, vu.VCF_R1_MINUS_MED_NM_ID, vu.VCF_R2_PLUS_MED_NM_ID, vu.VCF_R2_MINUS_MED_NM_ID)
)

tempfile.tempdir = DEFAULT_TEMPDIR

logger = logging.getLogger(__name__)


class VariantCaller(object):
    """Class for calling variants across target regions."""

    VARIANT_CALL_REFERENCE_DIR = "/tmp/references"
    VARIANT_CALL_ENSEMBL_ID = None
    VARIANT_CALL_REF = None
    VARIANT_CALL_GFF = None
    VARIANT_CALL_GFF_REF = None
    VARIANT_CALL_TARGET = None
    VARIANT_CALL_PRIMERS = None
    VARIANT_CALL_DEDUP = rp.DEDUP_FLAG
    VARIANT_CALL_CDEDUP = rp.CDEDUP_FLAG
    VARIANT_CALL_OUTDIR = "."
    VARIANT_CALL_PREFIX = "./out"
    VARIANT_CALL_STATS = False

    VARIANT_CALL_MIN_BQ = 30
    VARIANT_CALL_MIN_DP = 2
    VARIANT_CALL_MAX_NM = 10
    VARIANT_CALL_NORM_DP = 1000000
    VARIANT_CALL_COV_EXT = "cov.bedgraph"
    VARIANT_CALL_REF_CANDIDATE_EXT = "var.cand.vcf"
    VARIANT_CALL_MAX_MNP_WINDOW = 3

    DEFAULT_NTHREADS = 0
    _STATS_DELIM = ","

    R1_PLUS_INDEX = 0
    R1_MINUS_INDEX = 1
    R2_PLUS_INDEX = 2
    R2_MINUS_INDEX = 3
    R_COUNTS_INDEX = 0
    R_BQ_INDEX = 1
    R_RP_INDEX = 2
    R_NM_INDEX = 3

    def __init__(self, am, ref, trx_gff, gff_ref, targets=VARIANT_CALL_TARGET, primers=VARIANT_CALL_PRIMERS,
                 output_dir=VARIANT_CALL_OUTDIR, nthreads=DEFAULT_NTHREADS, mut_sig=DEFAULT_MUT_SIG):
        r"""Constructor for VariantCaller.

        :param str am: SAM/BAM file to enumerate variants in
        :param str ref: path to reference FASTA used in alignment. Must be samtools faidx indexed.
        :param str trx_gff: GFF file containing transcript metafeatures and exon features, in 5' to 3' order, \
        regardless of strand. Ordering is essential.
        :param str gff_ref: reference FASTA corresponding to the GFF features
        :param str targets: BED or GFF containing target regions to call variants in
        :param str | None primers: BED or GFF file containing primers. For masking synthetic sequences for \
        accurate AF of variants overlapping primer regions. This feature file must contain the strand of the primer. \
        Set to None for no masking.
        :param str | None output_dir: output dir to use for output files; if None, will create a tempdir.
        :param int nthreads: number threads to use for SAM/BAM file manipulations. Default 0 (autodetect).
        :param str mut_sig: mutagenesis signature- one of {NNN, NNK, NNS}. Default NNK.
        :raises RuntimeError: if no alignments are found in the input BAM
        """

        logger.info("Initializing %s" % self.__class__.__name__)

        self.am = am
        self.ref = ref
        self.transcript_gff = trx_gff
        self.gff_reference = gff_ref
        self.targets = targets
        self.primers = primers
        self.output_dir = output_dir
        self.nthreads = nthreads
        self.mut_sig = mut_sig

        if output_dir is None:
            self.output_dir = tempfile.mkdtemp(suffix=__class__.__name__)

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

        # Preprocess the alignments and setup the output directory
        self.vc_preprocessor = rp.VariantCallerPreprocessor(am=am, ref=ref, output_dir=output_dir, nthreads=nthreads)

        # This assumes the BAM is coordinate sorted and indexed
        with pysam.AlignmentFile(self.vc_preprocessor.in_bam, "rb") as rs_af:
            self.total_mapped = rs_af.mapped
            self.contigs = rs_af.references

        if int(self.total_mapped) == 0:
            raise RuntimeError("No alignments to process.")

        # Divide the mapped reads by 2 to approximate pairs
        self.norm_factor = self.VARIANT_CALL_NORM_DP / (self.total_mapped / 2)

        logger.info("Loading transcript CDS annotations for AA change determination.")
        self.amino_acid_mapper = cm.AminoAcidMapper(
            gff=self.transcript_gff, ref=self.gff_reference, mut_sig=mut_sig, outdir=output_dir)

        # Get the list of unique contigs for constructing VCF headers

        if self.targets is not None:
            self.target_contigs = {str(f.chrom) for f in pybedtools.BedTool(self.targets)}
            if len(self.target_contigs) != 1:
                raise RuntimeError("Currently only one target contig/transcript is supported.")

        # Keep total R1 + R2 counts at each position for frequency calculations
        self.coordinate_counts = collections.defaultdict(int)

        # TODO: store multiple np arrays for supporting read data- one for each stat to be summarized
        # This would be p arrays of length m x n, where p is number of stats to be summarized (e.g. BQ), m is the
        # variant call key (contig:pos:ref:alt), and n is 4 for PE data (R1,+ ; R1,- ; R2,+ ; R2,-)

        # Keeps counts and stats for non-reference base supporting reads
        self.variant_counts = collections.OrderedDict()

    @staticmethod
    def _is_indel(aligned_pair):
        """Determines if a base has an InDel operation.

        :param tuple aligned_pair: query_pos, reference_pos, reference_base
        :return bool: whether or not the base is in an InDel operation
        """

        query_pos, reference_pos, reference_base = aligned_pair

        if query_pos is None or reference_pos is None:
            return True

        return False

    def _get_haplotype_ref(self, contig, first_pos, last_pos, mm_pos_set):
        r"""Determines the variant REF and local mismatch positions given the first and last MM_TUPLE positions.

        :param str contig: contig name
        :param int first_pos: 1-based position of the first mismatch
        :param int last_pos: 1-based position of the last mismatch
        :param set mm_pos_set: set of all 1-based mismatch coordinate positions for the haplotype
        :return tuple: (str, set) REF field and 0-based positions of mismatches in the REF field
        """

        ref_nts = su.extract_seq(contig, first_pos, last_pos, self.ref)
        ref_pos = list(range(first_pos, last_pos + 1))
        ref_mm_pos = [i for i, e in enumerate(ref_pos) if e in mm_pos_set]
        return ref_nts, ref_mm_pos

    def _generate_call_tuple(self, mmts):
        """Generates a CALL_TUPLE from MM_TUPLEs.

        :param iter mmts: MM_TUPLEs
        :return tuple: (collections.namedtuple, set) analysis.variant_caller.CALL_TUPLE and 1-based reference positions
        """

        first_mmt = mmts[0]
        last_mmt = mmts[-1]
        mm_pos_set = {mmt.pos for mmt in mmts}

        # Generate the REF field based on the span of the haplotype
        ref_nts, ref_mm_pos = self._get_haplotype_ref(first_mmt.contig, first_mmt.pos, last_mmt.pos, mm_pos_set)

        # Generate the ALT field based on the 0-based mismatch positions in REF
        alt_nts = list(copy.copy(ref_nts))
        for mmt, mm_pos in zip(mmts, ref_mm_pos):
            alt_nts[mm_pos] = mmt.alt

        call_tuple = CALL_TUPLE(first_mmt.contig, first_mmt.pos, ref_nts, "".join(alt_nts), None, None, None)

        return call_tuple, mm_pos_set

    def _get_haplotype_dict(self, filt_r_mms, max_mnp_window=VARIANT_CALL_MAX_MNP_WINDOW):
        r"""Finds MNPs/haplotypes within a MNP window size.

        :param list filt_r_mms: list of concordant R1 or R2 MM_TUPLEs
        :param int max_mnp_window: max number of consecutive nucleotides to search for haplotypes. Default 3.
        :return tuple: (collections.defaultdict, set) of ({collections.namedtuple: set} dict keyed by \
        CALL_TUPLE and values the set of coordinate positions included in the haplotype, and the \
        position_blacklist for mismatches involved in haplotypes)
        """

        haplotypes = collections.defaultdict(set)
        position_blacklist = set()

        filt_r_mms_len = len(filt_r_mms)
        break_index = filt_r_mms_len - 2

        # For the second to last mismatch, add a buffer MM_TUPLE with pos > the window so that i + 2 is valid
        filt_r_mms_ext = filt_r_mms + [
            MM_TUPLE(contig=None, pos=filt_r_mms[-2].pos + max_mnp_window + 1, ref=None, alt=None, bq=None, read_pos=None)]

        for i, mm_tuple in enumerate(filt_r_mms_ext):

            # If we are at the second to last mismatch, break to avoid IndexErrors
            if i > break_index and i > 0:
                break

            if mm_tuple.pos in position_blacklist:
                continue

            if filt_r_mms_ext[i + 1].pos - filt_r_mms_ext[i].pos >= max_mnp_window:
                continue

            if filt_r_mms_len == 2:
                # Make an isolated di-nt MNP call; need this first otherwise we get an IndexError on the next line
                call_tuple, mm_pos_set = self._generate_call_tuple(filt_r_mms)
                haplotypes[call_tuple] |= mm_pos_set
                position_blacklist |= mm_pos_set
                continue

            # Here we know that the next mismatch is within the window, but we don't know if its position is +1 or +2
            if filt_r_mms_ext[i + 2].pos - filt_r_mms_ext[i].pos < max_mnp_window:
                # Here we have a tri-nt MNP that spans 3 consecutive nts
                mmts = [filt_r_mms[i], filt_r_mms[i + 1], filt_r_mms[i + 2]]
                call_tuple, mm_pos_set = self._generate_call_tuple(mmts)
                haplotypes[call_tuple] |= mm_pos_set
                position_blacklist |= mm_pos_set
                continue

            # if i + 1 and i + 2 are consecutive continue as i + 1 could initiate a tri-nt MNP
            # if they are not consecutive call a di-nt MNP
            if (filt_r_mms_ext[i + 2].pos - filt_r_mms_ext[i + 1].pos) != 1:
                mmts = [filt_r_mms[i], filt_r_mms[i + 1]]
                call_tuple, mm_pos_set = self._generate_call_tuple(mmts)
                haplotypes[call_tuple] |= mm_pos_set
                position_blacklist |= mm_pos_set

        return haplotypes, position_blacklist

    def _call_haplotypes(self, filt_r_mms, max_mnp_window=VARIANT_CALL_MAX_MNP_WINDOW):
        """Calls haplotypes.

        :param list filt_r_mms: list of filtered MM_TUPLEs for either R1 or R2
        :param int max_mnp_window: max number of consecutive nucleotides to search for haplotypes. Default 3.
        :return tuple: (haplotype dict, position_blacklist), or if no haplotypes were called, None.
        """

        mm_positions = [mm.pos for mm in filt_r_mms]
        n_mm_positions = len(mm_positions)

        # We can make a simple check for reads with just two mms
        if n_mm_positions == 1:
            return None
        elif n_mm_positions == 2:
            # We can also make a simple check for reads with just two mms
            min_pos = min(mm_positions)
            max_pos = max(mm_positions)
            if max_pos - min_pos >= max_mnp_window:
                return None

        haplotypes, position_blacklist = self._get_haplotype_dict(filt_r_mms, max_mnp_window)
        return haplotypes, position_blacklist

    @staticmethod
    def _call_snps(filt_r_mms, position_blacklist, haplotypes=None):
        r"""Calls SNPs if the mismatches are not included in existing MNPs/haplotypes.

        :param list filt_r_mms: list of filtered R1 or R2 MM_TUPLEs
        :param set position_blacklist: set of mismatches involved in haplotypes, for which we don't want to \
        call SNPs
        :param dict | None haplotypes: {collections.namedtuple: set} dict keyed by CALL_TUPLE with values the set of \
        coordinate positions included in the haplotype.
        :return dict: {collections.namedtuple: set} keyed by CALL_TUPLE with values a set of coordinate positions
        for the SNP (just one, but consistent with haplotypes dict format)
        """

        snp_dict = collections.defaultdict(set)

        for mmt in filt_r_mms:
            if haplotypes is None or mmt.pos not in position_blacklist:
                snp_dict[CALL_TUPLE(mmt.contig, mmt.pos, mmt.ref, mmt.alt, None, None, None)] |= {mmt.pos}

        return snp_dict

    def _enumerate_mismatches(self, align_seg, min_bq=VARIANT_CALL_MIN_BQ):
        """Enumerates the mismatches in a read and excludes detection of InDels.

        :param pysam.AlignedSegment align_seg: read object
        :param int min_bq: min base quality for a variant to be considered as supporting a variant call
        :return list: list of MM_TUPLEs

        We require a MD tag for pysam.AlignedSegment.get_aligned_pairs(), which forms the basis of read-level calling.
        The output of get_aligned_pairs is a 3-tuple for each base.

        if query_pos is None, we are in a del operation
        if reference_pos is None, we are in an ins operation
        if reference_base is lowercase that indicates a mismatch to the reference
        """

        read_strand = su.Strand(align_seg.is_reverse)
        mms = []

        for apt in align_seg.get_aligned_pairs(with_seq=True):

            # Make sure to skip the current position if we are in an InDel operation- this method should only call
            # SNPs and MNPs
            if self._is_indel(apt):
                continue

            # Determine attrs of the current base
            query_pos, reference_pos, reference_base = apt

            # Skip mismatches with low quality
            var_qual = align_seg.query_qualities[query_pos]
            if var_qual < min_bq:
                continue

            # If we pass this, we have a mismatch to the reference
            if str(reference_base).islower():

                ref_bp = str(reference_base).upper()
                alt_bp = align_seg.query_sequence[query_pos].upper()

                if alt_bp == su.UNKNOWN_BASE:
                    continue

                # We want the variant position 1-based with respect to the 5' end
                var_pos = query_pos + 1 if read_strand == su.Strand.PLUS else align_seg.query_length - query_pos

                # Append the candidate mismatch
                mms.append(MM_TUPLE(
                    contig=align_seg.reference_name, pos=reference_pos + 1, ref=ref_bp, alt=alt_bp,
                    bq=var_qual, read_pos=var_pos))

        return mms

    @staticmethod
    def _intersect_mismatches(r1_mms, r2_mms):
        """Intersects the mismatches between read pairs.

        :param list r1_mms: list of MM_TUPLEs for R1
        :param list r2_mms: list of MM_TUPLEs for R2
        :return tuple: filtered list of R1 and R2 MM_TUPLES
        """

        r1_mm_bases = {(r1_mm.contig, r1_mm.pos, r1_mm.alt,) for r1_mm in r1_mms}
        r2_mm_bases = {(r2_mm.contig, r2_mm.pos, r2_mm.alt,) for r2_mm in r2_mms}
        r_intersect = r1_mm_bases & r2_mm_bases

        # Now we can iterate back through the list and filter it
        # First unpack the fields that will be intersected
        filtered_r1_mms = []
        filtered_r2_mms = []

        # Must use itertools.zip_longest in case differing number of mismatches are found in the pair
        for r1_mm, r2_mm in itertools.zip_longest(r1_mms, r2_mms):

            if r1_mm is not None and (r1_mm.contig, r1_mm.pos, r1_mm.alt,) in r_intersect:
                filtered_r1_mms.append(r1_mm)

            if r2_mm is not None and (r2_mm.contig, r2_mm.pos, r2_mm.alt,) in r_intersect:
                filtered_r2_mms.append(r2_mm)

        return filtered_r1_mms, filtered_r2_mms

    def _unpack_stats(self, supporting_positions, filt_r1_mms, filt_r2_mms):
        r"""Unpacks data for each position supporting the call.

        :param set supporting_positions: coordinate positions supporting the call
        :param list filt_r1_mms: list of filtered MM_TUPLEs for R1
        :param list filt_r2_mms: list of filtered MM_TUPLEs for R2
        :return tuple: (refs, alts, positions, analysis.variant_caller.PER_BP_STATS) of types \
        (str, str, str, collections.namedtuple)
        """
        
        # Extract the per-bp stats for the supporting positions
        mms_len = len(supporting_positions)
        refs = np.empty(shape=mms_len, dtype=np.str)
        alts = np.empty(shape=mms_len, dtype=np.str)
        positions = np.zeros(shape=mms_len, dtype=np.int32)
        r1_read_positions = np.zeros(shape=mms_len, dtype=np.int32)
        r2_read_positions = np.zeros(shape=mms_len, dtype=np.int32)
        r1_bqs = np.zeros(shape=mms_len, dtype=np.int32)
        r2_bqs = np.zeros(shape=mms_len, dtype=np.int32)

        # Collects stats from the lists of MM_TUPLEs for those that support the variant in question
        # (many variants could be called from the same pair)
        pos_counter = 0
        for r1_mm, r2_mm in zip(filt_r1_mms, filt_r2_mms):
            
            # We don't need to check R2 because we have already intersected; we also only need to collect the
            # REF and ALT for R1, as they are the same for R2.
            if r1_mm.pos in supporting_positions:

                refs[pos_counter] = r1_mm.ref
                alts[pos_counter] = r1_mm.alt
                positions[pos_counter] = r1_mm.pos
                r1_read_positions[pos_counter] = r1_mm.read_pos
                r1_bqs[pos_counter] = r1_mm.bq
                r2_read_positions[pos_counter] = r2_mm.read_pos
                r2_bqs[pos_counter] = r2_mm.bq
                pos_counter += 1

        # These are now static and associated with each variant call primary POS, REF, and ALT
        refs = self._STATS_DELIM.join(refs)
        alts = self._STATS_DELIM.join(alts)
        positions = self._STATS_DELIM.join(map(str, positions))

        per_bp_stats = PER_BP_STATS(r1_bqs=r1_bqs, r2_bqs=r2_bqs,
                                    r1_read_pos=r1_read_positions, r2_read_pos=r2_read_positions)

        return refs, alts, positions, per_bp_stats

    def _assign_stats(self, r_index, call_tuple, per_bp_stats):
        """Adds the BQ and read position stats to the variant counts dict.

        :param int r_index: variant counts dict read index
        :param collections.namedtuple call_tuple: CALL_TUPLE specifying the variant
        :param collections.namedtuple per_bp_stats: PER_BP_STATS of quality data for bases participating in the call
        """

        # Join the stats elements together for all positions in a call
        if r_index in {self.R1_PLUS_INDEX, self.R1_MINUS_INDEX}:
            self.variant_counts[call_tuple][r_index][self.R_BQ_INDEX].append(
                self._STATS_DELIM.join(list(map(str, per_bp_stats.r1_bqs))))

            self.variant_counts[call_tuple][r_index][self.R_RP_INDEX].append(
                self._STATS_DELIM.join(list(map(str, per_bp_stats.r1_read_pos))))

        elif r_index in {self.R2_PLUS_INDEX, self.R2_MINUS_INDEX}:
            self.variant_counts[call_tuple][r_index][self.R_BQ_INDEX].append(
                self._STATS_DELIM.join(list(map(str, per_bp_stats.r2_bqs))))

            self.variant_counts[call_tuple][r_index][self.R_RP_INDEX].append(
                self._STATS_DELIM.join(list(map(str, per_bp_stats.r2_read_pos))))

    def _add_counts_and_stats(self, call_tuple, per_bp_stats, r1_nm, r2_nm, r1_strand, r2_strand):
        """Adds counts and BQ, NM, read position stats to the call dictionary.

        :param collections.namedtuple call_tuple: CALL_TUPLE specifying the variant
        :param collections.namedtuple per_bp_stats: PER_BP_STATS of quality data for bases participating in the call
        :param int r1_nm: R1 edit distance
        :param int r2_nm: R2 edit distance
        :param analysis.seq_utils.Strand r1_strand: R1 strand
        :param analysis.seq_utils.Strand r2_strand: R2 strand
        """

        # Enumerate the supporting counts and stats
        # We will use a nested list for collection. The first list has 4 elements which correspond to each read and
        # strand pair. The inner lists will hold counts and individual base stats (e.g. BQ, NM)
        # This solution is ugly, and needs refactoring in the future. To save on memory use numpy N-dimensional arrays?
        if call_tuple not in self.variant_counts:

            # 1st position holds counts
            # Positions 2-4 contain BQ, NM, and read position data; set as lists because the number of supporting pairs
            # is unknown for each variant call candidate
            self.variant_counts[call_tuple] = [
                [0, [], [], []],  # R1, +
                [0, [], [], []],  # R1, -
                [0, [], [], []],  # R2, +
                [0, [], [], []]   # R2, -
            ]

        r1_strand_index = self.R1_MINUS_INDEX
        if r1_strand == su.Strand.PLUS:
            r1_strand_index = self.R1_PLUS_INDEX

        self._assign_stats(r1_strand_index, call_tuple, per_bp_stats)
        self.variant_counts[call_tuple][r1_strand_index][self.R_COUNTS_INDEX] += 1
        self.variant_counts[call_tuple][r1_strand_index][self.R_NM_INDEX].append(r1_nm)

        r2_strand_index = self.R2_MINUS_INDEX
        if r2_strand == su.Strand.PLUS:
            r2_strand_index = self.R2_PLUS_INDEX

        self._assign_stats(r2_strand_index, call_tuple, per_bp_stats)
        self.variant_counts[call_tuple][r2_strand_index][self.R_COUNTS_INDEX] += 1
        self.variant_counts[call_tuple][r2_strand_index][self.R_NM_INDEX].append(r2_nm)

    def _update_counts(self, collective_variants, filt_r1_mms, filt_r2_mms, r1_nm, r2_nm, r1_strand, r2_strand):
        """Updates the global dict with variant call counts and supporting read statistics.

        :param dict collective_variants: dict keyed by CALL_TUPLE and valued by set of mismatch coordinate positions
        :param list filt_r1_mms: list of filtered MM_TUPLEs for R1
        :param list filt_r2_mms: list of filtered MM_TUPLEs for R2
        :param int r1_nm: R1 edit distance
        :param int r2_nm: R2 edit distance
        :param analysis.seq_utils.Strand r1_strand: R1 strand
        :param analysis.seq_utils.Strand r2_strand: R2 strand
        """

        # Here k is a CALL_TUPLE and v is a set of coordinate positions supporting the call
        for call_tuple, supporting_positions in collective_variants.items():

            # First unpack useful statistics for each base participating in the variant call
            refs, alts, positions, per_bp_stats = self._unpack_stats(supporting_positions, filt_r1_mms, filt_r2_mms)

            # Update the call tuple key to put the component REFs, ALTs, and positions
            new_call_tuple = CALL_TUPLE(
                contig=call_tuple.contig, pos=call_tuple.pos, ref=call_tuple.ref, alt=call_tuple.alt,
                refs=refs, alts=alts, positions=positions)

            # Then store the counts and stats into a dict to facilitate summation across reads
            self._add_counts_and_stats(new_call_tuple, per_bp_stats, r1_nm, r2_nm, r1_strand, r2_strand)

    @staticmethod
    def _get_unmasked_positions(align_seg):
        """Gets the aligned reference positions of the read that are not BQ-masked.

        :param pysam.AlignedSegment align_seg: read object
        :return set: set of reference positions
        """

        unmasked_positions = set()

        for ref_pos, bq in zip(align_seg.get_reference_positions(), align_seg.query_alignment_qualities):

            if bq == su.MASKED_BQ:
                continue

            unmasked_positions.add(ref_pos)

        return unmasked_positions

    @staticmethod
    def _reads_overlap(r1, r2):
        """Determines if paired reads overlap.

        :param pysam.AlignedSegment r1: R1 read object
        :param pysam.AlignedSegment r2: R2 read object
        :return bool: whether or not the reads overlap (require >= 1 nt overlap)
        """

        r1_ref_pos = set(r1.get_reference_positions())
        r2_ref_pos = set(r2.get_reference_positions())

        if len(r1_ref_pos & r2_ref_pos) > 0:
            return True

        return False

    def _update_pos_dp(self, r1, r2):
        """Updates the reference position dict for fragment coverage.

        :param pysam.AlignedSegment r1: R1 read object
        :param pysam.AlignedSegment r2: R2 read object
        :return None: if no unmasked positions exist
        """

        # Only update the DP for positions with nonzero BQ
        ref_positions = self._get_unmasked_positions(r1)
        ref_positions |= self._get_unmasked_positions(r2)

        if len(ref_positions) == 0:
            return None

        # Reference positions are 0-based
        # pysam says reference_end (last element in positions) points to one past the last
        # aligned base; but this is not true, so add 1 to the range end
        min_ref_pos = min(ref_positions)
        max_ref_pos = max(ref_positions)
        fragment_pos = range(min_ref_pos, max_ref_pos + 1)

        for pos in fragment_pos:
            # Enumerate the read pair for the DP denominator to frequency
            self.coordinate_counts[COORDINATE_KEY(r1.reference_name, pos + 1)] += 1

    def _iterate_over_reads(self, af1, af2, min_bq=VARIANT_CALL_MIN_BQ, max_nm=VARIANT_CALL_MAX_NM,
                            max_mnp_window=VARIANT_CALL_MAX_MNP_WINDOW):
        """Iterates over read pairs to enumerate variants.

        :param pysam.AlignmentFile af1: object corresponding to the R1 BAM
        :param pysam.AlignmentFile af2: object corresponding to the R2 BAM
        :param int min_bq: min base quality
        :param int max_nm: max edit distance (NM tag) to consider a read for variant calls
        :param int max_mnp_window: max number of consecutive nucleotides to search for haplotypes
        """

        # Because we must start with qname-sorted BAMs, we can't extract reads particular to a genomic region without
        # first intersecting
        for r1, r2 in zip(af1.fetch(until_eof=True), af2.fetch(until_eof=True)):

            # Sanity check to make sure we are always paired
            # Split to handle consensus deduplicated input which include the mate ID along with the group ID
            if r1.query_name != r2.query_name:
                raise RuntimeError(
                    "Improper pairing of reads. R1 was %s and R2 was %s." % (r1.query_name, r2.query_name))

            # Only call variants for pairs that pass filters
            r1_nm = su.get_edit_distance(r1)
            r2_nm = su.get_edit_distance(r2)
            if r1_nm > max_nm or r2_nm > max_nm:
                continue

            if not self._reads_overlap(r1, r2):
                continue

            # Compute fragment coverage/depth; note only filtered read pairs contribute to depth
            self._update_pos_dp(r1, r2)

            # Enumerate mismatch positions for each read in the pair
            r1_mms = self._enumerate_mismatches(r1, min_bq)
            r2_mms = self._enumerate_mismatches(r2, min_bq)

            # Now find intersections between the mismatches
            filt_r1_mms, filt_r2_mms = self._intersect_mismatches(r1_mms, r2_mms)

            # If we found no intersected mismatches, continue to the next read pair
            if len(filt_r1_mms) == 0:
                continue

            # Call MNPs and haplotypes
            haplotype_res = self._call_haplotypes(filt_r1_mms, max_mnp_window)

            if isinstance(haplotype_res, tuple):
                haplotypes, position_blacklist = haplotype_res
                if len(position_blacklist) == 0:
                    # In this case we had >= 3 mismatches but they were not within the window and position_blacklist
                    # is an empty set
                    haplotypes = None
            else:
                # In this case we had either 1 mismatch or 2 mismatches that were not within the window
                haplotypes = None
                position_blacklist = set()

            # Call SNPs not captured in a MNP or haplotype
            snps = self._call_snps(filt_r1_mms, position_blacklist, haplotypes)

            collective_variants = snps
            if haplotypes is not None:
                collective_variants = {**haplotypes, **snps}

            r1_strand = su.Strand(r1.is_reverse)
            r2_strand = su.Strand(r2.is_reverse)

            # Now that we have called haplotypes and SNPs for the pair, update the dict of counts and stats
            self._update_counts(collective_variants, filt_r1_mms, filt_r2_mms, r1_nm, r2_nm, r1_strand, r2_strand)

            # InDel enumeration can exist as a separate call here in the future

    def _summarize_stats(self, var_list, read_index, stat_index, pos_index):
        """Summarizes per-bp stats across reads for a particular contributing base to a variant call.

        :param list var_list: nested list containing data for each read/strand pair
        :param int read_index: which read/strand pair should be summarized?
        :param int stat_index: which stat should be summarized? Should be either R_BQ_INDEX or R_RP_INDEX
        :param int pos_index: index of the contributing mismatch to the variant call
        :return list: summarized stat for each position contributing to the variant
        """

        read_list = var_list[read_index]
        stat_list = read_list[stat_index]  # this is a list of comma-delimited stats for each read supporting the call
        position_stats = []  # this will collect a position's stat across multiple reads

        # Now iterate over the stats for each read and recollect stats across reads
        for read_stats in stat_list:

            # Cast the joined values back to int
            stat_split = list(map(int, read_stats.split(self._STATS_DELIM)))
            position_stats.append(stat_split[pos_index])

        # Once we have collected the stats for the specified position across reads, summarize them
        position_stat_summary = str(statistics.median(position_stats)) if len(position_stats) > 0 else su.R_COMPAT_NA

        return position_stat_summary

    def _get_read_pos_stats(self, var_list, pos_index):
        """Gets the median read positions for each read/strand pair.

        :param list var_list: nested list containing data for each read/strand pair
        :param int pos_index: index of the contributing mismatch to the variant call
        :return tuple: tuple of str (r1_plus_rp, r1_minus_rp, r2_plus_rp, r2_minus_rp)
        """

        read_indices = (self.R1_PLUS_INDEX, self.R1_MINUS_INDEX, self.R2_PLUS_INDEX, self.R2_MINUS_INDEX)
        rp_stats = [self._summarize_stats(var_list, ri, self.R_RP_INDEX, pos_index) for ri in read_indices]
        return tuple(rp_stats)

    def _get_read_bq_stats(self, var_list, pos_index):
        """Gets the median BQs for each read/strand pair.

        :param list var_list: nested list containing data for each read/strand pair
        :param int pos_index: index of the contributing mismatch to the variant call
        :return tuple: tuple of str (r1_plus_bq, r1_minus_bq, r2_plus_bq, r2_minus_bq)
        """

        read_indices = (self.R1_PLUS_INDEX, self.R1_MINUS_INDEX, self.R2_PLUS_INDEX, self.R2_MINUS_INDEX)
        bq_stats = [self._summarize_stats(var_list, ri, self.R_BQ_INDEX, pos_index) for ri in read_indices]
        return tuple(bq_stats)

    def _get_read_nm_stats(self, var_list):
        """Gets the median NM for each read/strand pair.

        :param list var_list: nested list containing data for each read/strand pair
        :return tuple: tuple of str (r1_plus_nm, r1_minus_nm, r2_plus_nm, r2_minus_nm)
        """

        read_indices = (self.R1_PLUS_INDEX, self.R1_MINUS_INDEX, self.R2_PLUS_INDEX, self.R2_MINUS_INDEX)

        nm_stats = [str(statistics.median(var_list[ri][self.R_NM_INDEX])) if len(var_list[ri][self.R_NM_INDEX]) > 0
                    else su.R_COMPAT_NA for ri in read_indices]

        return tuple(nm_stats)

    def _call_variants(self, min_supporting_qnames=VARIANT_CALL_MIN_DP):
        r"""Calls variants and summarizes stats for bases contributing to the variant calls.

        :param int min_supporting_qnames: min fragments with R1 and R2 supporting coverage for which to keep a variant
        :return collections.OrderedDict: dict keyed same as self.variant_calls, but containing only passing variants
        """

        # Keep track of concordant counts
        concordant_counts_dict = collections.OrderedDict()

        # TODO: Consider vectorizing or parallelizing per-variant operations
        # k = collections.namedtuple("CALL_TUPLE", "contig, pos, ref, alt, refs, alts, positions")
        # v is nested list: 4 elements for each read/strand pair, with internal 4-element lists containing
        # counts and lists of BQs, RPs, NMs with each element specifying a supporting read
        for k, v in self.variant_counts.items():

            # We just count one of the mates as we are getting fragment-based counts and frequencies
            s1_cao = v[self.R1_PLUS_INDEX][self.R_COUNTS_INDEX]
            s2_cao = v[self.R1_MINUS_INDEX][self.R_COUNTS_INDEX]
            cao = s1_cao + s2_cao

            # Only call variants if they exceed the threshold
            if cao < min_supporting_qnames:
                continue

            # Unpack the component positions, REF bases, and ALT bases in each primary POS, REF, ALT call
            var_refs = k.refs.split(self._STATS_DELIM)
            var_alts = k.alts.split(self._STATS_DELIM)
            var_pos = list(map(int, k.positions.split(self._STATS_DELIM)))
            var_nm_stats = self._get_read_nm_stats(v)

            # For MNPs use the floor of the contributed DP to fix the denominator
            ref_len = len(k.ref)
            if ref_len >= 2:

                dp_list = [self.coordinate_counts[COORDINATE_KEY(k.contig, pos)]
                           for pos in range(k.pos, k.pos + ref_len) if pos in set(var_pos)]

                dp = min(dp_list)
            else:
                dp = self.coordinate_counts[COORDINATE_KEY(k.contig, k.pos)]

            norm_cao = cao * self.norm_factor
            caf = cao / dp

            # Report stats for each component base in variant (for delineating contributing bases to MNPs/haplotypes)
            for i in range(len(var_pos)):

                # We will split out MNPs into multiple records, for each component base
                new_key = VARIANT_CALL_KEY_TUPLE(contig=k.contig, pos=k.pos, ref=k.ref, alt=k.alt, index=i)

                # Extract and summarize the stats for each contributing base to the call
                # For read/strand pairs with no support of the call, the values will be NA
                var_rp_stats = self._get_read_pos_stats(v, i)
                var_bq_stats = self._get_read_bq_stats(v, i)

                upstream_base = su.extract_seq(
                    contig=k.contig, start=var_pos[i] - 1, stop=var_pos[i] - 1, ref=self.ref)

                downstream_base = su.extract_seq(
                    contig=k.contig, start=var_pos[i] + 1, stop=var_pos[i] + 1, ref=self.ref)

                # Finally store the summarized data with the original coordinate and REF-ALT specific key
                vcst = VARIANT_CALL_SUMMARY_TUPLE(
                    POS_NT=var_pos[i], REF_NT=var_refs[i], ALT_NT=var_alts[i],
                    UP_REF_NT=upstream_base, DOWN_REF_NT=downstream_base,
                    DP=dp, CAO=cao, CAF=caf, NORM_CAO=norm_cao,
                    R1_PLUS_AO=v[self.R1_PLUS_INDEX][self.R_COUNTS_INDEX],
                    R1_MINUS_AO=v[self.R1_MINUS_INDEX][self.R_COUNTS_INDEX],
                    R2_PLUS_AO=v[self.R2_PLUS_INDEX][self.R_COUNTS_INDEX],
                    R2_MINUS_AO=v[self.R2_MINUS_INDEX][self.R_COUNTS_INDEX],
                    R1_PLUS_MED_RP=var_rp_stats[self.R1_PLUS_INDEX],
                    R1_MINUS_MED_RP=var_rp_stats[self.R1_MINUS_INDEX],
                    R2_PLUS_MED_RP=var_rp_stats[self.R2_PLUS_INDEX],
                    R2_MINUS_MED_RP=var_rp_stats[self.R2_MINUS_INDEX],
                    R1_PLUS_MED_BQ=var_bq_stats[self.R1_PLUS_INDEX],
                    R1_MINUS_MED_BQ=var_bq_stats[self.R1_MINUS_INDEX],
                    R2_PLUS_MED_BQ=var_bq_stats[self.R2_PLUS_INDEX],
                    R2_MINUS_MED_BQ=var_bq_stats[self.R2_MINUS_INDEX],
                    R1_PLUS_MED_NM=var_nm_stats[self.R1_PLUS_INDEX],
                    R1_MINUS_MED_NM=var_nm_stats[self.R1_MINUS_INDEX],
                    R2_PLUS_MED_NM=var_nm_stats[self.R2_PLUS_INDEX],
                    R2_MINUS_MED_NM=var_nm_stats[self.R2_MINUS_INDEX]
                    )

                concordant_counts_dict[new_key] = vcst

        return concordant_counts_dict

    def _write_results(self, concordant_counts, cov_fh, reference_candidates_fh):
        r"""Writes output VCF and coverage BED.

        :param collections.OrderedDict concordant_counts: dict containing only passing variants
        :param file cov_fh: output BED file handle to write fragment coverage
        :param pysam.VariantFile reference_candidates_fh: output VCF file handle to write candidate variants
        """

        # Dump the concordant counts
        for k, v in concordant_counts.items():
            # k is VARIANT_CALL_KEY_TUPLE, v is VARIANT_CALL_SUMMARY_TUPLE
            self._write_concordant_variants(k, v, reference_candidates_fh)

        # And the coverage file in BED6 format
        for k, v in self.coordinate_counts.items():
            # total_counts/BED score: bases mapped for R1 and R2 (no InDels)
            total_counts = self.coordinate_counts[COORDINATE_KEY(k.contig, k.pos)]
            cov_line = [k.contig, k.pos - 1, k.pos, total_counts]
            cov_fh.write(fu.FILE_DELIM.join(list(map(str, cov_line))) + fu.FILE_NEWLINE)

    def _write_concordant_variants(self, vckt, vcst, reference_candidates_fh):
        """Writes a concordant variant to VCF.

        :param collections.namedtuple vckt: VARIANT_CALL_KEY_TUPLE
        :param collections.namedtuple vcst: VARIANT_CALL_SUMMARY_TUPLE
        :param pysam.VariantFile reference_candidates_fh: output VCF file handle to write candidate variants
        """

        trx_id = vckt.contig.split(APPRIS_CONTIG_DELIM)[APPRIS_TRX_INDEX]

        # Determine the codon, AA change(s) and whether or not the variant matches the mutagenesis signature
        var_location, var_wt_codons, var_mut_codons, var_wt_aas, var_mut_aas, var_aa_changes, var_aa_positions, \
        var_matches_mut_sig = self.amino_acid_mapper.get_codon_and_aa_changes(
            trx_id=trx_id, pos=vckt.pos, ref=vckt.ref, alt=vckt.alt)

        # Create an INFO field for the VCF
        ri = {
            vu.VCF_POS_NT_ID: vcst.POS_NT,
            vu.VCF_REF_NT_ID: vcst.REF_NT,
            vu.VCF_ALT_NT_ID: vcst.ALT_NT,
            vu.VCF_UP_REF_NT_ID: vcst.UP_REF_NT,
            vu.VCF_DOWN_REF_NT_ID: vcst.DOWN_REF_NT,
            vu.VCF_DP_ID: vcst.DP,
            vu.VCF_CAO_ID: vcst.CAO,
            vu.VCF_NORM_CAO_ID: vcst.NORM_CAO,
            vu.VCF_CAF_ID: vcst.CAF,
            vu.VCF_R1_PLUS_AO_ID: vcst.R1_PLUS_AO,
            vu.VCF_R1_MINUS_AO_ID: vcst.R1_MINUS_AO,
            vu.VCF_R2_PLUS_AO_ID: vcst.R2_PLUS_AO,
            vu.VCF_R2_MINUS_AO_ID: vcst.R2_MINUS_AO,
            vu.VCF_R1_PLUS_MED_POS_ID: vcst.R1_PLUS_MED_RP,
            vu.VCF_R1_MINUS_MED_POS_ID: vcst.R1_MINUS_MED_RP,
            vu.VCF_R2_PLUS_MED_POS_ID: vcst.R2_PLUS_MED_RP,
            vu.VCF_R2_MINUS_MED_POS_ID: vcst.R2_MINUS_MED_RP,
            vu.VCF_R1_PLUS_MED_BQ_ID: vcst.R1_PLUS_MED_BQ,
            vu.VCF_R1_MINUS_MED_BQ_ID: vcst.R1_MINUS_MED_BQ,
            vu.VCF_R2_PLUS_MED_BQ_ID: vcst.R2_PLUS_MED_BQ,
            vu.VCF_R2_MINUS_MED_BQ_ID: vcst.R2_MINUS_MED_BQ,
            vu.VCF_R1_PLUS_MED_NM_ID: vcst.R1_PLUS_MED_NM,
            vu.VCF_R1_MINUS_MED_NM_ID: vcst.R1_MINUS_MED_NM,
            vu.VCF_R2_PLUS_MED_NM_ID: vcst.R2_PLUS_MED_NM,
            vu.VCF_R2_MINUS_MED_NM_ID: vcst.R2_MINUS_MED_NM,
            vu.VCF_AAM_LOCATION_ID: var_location,
            vu.VCF_AAM_CODON_REF_ID: var_wt_codons,
            vu.VCF_AAM_CODON_ALT_ID: var_mut_codons,
            vu.VCF_AAM_AA_REF_ID: var_wt_aas,
            vu.VCF_AAM_AA_ALT_ID: var_mut_aas,
            vu.VCF_AAM_AA_CHANGE_ID: var_aa_changes,
            vu.VCF_AAM_AA_POS_ID: var_aa_positions,
            vu.VCF_MUT_SIG_MATCH: str(var_matches_mut_sig)
        }

        # Generate a VCF record from scracth
        new_variant_record = reference_candidates_fh.new_record(
            contig=vckt.contig, start=vckt.pos - 1, alleles=(vckt.ref, vckt.alt), info=ri)

        reference_candidates_fh.write(new_variant_record)

    def _create_vcf_header(self):
        """Creates a VCF header for the candidate variants.

        :return pysam.VariantHeader: VCF header object
        """

        vcf_header = pysam.VariantHeader()

        class_module_path = ".".join(
            [PROJECT_LAB, PROJECT_AUTHOR, os.path.basename(PROJECT_ROOT), __name__, self.__class__.__name__])

        vcf_header.add_line(vu.VCF_METADATA_CHAR + "source=" + class_module_path)

        # Add the contig ID(s)
        for contig in self.contigs:
            vcf_header.add_line(vu.VCF_CONTIG_HEADER_FORMAT.format(contig))

        aa_mapper_module_path = ".".join(
            [PROJECT_LAB, PROJECT_AUTHOR, os.path.basename(PROJECT_ROOT), cm.__name__, cm.AminoAcidMapper.__name__])

        # Tuples of (ID, Number, Type, Description)
        # Consider putting these in a config file
        info_metadata = [
            (vu.VCF_POS_NT_ID, 1, "Integer", "Coordinate position of component nucleotide."),
            (vu.VCF_REF_NT_ID, ".", "String", "Component reference nucleotide."),
            (vu.VCF_ALT_NT_ID, ".", "String", "Component alternate nucleotide."),
            (vu.VCF_UP_REF_NT_ID, 1, "String", "-1 upstream reference nucleotide."),
            (vu.VCF_DOWN_REF_NT_ID, 1, "String", "+1 downstream reference nucleotide."),
            (vu.VCF_DP_ID, 1, "Integer", "Fragment-based depth of coverage after quality filters (pair overlap, edit distance)."),
            (vu.VCF_CAO_ID, 1, "Integer", "Concordant alternate observations- alternate found in both mates."),
            (vu.VCF_NORM_CAO_ID, 1, "Float", "Mate-concordant observations per %i pairs." % self.VARIANT_CALL_NORM_DP),
            (vu.VCF_CAF_ID, 1, "Float", "Concordant allele frequency in range (0,1). Calculated as CAO/DP."),
            (vu.VCF_R1_PLUS_AO_ID, 1, "Integer", "Read 1 alternate observations on (+) strand."),
            (vu.VCF_R1_MINUS_AO_ID, 1, "Integer", "Read 1 alternate observations on (-) strand."),
            (vu.VCF_R2_PLUS_AO_ID, 1, "Integer", "Read 2 alternate observations on (+) strand."),
            (vu.VCF_R2_MINUS_AO_ID, 1, "Integer", "Read 2 alternate observations on (-) strand."),
            (vu.VCF_R1_PLUS_MED_POS_ID, ".", "String", "Read 1 (+) strand median read position supporting call."),
            (vu.VCF_R1_MINUS_MED_POS_ID, ".", "String", "Read 1 (-) strand median read position supporting call."),
            (vu.VCF_R2_PLUS_MED_POS_ID, ".", "String", "Read 2 (+) strand median read position supporting call."),
            (vu.VCF_R2_MINUS_MED_POS_ID, ".", "String", "Read 2 (-) strand median read position supporting call."),
            (vu.VCF_R1_PLUS_MED_BQ_ID, ".", "String", "Read 1 (+) strand median Phred base quality supporting call."),
            (vu.VCF_R1_MINUS_MED_BQ_ID, ".", "String", "Read 1 (-) strand median Phred base quality supporting call."),
            (vu.VCF_R2_PLUS_MED_BQ_ID, ".", "String", "Read 2 (+) strand median Phred base quality supporting call."),
            (vu.VCF_R2_MINUS_MED_BQ_ID, ".", "String", "Read 2 (-) strand median Phred base quality supporting call."),
            (vu.VCF_R1_PLUS_MED_NM_ID, ".", "String", "Read 1 (+) strand median edit distance supporting call."),
            (vu.VCF_R1_MINUS_MED_NM_ID, ".", "String", "Read 1 (-) strand median edit distance supporting call."),
            (vu.VCF_R2_PLUS_MED_NM_ID, ".", "String", "Read 2 (+) strand median edit distance supporting call."),
            (vu.VCF_R2_MINUS_MED_NM_ID, ".", "String", "Read 2 (-) strand median edit distance supporting call."),
            (vu.VCF_AAM_LOCATION_ID, ".", "String",
             "%s location of the variant in the transcript. One of {CDS, 5_UTR, 3_UTR, intergenic, untranslated}."
             % aa_mapper_module_path),
            (vu.VCF_AAM_CODON_REF_ID, ".", "String",
             "%s comma-delimited reference codon(s). NA if the variant is out of CDS bounds." % aa_mapper_module_path),
            (vu.VCF_AAM_CODON_ALT_ID, ".", "String",
             "%s comma-delimited alternate codon(s). NA if the variant is out of CDS bounds." % aa_mapper_module_path),
            (vu.VCF_AAM_AA_REF_ID, ".", "String",
             "%s comma-delimited reference amino acid(s). NA if the variant is out of CDS bounds." % aa_mapper_module_path),
            (vu.VCF_AAM_AA_ALT_ID, ".", "String",
             "%s comma-delimited alternate amino acid(s). NA if the variant is out of CDS bounds." % aa_mapper_module_path),
            (vu.VCF_AAM_AA_CHANGE_ID, ".", "String",
             "%s comma-delimited amino acid change(s). NA if the variant is out of CDS bounds." % aa_mapper_module_path),
            (vu.VCF_AAM_AA_POS_ID, ".", "String",
             "%s comma-delimited amino acid position(s). NA if the variant is out of CDS bounds." % aa_mapper_module_path),
            (vu.VCF_MUT_SIG_MATCH, ".", "String", "Whether or not the variant matches the mutagenesis signature.")
        ]

        # Add the INFO fields to be populated in the variant records
        for info_id in info_metadata:
            vcf_header.add_line(vu.VCF_INFO_HEADER_FORMAT.format(*info_id))

        return vcf_header

    def workflow(self, min_bq=VARIANT_CALL_MIN_BQ, max_nm=VARIANT_CALL_MAX_NM,
                 min_supporting_qnames=VARIANT_CALL_MIN_DP, max_mnp_window=VARIANT_CALL_MAX_MNP_WINDOW,
                 out_prefix=VARIANT_CALL_PREFIX):
        """Executes the variant calling workflow with specified quality parameters and count thresholds.

        :param int min_bq: min base quality
        :param int max_nm: max edit distance (NM tag) to consider a read for variant calls
        :param int min_supporting_qnames: min number of fragments with R1-R2 concordant coverage to keep a variant
        :param int max_mnp_window: max number of consecutive nucleotides to search for MNPs; must be between 1 and 3.
        :param str out_prefix: output directory and filename prefix to write results to.
        :return tuple: (VCF, BED) filepaths
        :raises NotImplementedError: if min_bq is 0 while primers are provided, or the max_mnp_window is < 3
        """

        logger.info("Starting variant calling workflow.")

        if self.primers is not None and min_bq == 0:
            raise NotImplementedError("If primers are provided, min_bq must be >= 1 so that synthetic sequences "
                                      "can be detected.")

        if max_mnp_window not in {1, 2, 3}:
            raise NotImplementedError("--max_mnp_window must be one of {1,2,3}.")

        # Create some temp VCFs to use in patch for removing pysam's obligatory END INFO tag addition,
        # which interferes with IGV visualization and is only applicable for structural variants in VCF.
        patch_reference = tempfile.NamedTemporaryFile(suffix=".patch.ref.vcf", delete=False).name

        reference_vcf = fu.add_extension(out_prefix, self.VARIANT_CALL_REF_CANDIDATE_EXT)
        reference_bed = fu.add_extension(out_prefix, self.VARIANT_CALL_COV_EXT)

        # The headers are specific to the type because they include the contig names.
        reference_vcf_header = self._create_vcf_header()

        # Quick fix to ignore internal htslib errors
        # see https://github.com/pysam-developers/pysam/issues/939
        # verbosity_save = pysam.set_verbosity(0)

        with pysam.AlignmentFile(self.vc_preprocessor.r1_calling_bam, "rb", check_sq=False) as af1, \
                pysam.AlignmentFile(self.vc_preprocessor.r2_calling_bam, "rb", check_sq=False) as af2, \
                open(reference_bed, "w") as cov_fh, \
                pysam.VariantFile(patch_reference, "w", header=reference_vcf_header) as reference_candidates_fh:

            logger.info("Collecting read mismatch data. This may take some time...")
            self._iterate_over_reads(af1, af2, min_bq, max_nm, max_mnp_window)

            logger.info("Calling variants.")
            concordant_counts = self._call_variants(min_supporting_qnames)

            logger.info("Writing results.")
            self._write_results(concordant_counts, cov_fh, reference_candidates_fh)

        # Run the patch to remove the INFO END tag, which interferes with visualization of VCFs in IGV
        temp_files = [patch_reference]
        if self.targets is None:
            vu.remove_end_info_tag(in_vcf=patch_reference, out_vcf=reference_vcf)
        else:
            patch_temp = vu.remove_end_info_tag(in_vcf=patch_reference)
            # Need to write the header as this can cause issues with opening the file in table_from_vcf()
            intersect_features(ff1=patch_temp, ff2=self.targets, outfile=reference_vcf, as_bedtool=False, header=True)
            temp_files.append(patch_temp)

        # Create a summary table for the VCF
        vu.table_from_vcf(vcf=reference_vcf)

        fu.safe_remove(tuple(temp_files))
        # pysam.set_verbosity(verbosity_save)

        logger.info("Completed variant calling workflow.")

        return reference_vcf, reference_bed
