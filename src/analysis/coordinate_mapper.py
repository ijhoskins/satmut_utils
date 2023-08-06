#!/usr/bin/env python3
"""Objects for mapping between transcriptomic and protein coordinates."""

import collections
import gzip
import logging
import pickle
import pybedtools
import re
import tempfile
import warnings

import analysis.seq_utils as su
import core_utils.feature_file_utils as ffu
import core_utils.file_utils as fu
import core_utils.vcf_utils as vu
from satmut_utils.definitions import *

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

CODONS = ("GCC", "GCT", "GCA", "GCG",
          "TGC", "TGT",
          "GAC", "GAT",
          "GAG", "GAA",
          "TTC", "TTT",
          "GGC", "GGG", "GGA", "GGT",
          "CAC", "CAT",
          "ATC", "ATT", "ATA",
          "AAG", "AAA",
          "CTG", "CTC", "TTG", "CTT", "CTA", "TTA",
          "ATG",
          "AAC", "AAT",
          "CCC", "CCT", "CCA", "CCG",
          "CAG", "CAA",
          "CGC", "AGG", "CGG", "AGA", "CGA", "CGT",
          "AGC", "TCC", "TCT", "AGT", "TCA", "TCG",
          "ACC", "ACA", "ACT", "ACG",
          "GTG", "GTC", "GTT", "GTA",
          "TGG",
          "TAC", "TAT",
          "TGA", "TAA", "TAG")

CODON_AAS = ["A"] * 4 + ["C"] * 2 + ["D"] * 2 + ["E"] * 2 + ["F"] * 2 + ["G"] * 4 + ["H"] * 2 + ["I"] * 3 + ["K"] * 2 + \
["L"] * 6 + ["M"] + ["N"] * 2 + ["P"] * 4 + ["Q"] * 2 + ["R"] * 6 + ["S"] * 6 + ["T"] * 4 + ["V"] * 4 + ["W"] + \
            ["Y"] * 2 + ["*"] * 3

CODON_AA_DICT = dict(zip(CODONS, CODON_AAS))
STOP_CODONS = {"TGA", "TAA", "TAG"}

HGVS_AA_FORMAT = "p.{}{}{}"
MUT_SIG_UNEXPECTED_WOBBLE_BPS = {"NNN": set(), "NNK": {"A", "C"}, "NNS": {"A", "T"}}

EXON_COORDS_TUPLE = collections.namedtuple("EXON_COORDS_TUPLE", "exon_id, contig, start, stop, exon_len, strand")

MUT_INFO_TUPLE = collections.namedtuple(
    "MUT_INFO_TUPLE", "location, wt_codons, mut_codons, wt_aas, mut_aas, aa_changes, aa_positions, matches_mut_sig")

CODON_TUPLE = collections.namedtuple("CODON_TUPLE", "codon_pos, codon, base_index")

tempfile.tempdir = DEFAULT_TEMPDIR
logger = logging.getLogger(__name__)


class TranscriptNotFound(Exception):
    """Exception for when a transcript ID was not found in the CoordinateMapper."""
    pass


class TranscriptomicCoordNotFound(Warning):
    """Warning for when a coordinate is outside of transcript bounds."""
    pass


class MapperBase(object):
    """Base class for GTF/GFF related data parsing."""

    # Some GFF/GTF type fields
    # Common feature types
    GENE_ID = "gene"
    TRX_ID = "transcript"
    PRIM_TRX_ID = "primary_transcript"
    MRNA_ID = "mRNA"
    CDS_ID = "CDS"
    EXON_ID = "exon"
    START_CODON_ID = "start_codon"
    STOP_CODON_ID = "stop_codon"
    UTR_ID = "UTR"

    ENSEMBL_GENE_PREFIX = "ENSG"
    ENSEMBL_TRX_PREFIX = "ENST"
    GFF_FEATURE_FIELD = 2

    # RNA feature types
    R_RNA_ID = "rRNA"
    T_RNA_ID = "tRNA"
    MI_RNA_ID = "miRNA"
    SN_RNA_ID = "snRNA"
    SNO_RNA_ID = "snoRNA"
    NC_RNA_ID = "ncRNA"
    LNC_RNA_ID = "lnc_RNA"
    TELO_RNA_ID = "telomerase_RNA"
    DLOOP_ID = "D_loop"
    VAULT_RNA_ID = "vault_RNA"
    ANTIS_RNA_ID = "antisense_RNA"
    Y_RNA = "Y_RNA"
    RNASE_MRP_RNA_ID = "RNase_MRP_RNA"
    RNASE_P_RNA_ID = "RNase_P_RNA"
    SRP_RNA_ID = "SRP_RNA"

    # DNA feature types
    PROMOTER_ID = "promoter"
    ENHANCER_ID = "enhancer"

    # Misc. feature types
    REPEAT_ID = "repeat_region"
    REGION_ID = "region"
    SEQ_FEAT_ID = "sequence_feature"
    MATCH_ID = "match"
    CDNA_MATCH_ID = "cDNA_match"
    SELENOCYSTEINE_ID = "Selenocysteine"

    # Immunoglobulin locus feature types
    V_GENE_SEG_ID = "V_gene_segment"
    D_GENE_SEG_ID = "D_gene_segment"
    J_GENE_SEG_ID = "J_gene_segment"
    C_GENE_SEG_ID = "C_gene_segment"

    # Omit transcripts on non-NC contigs
    CONTIGS_TO_OMIT = re.compile("|".join(["NT", "NW"]))

    FEATURES_TO_OMIT = [R_RNA_ID, T_RNA_ID, MI_RNA_ID, SN_RNA_ID, SNO_RNA_ID, NC_RNA_ID, LNC_RNA_ID, TELO_RNA_ID,
                        DLOOP_ID, VAULT_RNA_ID, ANTIS_RNA_ID, Y_RNA, RNASE_MRP_RNA_ID, RNASE_P_RNA_ID, SRP_RNA_ID,
                        PROMOTER_ID, ENHANCER_ID, REPEAT_ID, REGION_ID, SEQ_FEAT_ID, MATCH_ID, CDNA_MATCH_ID,
                        V_GENE_SEG_ID, D_GENE_SEG_ID, J_GENE_SEG_ID, C_GENE_SEG_ID, SELENOCYSTEINE_ID]


class AminoAcidMapper(MapperBase):
    """Class for determining the codon and AA change(s) given a transcriptomic variant."""

    FEATURES_TO_OMIT = re.compile("|".join(MapperBase.FEATURES_TO_OMIT +
                                           [MapperBase.GENE_ID, MapperBase.TRX_ID, MapperBase.MRNA_ID,
                                            MapperBase.UTR_ID, MapperBase.START_CODON_ID]))

    VALID_CDS_FEATURES = {MapperBase.CDS_ID, MapperBase.STOP_CODON_ID}

    TRX_LEN = "trx_len"
    CDS_START_OFFSET = "cds_start_offset"  # local transcript offset from start of transcript to start of start codon
    CDS_STOP_OFFSET = "cds_stop_offset"  # local transcript offset from start of transcript to start of stop codon
    CDS_NT_SEQ = "cds_nt_seq"
    FIVEPRIME_UTR = "5_UTR"
    THREEPRIME_UTR = "3_UTR"
    INTERGENIC = "intergenic"
    UNTRANSLATED = "untranslated"
    MUT_INFO_DELIM = ","
    DEFAULT_OUTDIR = "."
    MUT_SIG_ANY = "NNN"
    CDS_INFO_TRX_LEN_INDEX = 0
    CDS_INFO_CDS_START_INDEX = 1
    CDS_INFO_CDS_STOP_INDEX = 2
    CDS_INFO_TRX_SEQ_INDEX = 3
    CDS_INFO_CDS_SEQ_INDEX = 4
    DEFAULT_ARGS = MUT_INFO_TUPLE._fields[1:]
    DEFAULT_KWARGS = dict(zip(DEFAULT_ARGS, [su.R_COMPAT_NA] * len(DEFAULT_ARGS)))
    NONSTOP_CHAR = "X"

    def __init__(self, gff, ref, outdir=DEFAULT_OUTDIR, use_pickle=True, make_pickle=True, overwrite_pickle=False,
                 mut_sig=MUT_SIG_ANY, filter_unexpected=False):
        r"""Constructor for AminoAcidMapper.

        :param str gff: GFF/GTF to create a mapper for; must have "transcript_id" and "CDS", "stop_codon" features.
        :param str ref: reference FASTA corresponding to GFF features
        :param str outdir: output directory to write pickles to, if make_pickle=True.
        :param bool use_pickle: Use a previously generated pickle for annotations if one exists? Default True. Otherwise \
        re-create the data structure.
        :param bool make_pickle: Should a new pickle be made? Default True.
        :param bool overwrite_pickle: if make_pickle=True and use_pickle=True, should the existing pickle be \
        overwritten? Default False. This is useful if one makes dynamic edits to an existing GFF. Otherwise, the older \
        pickle will be used.
        :param str mut_sig: mutagenesis signature- one of {NNN, NNK, NNS}. Default NNN.
        :param bool filter_unexpected: Filter changes that did not match an expected mutagenesis signature? Default False.

        Warning! The exon feature GFF attributes must have transcript_id and exon_number or exon_ID, key-value pairs, and \
        features should be from 5' to 3', regardless of strand! For example, a transcript on the (-) will NOT be sorted \
        by coordinate if features are ordered 5' to 3'.
        """

        file_ext = fu.get_extension(gff)
        if not (re.match(ffu.GFF_FILETYPE, file_ext) or re.match(ffu.GTF_FILETYPE, file_ext)):
            raise NotImplementedError("Input file must be a GFF/GTF filetype.")

        self.gff = gff
        self.ref = ref
        self.outdir = outdir
        self.use_pickle = use_pickle
        self.make_pickle = make_pickle
        self.overwrite_pickle = overwrite_pickle
        self.filter_unexpected = filter_unexpected
        self.mut_sig = mut_sig

        self.output_dir = outdir
        if outdir is None:
            self.output_dir = tempfile.mkdtemp(suffix=__class__.__name__)

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

        self.pkl_filepath = os.path.join(self.output_dir, fu.replace_extension(os.path.basename(self.gff), "cds.pkl.gz"))
        pkl_exists = os.path.exists(self.pkl_filepath)

        if self.use_pickle and pkl_exists:
            logger.info("Loading pickled transcript CDS annotations.")
            with open(self.pkl_filepath, "rb") as info_pkl:
                self.cds_info = pickle.load(info_pkl)
        else:
            self.gff_bedtool = pybedtools.BedTool(self.gff)

            logger.info("Collecting transcript CDS annotations.")
            self.cds_info = self._get_cds_info()

        if (self.make_pickle and not pkl_exists) or self.overwrite_pickle:
            logger.info("Pickling transcript CDS annotations to %s for %s" % (self.pkl_filepath, self.gff))
            self.generate_pickle()

    @staticmethod
    def _add_cds_info(trx_id, transcript_cds_info, trx_seq, cds_seq):
        """Adds CDS information for a transcript.

        :param str trx_id: transcript ID
        :param dict transcript_cds_info: dict for storing information on each transcript
        :param str trx_seq: transcript sequence, from 5'-3'
        :param str cds_seq: CDS sequence, including the stop codon, from 5'-3'
        """

        # The information we need is the CDS sequence (including stop codon), the total transcript length
        # and the index into the transcript where the CDS starts (from these we know the location of UTRs)
        trx_len = len(trx_seq)

        # If we have a pseudogene or other non-translated transcript, we must set the CDS offsets to None
        cds_start_offset = None
        cds_stop_offset = None
        if len(cds_seq) > 0:
            cds_start_offset = trx_seq.index(cds_seq)  # 0-based start (BED)
            cds_stop_offset = cds_start_offset + len(cds_seq)  # 1-based stop (BED)

        transcript_cds_info[trx_id] = (trx_len, cds_start_offset, cds_stop_offset, trx_seq, cds_seq)

    def _get_cds_info(self):
        """Stores CDS information for each transcript."""

        transcript_cds_info = {}
        trx_id = None
        last_trx_id = None
        trx_seq = ""
        cds_seq = ""

        # feature will be a pybedtools.Interval object
        for feature in self.gff_bedtool:

            if self.CONTIGS_TO_OMIT.match(feature.chrom):
                continue

            # Filter the feature file of certain features so we have exon, CDS, and stop codon features only
            if self.FEATURES_TO_OMIT.match(feature.fields[ffu.GFF_FEATURE_TYPE_FIELD]):
                continue

            # Assuming we can identify the transcript, store the CDS data
            trx_id = feature.attrs[ffu.GFF_ATTR_TRANSCRIPT_ID] if ffu.GFF_ATTR_TRANSCRIPT_ID in feature.attrs else None
            if trx_id is not None and trx_id not in ffu.PYBEDTOOLS_NULL_CHARS:

                feature_strand = su.Strand(feature.strand)

                # Store the results of the last transcript and reset the sequences every time we switch between
                # blocks of features belonging to the same transcript
                if last_trx_id is not None and trx_id != last_trx_id:

                    self._add_cds_info(
                        trx_id=last_trx_id, transcript_cds_info=transcript_cds_info, trx_seq=trx_seq, cds_seq=cds_seq)

                    # Once we have added the information for the transcript in question, reset the sequence strings
                    trx_seq = ""
                    cds_seq = ""

                # Get the sequence of the entire transcript; we might normally just extract this from the transcriptome
                # reference FASTA, but we want this to be as extensible as possible (e.g. to store novel isoforms)
                if feature.fields[ffu.GFF_FEATURE_TYPE_FIELD] == self.EXON_ID:

                    exon_seq = su.extract_seq(
                        contig=str(feature.chrom), start=feature.start + 1, stop=feature.stop, ref=self.ref)
                    exon_seq = su.reverse_complement(exon_seq) if feature_strand == su.Strand.MINUS else exon_seq
                    trx_seq += exon_seq

                elif feature.fields[ffu.GFF_FEATURE_TYPE_FIELD] in self.VALID_CDS_FEATURES:

                    cds_exon_seq = su.extract_seq(
                        contig=str(feature.chrom), start=feature.start + 1, stop=feature.stop, ref=self.ref)
                    cds_exon_seq = su.reverse_complement(cds_exon_seq) if feature_strand == su.Strand.MINUS else cds_exon_seq
                    cds_seq += cds_exon_seq

                last_trx_id = trx_id

        # Make sure to add the information for the last transcript block
        if trx_id is not None:
            self._add_cds_info(
                trx_id=trx_id, transcript_cds_info=transcript_cds_info, trx_seq=trx_seq, cds_seq=cds_seq)

        return transcript_cds_info

    def _concat_multi_mut_info(self, mut_info_tuple):
        """Concatenates multi-position mutation changes into comma-delimited strings.

        :param collections.namedtuple mut_info_tuple: MUT_INFO_TUPLE with multi-position mutations in iterable
        :return collections.namedtuple MUT_INFO_TUPLE: MUT_INFO_TUPLE with multi-position mutations in str
        """

        wt_codons = self.MUT_INFO_DELIM.join(mut_info_tuple.wt_codons)
        mut_codons = self.MUT_INFO_DELIM.join(mut_info_tuple.mut_codons)
        wt_aas = self.MUT_INFO_DELIM.join(mut_info_tuple.wt_aas)
        mut_aas = self.MUT_INFO_DELIM.join(mut_info_tuple.mut_aas)
        aa_changes = self.MUT_INFO_DELIM.join(mut_info_tuple.aa_changes)
        aa_pos_list = [re.sub("[,p.A-Z*]", "", aa) for aa in mut_info_tuple.aa_changes]
        aa_positions = self.MUT_INFO_DELIM.join(aa_pos_list)
        matches_mut_sig = self.MUT_INFO_DELIM.join(list(map(str, mut_info_tuple.matches_mut_sig)))

        new_mut_info_tuple = MUT_INFO_TUPLE(
            location=mut_info_tuple.location, wt_codons=wt_codons, mut_codons=mut_codons, wt_aas=wt_aas, mut_aas=mut_aas,
            aa_changes=aa_changes, aa_positions=aa_positions, matches_mut_sig=matches_mut_sig)

        return new_mut_info_tuple

    @staticmethod
    def _get_pos_codon_dict(cds_seq):
        """Gets a mapping between index in the coding sequence and codon.

        :param str cds_seq: coding sequence
        :return dict: index in the CDS: CODON_TUPLE (codon, index of base in the codon)
        """

        cds_dict = {}
        wobble_positions = set(range(2, len(cds_seq), 3))
        codon_pos = 0

        for i in range(0, len(cds_seq)):
            if i % 3 == 0:
                codon_pos += 1
                curr_codon = cds_seq[i:i + 3]
                cds_dict[i] = CODON_TUPLE(codon_pos, curr_codon, 0)
                last_codon = curr_codon
            else:
                if i in wobble_positions:
                    cds_dict[i] = CODON_TUPLE(codon_pos, last_codon, 2)
                else:
                    cds_dict[i] = CODON_TUPLE(codon_pos, last_codon, 1)

        return cds_dict

    def _get_ins_downstream_remainder(self, trx_seq, pos, alt, curr_codon, base_index):
        """Gets the remaining codons and amino acids downstream of an insertion.

        :param str trx_seq: transcript sequence
        :param int pos: 1-based position of the variant within the transcript
        :param str alt: alternate bases
        :param str curr_codon: current reference codon
        :param int base_index: index of variant position in codon
        :return tuple: (remaining codons, remaining amino acids)
        """

        alt_plus_downstream = curr_codon[:base_index] + alt + trx_seq[pos:]

        alt_plus_downstream_codons = [
            alt_plus_downstream[i:i + 3] for i in range(0, len(alt_plus_downstream), 3)]

        remaining_codons = []
        remaining_aas = []

        for codon in alt_plus_downstream_codons:

            if len(codon) != 3:
                # Very unlikely this would happen before we break (if nonsense codon is not seen before the end of
                # the transcript), but need to handle
                remaining_codons.append(codon)
                remaining_aas.append(self.NONSTOP_CHAR)
                break

            remaining_codons.append(codon)
            remaining_aas.append(translate(codon))

            if codon in STOP_CODONS:
                break

        return tuple(remaining_codons), tuple(remaining_aas)

    def _get_del_downstream_remainder(self, trx_seq, pos, ref_len, curr_codon, base_index):
        """Gets the remaining codons and amino acids downstream of a deletion.

        :param str trx_seq: transcript sequence
        :param int pos: 1-based position of the variant within the transcript
        :param int ref_len: length of REF
        :param str curr_codon: current reference codon
        :param int base_index: 0-based index of variant position in codon
        :return tuple: (remaining codons, remaining amino acids)
        """

        alt_plus_downstream = curr_codon[:base_index + 1] + trx_seq[pos + ref_len - 1:]

        alt_plus_downstream_codons = [
            alt_plus_downstream[i:i + 3] for i in range(0, len(alt_plus_downstream), 3)]

        remaining_codons = []
        remaining_aas = []

        for codon in alt_plus_downstream_codons:

            if len(codon) != 3:
                # Very unlikely this would happen before we break (if nonsense codon is not seen before the end of
                # the transcript), but need to handle
                remaining_codons.append(codon)
                remaining_aas.append(self.NONSTOP_CHAR)
                break

            remaining_codons.append(codon)
            remaining_aas.append(translate(codon))

            if codon in STOP_CODONS:
                break

        return tuple(remaining_codons), tuple(remaining_aas)

    def _annotate_snp(self, ref_codon_dict, alt_codon_dict, start_index):
        """Gets annotations for a single-codon change.

        :param dict ref_codon_dict: index in the REF CDS: CODON_TUPLE (codon, index of base in the codon)
        :param dict alt_codon_dict: index in the ALT CDS: CODON_TUPLE (codon, index of base in the codon)
        :param int start_index: 0-based position of the variant in the CDS
        :return tuple: (ref_codons, alt_codons, ref_aas, alt_aas, aa_changes, matches_mut_sigs)
        """

        ref_codons = (ref_codon_dict[start_index].codon,)
        alt_codons = (alt_codon_dict[start_index].codon,)
        ref_aas = (translate(ref_codons[0]),)
        alt_aas = (translate(alt_codons[0]),)
        aa_changes = (HGVS_AA_FORMAT.format(ref_aas[0], ref_codon_dict[start_index].codon_pos, alt_aas[0]),)
        wt_wobble = ref_codons[0][2]
        mut_wobble = alt_codons[0][2]
        matches_mut_sigs = (False if self.mut_sig != self.MUT_SIG_ANY and wt_wobble != mut_wobble and
                                     mut_wobble in MUT_SIG_UNEXPECTED_WOBBLE_BPS[self.mut_sig] else True,)

        return ref_codons, alt_codons, ref_aas, alt_aas, aa_changes, matches_mut_sigs

    def _annotate_mnp(self, ref_codon_dict, alt_codon_dict, start_index, end_index):
        """Gets annotations for a multi-codon change.

        :param dict ref_codon_dict: index in the REF CDS: CODON_TUPLE (codon, index of base in the codon)
        :param dict alt_codon_dict: index in the ALT CDS: CODON_TUPLE (codon, index of base in the codon)
        :param int start_index: 0-based position of the first mismatch in the CDS
        :param int end_index: 0-based position of the last mismatch in the CDS
        :return tuple: (ref_codons, alt_codons, ref_aas, alt_aas, aa_changes, matches_mut_sigs)
        """

        ref_codons = []
        alt_codons = []
        ref_aas = []
        alt_aas = []
        aa_changes = []
        matches_mut_sigs = []

        for i in range(start_index, end_index + 1, 3):
            ref_codon = ref_codon_dict[i].codon
            alt_codon = alt_codon_dict[i].codon

            # Skip intervening codons that match
            if ref_codon == alt_codon:
                continue

            ref_aa = translate(ref_codon)
            alt_aa = translate(alt_codon)
            aa_change = HGVS_AA_FORMAT.format(ref_aa, ref_codon_dict[i].codon_pos, alt_aa)
            wt_wobble = ref_codon[2]
            mut_wobble = alt_codon[2]
            matches_mut_sig = False if self.mut_sig != self.MUT_SIG_ANY and wt_wobble != mut_wobble and \
                                       mut_wobble in MUT_SIG_UNEXPECTED_WOBBLE_BPS[self.mut_sig] else True

            ref_codons.append(ref_codon)
            alt_codons.append(alt_codon)
            ref_aas.append(ref_aa)
            alt_aas.append(alt_aa)
            aa_changes.append(aa_change)
            matches_mut_sigs.append(matches_mut_sig)

        return tuple(ref_codons), tuple(alt_codons), tuple(ref_aas), tuple(alt_aas), \
               tuple(aa_changes), tuple(matches_mut_sigs)

    def _annotate_ins(self, trx_seq, ref_codon_dict, start_index, pos, alt):
        """Gets annotations for an insertion.

        :param str trx_seq: transcript sequence
        :param dict ref_codon_dict: index in the REF CDS: CODON_TUPLE (codon, index of base in the codon)
        :param int start_index: 0-based position of the first mismatch in the CDS
        :param int pos: 1-based position of the variant within the transcript
        :param str alt: alternate bases
        :return tuple: (ref_codons, alt_codons, ref_aas, alt_aas, aa_changes, matches_mut_sigs)
        """

        alt_len = len(alt)

        if ref_codon_dict[start_index].base_index == 2 and ((alt_len - 1) % 3) == 0:
            # Insertion is in-frame and does not change current/downstream codons
            alt_ins = alt[1:]
            alt_list = [alt_ins[i:i + 3] for i in range(0, len(alt_ins), 3)]
            ref_codons = (ref_codon_dict[start_index].codon,)
            alt_codons = (ref_codons[0], *alt_list,)
            ref_aas = (translate(ref_codons[0]),)
            alt_aas = (ref_aas[0], *[translate(e) for e in alt_list],)
            aa_changes = (HGVS_AA_FORMAT.format(ref_aas[0], ref_codon_dict[start_index].codon_pos, ",".join(alt_aas)),)
            matches_mut_sigs = (True,)
        else:
            # Out-of-frame insertion or shuffling of downstream codons
            codon_remainder, aa_remainder = self._get_ins_downstream_remainder(
                trx_seq, pos, alt, ref_codon_dict[start_index].codon, ref_codon_dict[start_index].base_index)

            ref_codons = (ref_codon_dict[start_index].codon,)
            alt_codons = codon_remainder
            ref_aas = (translate(ref_codons[0]),)
            alt_aas = aa_remainder
            aa_changes = (HGVS_AA_FORMAT.format(ref_aas[0], ref_codon_dict[start_index].codon_pos, ",".join(alt_aas)),)
            matches_mut_sigs = (False,)

        return ref_codons, alt_codons, ref_aas, alt_aas, aa_changes, matches_mut_sigs

    def _annotate_del(self, trx_seq, ref_codon_dict, start_index, pos, ref):
        """Gets annotations for an insertion.

        :param str trx_seq: transcript sequence
        :param dict ref_codon_dict: index in the REF CDS: CODON_TUPLE (codon, index of base in the codon)
        :param int start_index: 0-based position of the first mismatch in the CDS
        :param int pos: 1-based position of the variant within the transcript
        :param str ref: reference bases
        :return tuple: (ref_codons, alt_codons, ref_aas, alt_aas, aa_changes, matches_mut_sigs)
        """

        ref_len = len(ref)

        if ref_codon_dict[start_index].base_index == 2 and ((ref_len - 1) % 3) == 0:
            # Deletion is in-frame and does not change current/downstream codons
            ref_del = ref[1:]
            ref_list = [ref_del[i:i + 3] for i in range(0, len(ref_del), 3)]
            ref_codons = (ref_codon_dict[start_index].codon, *ref_list,)
            alt_codons = (ref_codons[0],)
            ref_aas = (translate(ref_codons[0]), *[translate(e) for e in ref_list],)
            alt_aas = (ref_aas[0],)
            aa_changes = (HGVS_AA_FORMAT.format(",".join(ref_aas), ref_codon_dict[start_index].codon_pos, alt_aas[0]),)
            matches_mut_sigs = (True,)
        else:
            # Out-of-frame deletion or shuffling of downstream codons
            codon_remainder, aa_remainder = self._get_del_downstream_remainder(
                trx_seq, pos, ref_len, ref_codon_dict[start_index].codon, ref_codon_dict[start_index].base_index)

            ref_codons = (ref_codon_dict[start_index].codon,)
            alt_codons = codon_remainder
            ref_aas = (translate(ref_codons[0]),)
            alt_aas = aa_remainder
            aa_changes = (HGVS_AA_FORMAT.format(ref_aas[0], ref_codon_dict[start_index].codon_pos, ",".join(alt_aas)),)
            matches_mut_sigs = (False,)

        return ref_codons, alt_codons, ref_aas, alt_aas, aa_changes, matches_mut_sigs

    def _get_mut_info(self, trx_seq, wt_cds_seq, mut_cds_seq, pos, ref, alt, cds_start_offset):
        """Finds codon and amino acid changes given a variant POS, REF, and ALT fields.

        :param str trx_seq: transcript sequence
        :param str wt_cds_seq: WT CDS sequence
        :param str mut_cds_seq: Mutant CDS sequence
        :param int pos: 1-based position of the variant within the transcript
        :param str ref: reference bases
        :param str alt: alternate bases
        :param int cds_start_offset: transcript position of the CDS start nucleotide
        :return collections.namedtuple: MUT_INFO_TUPLE
        """

        ref_len = len(ref)
        alt_len = len(alt)
        start_index = pos - cds_start_offset - 1
        end_index = start_index if (ref_len == alt_len and ref_len == 1) else start_index + ref_len
        ref_codon_dict = self._get_pos_codon_dict(wt_cds_seq)

        if ref_len == alt_len:
            # Here we have a SNP or MNP so we can simply index into the codon lists
            alt_codon_dict = self._get_pos_codon_dict(mut_cds_seq)

            if start_index == end_index:
                # Here we have a SNP
                ref_codons, alt_codons, ref_aas, alt_aas, aa_changes, matches_mut_sigs = self._annotate_snp(
                    ref_codon_dict, alt_codon_dict, start_index)
            else:
                # Here we have a single-codon MNP or a multi-codon change (haplotype)
                ref_codons, alt_codons, ref_aas, alt_aas, aa_changes, matches_mut_sigs = self._annotate_mnp(
                    ref_codon_dict, alt_codon_dict, start_index, end_index)

        elif ref_len < alt_len:
            # Here we have an insertion
            ref_codons, alt_codons, ref_aas, alt_aas, aa_changes, matches_mut_sigs = self._annotate_ins(
                trx_seq, ref_codon_dict, start_index, pos, alt)
        else:
            # Here we have a deletion
            ref_codons, alt_codons, ref_aas, alt_aas, aa_changes, matches_mut_sigs = self._annotate_del(
                trx_seq, ref_codon_dict, start_index, pos, ref)

        # Collect the data into a mutation info tuple
        mut_info_tuple = MUT_INFO_TUPLE(
            location=self.CDS_ID, wt_codons=ref_codons, mut_codons=alt_codons,
            wt_aas=ref_aas, mut_aas=alt_aas, aa_changes=aa_changes, aa_positions=[],
            matches_mut_sig=matches_mut_sigs)

        concat_mut_info_tuple = self._concat_multi_mut_info(mut_info_tuple)

        return concat_mut_info_tuple

    def get_codon_and_aa_changes(self, trx_id, pos, ref, alt):
        """Given a transcript ID and variant position, determines stats associated with a transcriptomic variant.

        :param str trx_id: transcript ID; must correspond to the gff used in this object's init
        :param int pos: 1-based position of the variant within the transcript
        :param str ref: reference bases
        :param str alt: alternate bases
        :return collections.namedtuple: MUT_INFO_TUPLE for variant. If the variant is outside the CDS but within \
        the transcript bounds, then either '5_UTR' or '3_UTR' is returned for location; if the variant is outside the \
        transcript bounds, 'intergenic' is returned; if no CDS exists for the transcript, 'untranslated' is returned.
        :raises TranscriptNotFound: if the transcript ID was not found in the mapper object
        :raises RuntimeError: if the REF field of the variant does not match the reference sequence
        """

        if trx_id not in self.cds_info:
            raise TranscriptNotFound("No transcript %s found in the annotations." % trx_id)

        # cds_start_offset is 0-based offset of the first base of the start codon
        # cds_stop_offset is 0-based offset of the base _after_ the stop codon
        trx_len, cds_start_offset, cds_stop_offset, trx_seq, cds_seq = self.cds_info[trx_id]

        if pos > trx_len or pos < 1:
            var_key = vu.VARIANT_FORMAT.format(trx_id, pos, ref, alt)
            msg_str = "Variant %s position is not within transcript bounds." % var_key
            warnings.warn(msg_str, TranscriptomicCoordNotFound)
            logger.warning(msg_str)
            return MUT_INFO_TUPLE(location=self.INTERGENIC, **self.DEFAULT_KWARGS)

        # One more sanity-check
        zbased_pos = pos - 1
        expected_ref = trx_seq[zbased_pos:zbased_pos + len(ref)]
        if ref.upper() != expected_ref.upper():
            var_key = vu.VARIANT_FORMAT.format(trx_id, pos, ref, alt)
            msg_str = "REF %s of %s does not match the expected reference sequence of the transcript: %s" \
                      % (ref, var_key, expected_ref)
            logger.exception(msg_str)
            raise RuntimeError(msg_str)

        # Now that we have ensured our variant is valid, we can check for its placement
        if cds_start_offset is None and cds_stop_offset is None:
            return MUT_INFO_TUPLE(location=self.UNTRANSLATED, **self.DEFAULT_KWARGS)
        elif zbased_pos < cds_start_offset:
            return MUT_INFO_TUPLE(location=self.FIVEPRIME_UTR, **self.DEFAULT_KWARGS)
        elif zbased_pos >= cds_stop_offset:
            return MUT_INFO_TUPLE(location=self.THREEPRIME_UTR, **self.DEFAULT_KWARGS)
        else:
            # We have a valid variant in the coding region; annotate AA changes
            cds_len = len(cds_seq)
            ref_len = len(ref)

            # Create the ALT within the CDS sequence
            cds_pos = zbased_pos - cds_start_offset
            mut_cds_seq = cds_seq[:cds_pos] + alt + cds_seq[cds_pos + ref_len:]

            # For SNPs and MNPs, ensure the WT and mutant CDS sequences are the same length
            # This will not call MNPs that span the CDS/3' UTR junction; this is a potential bug, but
            # nonstop variants are unlikely to be generated in most mutagenesis strategies
            # We should not slice back to the same length of the WT CDS for insertions; _get_mut_info will handle InDels
            if ref_len == len(alt):
                mut_cds_seq = mut_cds_seq[:cds_len]

            # Extract the information about the consequences of the variant change
            mut_info_tuple = self._get_mut_info(trx_seq, cds_seq, mut_cds_seq, pos, ref, alt, cds_start_offset)

            return mut_info_tuple

    def generate_pickle(self):
        """Makes a pickle of the object."""

        with open(self.pkl_filepath, "wb") as info_pkl:
            pickle.dump(obj=self.cds_info, file=info_pkl)

        gzip.GzipFile(self.pkl_filepath)


def translate(codon):
    """Translates a codon.

    :param str codon: codon string
    :return str: shorthand amino acid string
    :raises RuntimeError: if the codon is not in the recognized set of 64 codons
    """

    codon_upper = codon.upper()

    if codon_upper not in CODON_AA_DICT:
        raise RuntimeError("Unrecognized codon %s." % codon)

    aa = CODON_AA_DICT[codon_upper]
    return aa
