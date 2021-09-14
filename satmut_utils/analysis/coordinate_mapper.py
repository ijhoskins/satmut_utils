#!/usr/bin/env/python
"""Objects for mapping between genomic, transcriptomic, and protein coordinates."""

import aenum
import Bio.Seq
import collections
import gzip
import logging
import pickle
import pybedtools

import re
import tempfile
import warnings


import analysis.seq_utils as su
import prototype.variant_generator as vg
import core_utils.feature_file_utils as ffu
import core_utils.file_utils as fu
import core_utils.vcf_utils as vu
from definitions import *

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

EXON_COORDS_TUPLE = collections.namedtuple("EXON_COORDS_TUPLE", "exon_id, contig, start, stop, exon_len, strand")

MUT_INFO_TUPLE = collections.namedtuple(
    "MUT_INFO_TUPLE", "location, wt_codons, mut_codons, wt_aas, mut_aas, aa_changes, aa_positions, matches_mut_sig")

tempfile.tempdir = os.getenv("/tmp")
_logger = logging.getLogger(__name__)


class TranscriptNotFound(Exception):
    """Exception for when a transcript ID was not found in the CoordinateMapper."""
    pass


class GenomicCoordNotFound(Exception):
    """Exception for when a local index to a transcript is outside of exonic bounds."""
    pass


class TranscriptomicCoordNotFound(Warning):
    """Warning for when a genomic coordinate is outside of transcript bounds."""
    pass


class MapperBase(object):
    """Base class for GTF/GFF related data parsing."""

    # Some GFF/GTF type fields; some may not correspond to actual fields in your file

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


class CoordinateMapper(MapperBase):
    r"""Class for mapping transcriptomic coordinates to genomic coordinates.

    Warning! The exon feature GFF attributes must have transcript_id, and exon_number/exon_ID key-value pairs, and \
    features should be from 5' to 3', regardless of strand! For example, a transcript on the (-) will NOT be sorted \
    by coordinate if features are ordered 5' to 3'. This is typically the case for GFF files downloaded from NCBI and \
    other sources, but any coordinate sorting of the GFF prior to this object's init would lead to incorrect genomic \
    coordinates for (-) strand features.

    Note this object can map arbitrary metafeatures (e.g. transcript fusions, novel isoforms) that would otherwise be \
    impossible to map to genomic space using a liftOver approach.
    """

    # Depending on the source of the GFF/GTF, this may need to be modified to avoid inclusion of other metafeatures
    # While these other features are useful for various reasons, for the CoordinateMapper all we care about is the \
    # exon features (both UTR and CDS exons), so we can map a local transcript index from the transcriptome alignments \
    # to genomic space
    FEATURES_TO_OMIT = re.compile("|".join(MapperBase.FEATURES_TO_OMIT +
        [MapperBase.GENE_ID, MapperBase.TRX_ID, MapperBase.MRNA_ID, MapperBase.CDS_ID,
         MapperBase.UTR_ID, MapperBase.START_CODON_ID, MapperBase.STOP_CODON_ID]))

    def __init__(self, gff, use_pickle=True, make_pickle=True, overwrite_pickle=False):
        r"""Constructor for CoordinateMapper.

        :param str gff: GFF/GTF to create a mapper for; must have "transcript_id" and "exon_number"/"exon_id" attributes.
        :param bool use_pickle: Use a previously generated pickle for annotations if one exists? Default True. Otherwise \
        re-create the data structure.
        :param bool make_pickle: Should a new pickle be made? Default True.
        :param bool overwrite_pickle: if make_pickle=True and use_pickle=True, should the existing pickle be \
        overwritten? Default False. This is useful if one makes dynamic edits to an existing GFF. Otherwise, the older \
        pickle will be used.

        Warning! The exon feature GFF attributes must have transcript_id and exon_number or exon_ID, key-value pairs, and
        features should be from 5' to 3', regardless of strand! For example, a transcript on the (-) will NOT be sorted
        by coordinate if features are ordered 5' to 3'.
        """

        file_ext = fu.get_extension(gff)
        if not (re.match(ffu.GFF_FILETYPE, file_ext) or re.match(ffu.GTF_FILETYPE, file_ext)):
            raise NotImplementedError("Input file must be a GFF/GTF filetype.")

        self.gff = gff
        self.use_pickle = use_pickle
        self.make_pickle = make_pickle
        self.overwrite_pickle = overwrite_pickle

        self.base_dir = os.path.dirname(os.path.abspath(self.gff))
        self.pkl_filepath = os.path.join(self.base_dir, fu.replace_extension(os.path.basename(self.gff), "exon.pkl.gz"))
        pkl_exists = os.path.exists(self.pkl_filepath)

        if self.use_pickle and pkl_exists:
            _logger.info("Loading pickled transcript exon annotations.")
            with open(self.pkl_filepath, "rb") as info_pkl:
                self.exon_info = pickle.load(info_pkl)
        else:
            self.gff_bedtool = pybedtools.BedTool(self.gff)

            _logger.info("Collecting transcript exon annotations...")
            self.exon_info = self._get_exon_info()

        if (self.make_pickle and not pkl_exists) or self.overwrite_pickle:
            _logger.info("Pickling transcript exon annotations to %s for %s." % (self.pkl_filepath, self.gff))
            self.generate_pickle()

    def _get_exon_info(self):
        """Stores exon information for each transcript."""

        # Attempt to group exons by transcript; iter() method of BedTool returns Interval objects
        grouped_features = collections.defaultdict(list)

        for i, feature in enumerate(self.gff_bedtool):

            if self.CONTIGS_TO_OMIT.match(feature.chrom):
                continue

            # Filter the feature file of metafeatures (transcripts, genes, etc.) so we have just exons
            if self.FEATURES_TO_OMIT.match(feature.fields[ffu.GFF_FEATURE_TYPE_FIELD]):
                continue

            # Assuming we can identify the transcript and exon number of each exon feature, group
            trx_id = feature.attrs[ffu.GFF_ATTR_TRANSCRIPT_ID] if ffu.GFF_ATTR_TRANSCRIPT_ID in feature.attrs else None

            # Attempt to get a unique exon number or ID
            exon_id = feature.attrs[ffu.GFF_ATTR_EXON_ID] if ffu.GFF_ATTR_EXON_ID in feature.attrs else None
            if exon_id is None:
                # Last attempt to find a unique identifier for the exon
                exon_id = feature.attrs[ffu.GFF_ATTR_ID_ID] if ffu.GFF_ATTR_ID_ID in feature.attrs else None

            # We will not store information if we do not have a valid transcript or exon ID for the currently-
            # iterated feature
            if trx_id is not None and trx_id not in ffu.PYBEDTOOLS_NULL_CHARS and \
                    exon_id is not None and exon_id not in ffu.PYBEDTOOLS_NULL_CHARS:

                feature_strand = su.Strand(feature.strand)

                # Construct the genomic range of the exon coordinates; pybedtools convention is that start is 0-based
                # In keeping with GFF/GTF format, convert start to 1-based
                # Note start and stop here are strand aware
                if feature_strand == su.Strand.PLUS:
                    exon_juncs = [feature.start + 1, feature.stop]
                elif feature_strand == su.Strand.MINUS:
                    exon_juncs = [feature.stop, feature.start + 1]
                else:
                    raise RuntimeError("Found an invalid strand for feature on line %i" % (i + 1))

                # Remove version numbers on transcripts/chromosome names, as some annotation files do not contain one
                grouped_features[trx_id].append(
                    EXON_COORDS_TUPLE(
                        exon_id, fu.remove_extension(feature.chrom),
                        exon_juncs[0], exon_juncs[1], len(feature), feature_strand
                    )
                )

        return grouped_features

    @staticmethod
    def _get_genomic_coordinate(transcript_info, local_transcript_index):
        r"""Gets the genomic coordinate from a local transcript index.

        :param list transcript_info: transcript exon data
        :param int local_transcript_index: 0-based local transcript index
        :return tuple: (str, str, int) (exon ID, contig, 1-based genomic coordinate)
        :raises analysis.coordinate_mapper.GenomicCoordNotFound: if the local index is out of bounds

        Note: this logic is more complicated than just indexing into a list of the genomic positions for each \
        transcript, but it significantly decreases the memory required to store the annotations, and further allows \
        us to know what exons each local index is in.
        """

        cumulative_transcript_len = 0
        last_ect = transcript_info[0]

        for i, ect in enumerate(transcript_info):

            # Skip to the exon that contains the local index
            if local_transcript_index >= cumulative_transcript_len:
                cumulative_transcript_len += ect.exon_len
                last_ect = ect
                if i + 1 != len(transcript_info):
                    continue

            # Here we are at the exon containing the local index
            assoc_exon_id = last_ect.exon_id
            # Depending on the nature of the GFF reference, extract a str value
            genomic_contig = last_ect.contig.value if isinstance(last_ect.contig, aenum.Enum) else last_ect.contig
            cumulative_transcript_len = cumulative_transcript_len - last_ect.exon_len
            exon_offset = local_transcript_index - cumulative_transcript_len

            if last_ect.strand == su.Strand.PLUS:
                # Increment/decrement stop by one for python range behavior vs 1-based system
                genomic_coords = list(range(last_ect.start, last_ect.stop + 1, 1))
            else:
                genomic_coords = list(range(last_ect.start, last_ect.stop - 1, -1))

            # Now index into the genomic coordinates at the exon
            genomic_coord = genomic_coords[exon_offset]
            break
        else:
            raise GenomicCoordNotFound("No genomic coordinate for local index %i" % local_transcript_index)

        return assoc_exon_id, genomic_contig, genomic_coord

    def get_genomic_association(self, transcript_id, local_transcript_indices):
        """From a local transcript index, return the exon number and genomic coordinate.

        :param str transcript_id: transcript ID
        :param tuple local_transcript_indices: 0-based local transcript indices
        :return list: [(str, str, int), ...] (exon ID, contig, 1-based genomic coordinate)
        :raises analysis.coordinate_mapper.TranscriptNotFound: if the transcript was not found in the structure
        """

        if transcript_id in self.exon_info:
            transcript_info = self.exon_info[transcript_id]

            genomic_coords = [self._get_genomic_coordinate(transcript_info, local_index) for local_index in local_transcript_indices]
            return genomic_coords

        raise TranscriptNotFound(
            "No transcript %s found in the annotations. Do you need to recreate the data structure?" % transcript_id)

    def get_transcript_strand(self, transcript_id):
        """Gets the genomic strand of the transcript.

        :param str transcript_id: transcript ID
        :return analysis.seq_utils.Strand: strand Enum
        :raises analysis.coordinate_mapper.TranscriptNotFound: if the transcript was not found in the structure
        """

        # Note this returns the strand of the first exon; in rare cases where a metafeature contains exons on different
        # strands (e.g. transcript fusion), no warning will be issued
        if transcript_id in self.exon_info:
            transcript_info = self.exon_info[transcript_id]
            return transcript_info[0].strand

        raise TranscriptNotFound(
            "No transcript %s found in the annotations. Do you need to recreate the data structure?" % transcript_id)

    def get_transcript_length(self, transcript_id):
        """Gets the genomic strand of the transcript.

        :param str transcript_id: transcript ID
        :return int: length of the transcript
        :raises analysis.coordinate_mapper.TranscriptNotFound: if the transcript was not found in the structure
        """

        if transcript_id in self.exon_info:
            transcript_info = self.exon_info[transcript_id]

            transcript_length = 0
            for exon in transcript_info:
                transcript_length += exon.exon_len

            return transcript_length

        raise TranscriptNotFound(
            "No transcript %s found in the annotations. Do you need to recreate the data structure?" % transcript_id)

    def generate_pickle(self):
        """Makes a pickle of the object; for the human transcriptome, the pickle is roughly 50 Mb."""

        with open(self.pkl_filepath, "wb") as info_pkl:
            pickle.dump(obj=self.exon_info, file=info_pkl)

        gzip.GzipFile(self.pkl_filepath)


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
    MUT_SIG_ANY = "NNN"
    CDS_INFO_TRX_LEN_INDEX = 0
    CDS_INFO_CDS_START_INDEX = 1
    CDS_INFO_CDS_STOP_INDEX = 2
    CDS_INFO_TRX_SEQ_INDEX = 3
    CDS_INFO_CDS_SEQ_INDEX = 4
    DEFAULT_ARGS = MUT_INFO_TUPLE._fields[1:]
    DEFAULT_KWARGS = dict(zip(DEFAULT_ARGS, [su.R_COMPAT_NA] * len(DEFAULT_ARGS)))

    def __init__(self, gff, ref, use_pickle=True, make_pickle=True, overwrite_pickle=False,
                 mut_sig=MUT_SIG_ANY, filter_unexpected=False):
        r"""Constructor for AminoAcidMapper.

        :param str gff: GFF/GTF to create a mapper for; must have "transcript_id" and "CDS", "stop_codon" features.
        :param str ref: reference FASTA corresponding to GFF features
        :param bool use_pickle: Use a previously generated pickle for annotations if one exists? Default True. Otherwise \
        re-create the data structure.
        :param bool make_pickle: Should a new pickle be made? Default True.
        :param bool overwrite_pickle: if make_pickle=True and use_pickle=True, should the existing pickle be \
        overwritten? Default False. This is useful if one makes dynamic edits to an existing GFF. Otherwise, the older \
        pickle will be used.
        :param str mut_sig: mutagenesis signature- one of {NNN, NNK, NNS}. Default NNN.
        :param bool filter_unexpected: Filter codon/AA changes that did not match an expected NNK signature? This is a \
        chemistry-specific filter to enable filtering of false positive calls from POPcode-mutagenized data. However, \
        even if a variant does not match this signature, it could still be "real" in the sense that the variant was \
        really in the data due to a mutagensis primer synthesis error (despite it being unintended by design). Default False.

        Warning! The exon feature GFF attributes must have transcript_id and exon_number or exon_ID, key-value pairs, and \
        features should be from 5' to 3', regardless of strand! For example, a transcript on the (-) will NOT be sorted \
        by coordinate if features are ordered 5' to 3'.
        """

        file_ext = fu.get_extension(gff)
        if not (re.match(ffu.GFF_FILETYPE, file_ext) or re.match(ffu.GTF_FILETYPE, file_ext)):
            raise NotImplementedError("Input file must be a GFF/GTF filetype.")

        self.gff = gff
        self.ref = ref
        self.use_pickle = use_pickle
        self.make_pickle = make_pickle
        self.overwrite_pickle = overwrite_pickle
        self.filter_unexpected = filter_unexpected
        self.mut_sig = mut_sig

        self.base_dir = os.path.dirname(os.path.abspath(self.gff))
        self.pkl_filepath = os.path.join(self.base_dir, fu.replace_extension(os.path.basename(self.gff), "cds.pkl.gz"))
        pkl_exists = os.path.exists(self.pkl_filepath)

        if self.use_pickle and pkl_exists:
            _logger.info("Loading pickled transcript CDS annotations.")
            with open(self.pkl_filepath, "rb") as info_pkl:
                self.cds_info = pickle.load(info_pkl)
        else:
            self.gff_bedtool = pybedtools.BedTool(self.gff)

            _logger.info("Collecting transcript CDS annotations...")
            self.cds_info = self._get_cds_info()

        if (self.make_pickle and not pkl_exists) or self.overwrite_pickle:
            _logger.info("Pickling transcript CDS annotations to %s for %s" % (self.pkl_filepath, self.gff))
            self.generate_pickle()

    @staticmethod
    def _add_cds_info(trx_id, transcript_cds_info, trx_seq, cds_seq):
        """Adds CDS information for a transcript.

        :param str trx_id: transcript ID
        :param collections.defaultdict transcript_cds_info: dict for storing information on each transcript
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

        transcript_cds_info = collections.defaultdict(tuple)
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
                # reference FASTA, but we want this to be as extensible as possible (e.g. to map transcript fusions)
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

        :param collections.namedtuple mut_info_tuple: MUT_INFO_TUPLE with multi-position mutations in list
        :return collections.namedtuple MUT_INFO_TUPLE: MUT_INFO_TUPLE with multi-position mutations in str
        """

        wt_codons = self.MUT_INFO_DELIM.join(mut_info_tuple.wt_codons)
        mut_codons = self.MUT_INFO_DELIM.join(mut_info_tuple.mut_codons)
        wt_aas = self.MUT_INFO_DELIM.join(mut_info_tuple.wt_aas)
        mut_aas = self.MUT_INFO_DELIM.join(mut_info_tuple.mut_aas)
        aa_changes = self.MUT_INFO_DELIM.join(mut_info_tuple.aa_changes)
        aa_pos_list = [re.sub("[p.A-Z*]", "", aa) for aa in mut_info_tuple.aa_changes]
        aa_positions = self.MUT_INFO_DELIM.join(aa_pos_list)

        new_mut_info_tuple = MUT_INFO_TUPLE(
            location=mut_info_tuple.location, wt_codons=wt_codons, mut_codons=mut_codons, wt_aas=wt_aas, mut_aas=mut_aas,
            aa_changes=aa_changes, aa_positions=aa_positions, matches_mut_sig=mut_info_tuple.matches_mut_sig)

        return new_mut_info_tuple

    def _get_mut_info(self, wt_cds_seq, mut_cds_seq):
        """Finds codon and amino acid changes.

        :param str wt_cds_seq: WT CDS sequence
        :param str mut_cds_seq: Mutant CDS sequence
        :return collections.namedtuple: MUT_INFO_TUPLE
        """

        # If we weren't interested in reporting silent mutations, we might just translate the whole sequences then
        # compare AAs; but since we want this output for comparison to the Fritz Roth lab caller, translate each codon
        # Multi-base MNPs (haplotypes) make optimization logic a little tricky, as the MNP may span multiple codons
        # For now just iterate over all the codons/AAs
        codon_comparitors = zip(
            [wt_cds_seq[i:i + 3] for i in range(0, len(wt_cds_seq), 3)],
            [mut_cds_seq[i:i + 3] for i in range(0, len(mut_cds_seq), 3)]
        )

        # In the future this may need to be modified to handle InDels
        # Checking codons is needed to call synonymous AA changes
        wt_codons = []
        mut_codons = []
        wt_aas = []
        mut_aas = []
        mut_pos = set()
        aa_changes = []
        full_wt_aas = []
        full_mut_aas = []
        matches_mut_sig = True

        for i, (wt_codon, mut_codon) in enumerate(codon_comparitors):

            wt_aa = str(Bio.Seq.Seq(wt_codon).translate())
            mut_aa = str(Bio.Seq.Seq(mut_codon).translate())

            full_wt_aas.append(wt_aa)
            full_mut_aas.append(mut_aa)

            if wt_codon != mut_codon:

                wt_codons.append(wt_codon)
                mut_codons.append(mut_codon)

                # Either mark or filter changes that are unexpected by looking at the wobble bp
                wt_wobble = wt_codon[2]
                mut_wobble = mut_codon[2]

                # Determine if the mutation matches the mutagenesis signature
                if self.mut_sig != self.MUT_SIG_ANY and wt_wobble != mut_wobble and \
                        mut_wobble in vg.MUT_SIG_UNEXPECTED_WOBBLE_BPS[self.mut_sig]:
                    matches_mut_sig = False

                    if self.filter_unexpected:
                        matches_mut_sig = True  # reset for other AAs downstream
                        continue

                wt_aas.append(wt_aa)
                mut_aas.append(mut_aa)
                mut_pos.add(i + 1)
                aa_changes.append(vg.HGVS_AA_FORMAT.format(wt_aa, i + 1, mut_aa))

        # Collect the data into a mutation info tuple
        mut_info_tuple = MUT_INFO_TUPLE(
            location=self.CDS_ID, wt_codons=wt_codons, mut_codons=mut_codons,
            wt_aas=wt_aas, mut_aas=mut_aas, aa_changes=aa_changes, aa_positions=[],
            matches_mut_sig=matches_mut_sig)

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
        """

        if trx_id not in self.cds_info:
            raise TranscriptNotFound(
                "No transcript %s found in the annotations. Do you need to recreate the data structure?" % trx_id)

        # cds_start_offset is 0-based offset of the first base of the start codon
        # cds_stop_offset is 0-based offset of the base _after_ the stop codon
        trx_len, cds_start_offset, cds_stop_offset, trx_seq, cds_seq = self.cds_info[trx_id]

        if pos > trx_len or pos < 1:
            var_key = vu.VARIANT_FORMAT.format(trx_id, pos, ref, alt)
            warnings.warn("Variant %s position is not within transcript bounds." % var_key, TranscriptomicCoordNotFound)
            return MUT_INFO_TUPLE(location=self.INTERGENIC, **self.DEFAULT_KWARGS)

        # One more sanity-check
        zbased_pos = pos - 1
        expected_ref = trx_seq[zbased_pos:zbased_pos + len(ref)]
        if ref.upper() != expected_ref.upper():
            raise RuntimeError(
                "Provided REF field does not match the expected reference of the stored transcript: %s" % expected_ref)

        # Now that we have ensured our variant is valid, we can check for its placement
        if cds_start_offset is None and cds_stop_offset is None:
            return MUT_INFO_TUPLE(location=self.UNTRANSLATED, **self.DEFAULT_KWARGS)
        elif zbased_pos < cds_start_offset:
            return MUT_INFO_TUPLE(location=self.FIVEPRIME_UTR, **self.DEFAULT_KWARGS)
        elif zbased_pos >= cds_stop_offset:
            return MUT_INFO_TUPLE(location=self.THREEPRIME_UTR, **self.DEFAULT_KWARGS)
        else:
            # To determine the AA change, translate the mutant CDS sequence and find the differences
            cds_len = len(cds_seq)
            ref_len = len(ref)

            # Create the ALT within the CDS sequence
            cds_pos = zbased_pos - cds_start_offset
            mut_cds_seq = cds_seq[:cds_pos] + alt + cds_seq[cds_pos + ref_len:]
            mut_cds_seq = mut_cds_seq[:cds_len]

            # Extract the information about the consequences of the variant change
            mut_info_tuple = self._get_mut_info(wt_cds_seq=cds_seq, mut_cds_seq=mut_cds_seq)

            return mut_info_tuple

    def generate_pickle(self):
        """Makes a pickle of the object."""

        with open(self.pkl_filepath, "wb") as info_pkl:
            pickle.dump(obj=self.cds_info, file=info_pkl)

        gzip.GzipFile(self.pkl_filepath)

