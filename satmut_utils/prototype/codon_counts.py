#!/usr/bin/env/python
"""Objects for counting codons in an ORF."""

import collections
import logging
import tempfile

from analysis.coordinate_mapper import AminoAcidMapper, TranscriptNotFound
import core_utils.file_utils as fu
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
    "MUT_INFO_TUPLE", "location, wt_codons, mut_codons, wt_aas, mut_aas, aa_changes, aa_positions, matches_nnk")

tempfile.tempdir = os.getenv("/tmp")
_logger = logging.getLogger(__name__)


class CodonTrnaCounts(AminoAcidMapper):
    """Class for enumerating codons and tRNA counts for a CDS."""

    HG19_CODON_TRNA_COUNTS_FN = "nnk_codon_permutations/hg19-tRNAs-confidence-set.codons.counts.txt"
    HG19_CODON_TRNA_COUNTS = collections.defaultdict(int)

    with open(HG19_CODON_TRNA_COUNTS_FN, "r") as counts_fh:
        for line in counts_fh:
            k, v = line.split(fu.FILE_DELIM)
            HG19_CODON_TRNA_COUNTS[k] = int(v)

    CANDIDATE_PAUSE_CODONS = ("ATC", "GCC", "ACC", "CCC", "TCC", "GAT", "CTC", "CGC", "TGT", "GTC", "TTT", "AGT",
                              "AAT", "TAT", "GGT", "CTA", "CAT")

    def __init__(self, gff, ref, use_pickle=True, make_pickle=True, overwrite_pickle=False,
                 mut_sig=AminoAcidMapper.MUT_SIG_ANY, filter_unexpected=False):
        r"""Ctor for CodonTrnaCounts.

        :param str gff: GFF/GTF to create a mapper for; must have "transcript_id" and "CDS", "stop_codon" features
        :param str ref: reference FASTA corresponding to GFF features
        :param bool use_pickle: Use a previously generated pickle for annotations if one exists? Default True. Otherwise \
        re-create the data structure.
        :param bool make_pickle: Should a new pickle be made? Default True.
        :param bool overwrite_pickle: if make_pickle=True and use_pickle=True, should the existing pickle be \
        overwritten? Default False. This is useful if one makes dynamic edits to an existing GFF. Otherwise, the older \
        pickle will be used.
        :param str mut_sig: mutagenesis signature- one of {NNN, NNK, NNS}. Default NNN.
        :param bool filter_unexpected: Filter codon/AA changes that did not match an expected mut_sig? Default False.

        Warning! The exon feature GFF attributes must have transcript_id and exon_number or exon_ID, key-value pairs, and \
        features should be from 5' to 3', regardless of strand! For example, a transcript on the (-) will NOT be sorted \
        by coordinate if features are ordered 5' to 3'.
        """

        super(CodonTrnaCounts, self).__init__(
            gff, ref, use_pickle, make_pickle, overwrite_pickle, mut_sig, filter_unexpected)

    @classmethod
    def count_codons(cls, orf_seq, outfile=None):
        """Gets counts and tRNA copy numbers for codons in orf_seq.

        :param str orf_seq: open reading frame sequence
        :param str | None outfile: Optional filename to write results to
        :return dict: dictionary enumerating codons and respective tRNA copy numbers
        """

        orf_codons = [orf_seq[i:i + 3] for i in range(0, len(orf_seq), 3)]
        orf_codons_counter = collections.Counter(orf_codons)

        # Annotate with tRNA copy numbers
        orf_codons_annot = {k: (v, cls.HG19_CODON_TRNA_COUNTS[k]) for k, v in orf_codons_counter.items()}

        out_file = outfile
        if outfile is None:
            out_file = tempfile.NamedTemporaryFile(mode="w", suffix=".codon.counts.txt", delete=False)

        with open(out_file, "w") as out_fh:
            out_fh.write(fu.FILE_DELIM.join(["Codon", "ORF_count", "tRNA_count"]) + fu.FILE_NEWLINE)

            for k, (v1, v2) in orf_codons_annot.items():
                out_fh.write(fu.FILE_DELIM.join(list(map(str, [k, v1, v2]))) + fu.FILE_NEWLINE)

        return orf_codons_annot

    def get_codon_positions(self, trx_id, query_codons=CANDIDATE_PAUSE_CODONS):
        """Gets the transcript positions of codons.

        :param str trx_id: transcript ID to query
        :param tuple query_codons: codons to report positions for
        :return collections.defaultdict: dictionary keyed by codon, valued by list of 1-based start positions of codon
        """

        if trx_id not in self.cds_info:
            raise TranscriptNotFound(
                "No transcript %s found in the annotations. Do you need to recreate the data structure?" % trx_id)

        # cds_start_offset is 0-based offset of the first base of the start codon
        # cds_stop_offset is 0-based offset of the base _after_ the stop codon
        trx_len, cds_start_offset, cds_stop_offset, trx_seq, cds_seq = self.cds_info[trx_id]
        orf_codons = [cds_seq[i:i + 3] for i in range(0, len(cds_seq), 3)]

        query_codons = set(query_codons)
        codon_trx_positions = collections.defaultdict(list)

        for i, codon in enumerate(orf_codons):
            if codon in query_codons:
                codon_trx_positions[codon].append(cds_start_offset + (i * 3) + 1)

        return codon_trx_positions

    def get_codon_position_bed(self, trx_id, query_codons=CANDIDATE_PAUSE_CODONS, outfile=None):
        """Generates a BED file of positions with candidate pause codons.

        :param str trx_id: transcript ID to query
        :param tuple query_codons: codons to report positions for
        :param str | None outfile: Optional BED filename to write results to
        """

        # Transcript positions of the query codons
        codon_positions = self.get_codon_positions(trx_id, query_codons)

        # All codons in the CDS and their respective tRNA copy numbers
        cds_codon_counts = self.count_codons(self.cds_info[trx_id][self.CDS_INFO_CDS_SEQ_INDEX])

        out_bed = outfile
        if outfile is None:
            out_bed = tempfile.NamedTemporaryFile(mode="w", suffix="candidate.pause.codons.bed", delete=False).name

        with open(out_bed, "w") as out_bed_fh:
            for codon, codon_positions in codon_positions.items():
                for position in codon_positions:

                    codon_counts, trna_copy_number = cds_codon_counts[codon]

                    bed_rec = list(map(str, [
                        trx_id, position - 1, position + 2, ":".join([codon, str(trna_copy_number)])
                    ]))

                    out_bed_fh.write(fu.FILE_DELIM.join(bed_rec) + fu.FILE_NEWLINE)
