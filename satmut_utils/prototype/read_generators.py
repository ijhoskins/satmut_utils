#!/usr/bin/env python3
"""Objects for generating in silico NGS reads."""

import abc
import collections
import numpy
import os
import pybedtools
import random
import re
import tempfile

from analysis.seq_utils import fasta_to_fastq, extract_seq, reverse_complement, FASTA_HEADER_CHAR, \
    DEFAULT_READ_LEN, DEFAULT_FRAG_LEN, MIN_FRAG_LENGTH, SD_FROM_MEAN_FACTOR, COORD_FORMAT, Strand, HumanContig

from core_utils.feature_file_utils import slop_features, get_genome_file, DEFAULT_BP_SLOP, PYBEDTOOLS_NULL_CHARS, \
    GTF_FILETYPE, GFF_FILETYPE, GFF_ATTR_GENE_ID, GFF_ATTR_TRANSCRIPT_ID, GFF_FEATURE_TYPE_FIELD

from core_utils.file_utils import safe_remove, FILE_NEWLINE, replace_extension
from core_utils.string_utils import make_unique_ids, is_number

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.2"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

tempfile.tempdir = os.getenv("SCRATCH", "/tmp")

DEFAULT_NREADS = 30
MEAN_EXON_LEN = 170


class InsufficientFragLengthException(Exception):
    """Exception raised when a desired fragment length is too small to generate reads."""
    pass


class ReadGenerator(object):
    """Class for generation of in silico sequencing reads."""

    __metaclass__ = abc.ABCMeta

    DEFAULT_SNPS = False
    DEFAULT_INDELS = False

    def __init__(self, feature_file, ref, paired=True, rna=False,
                 read_length=DEFAULT_READ_LEN, frag_length=DEFAULT_FRAG_LEN, slop_length=DEFAULT_BP_SLOP,
                 make_amplicons=False):
        r"""Ctor for ReadGenerator.

        :param str feature_file: BED, GFF, GTF file containing features to make reads for; e.g. target BED
        :param str ref: genomic reference FASTA corresponding to the features
        :param bool paired: should paired-end reads be created? Default True.
        :param bool rna: should RNA reads be created? Assumes features are exons and metafeatures are transcripts. \
        Default False, make DNA reads.
        :param int read_length: length of reads to generate
        :param int frag_length: average fragment length to simulate
        :param int slop_length: number of bases to slop around targets for DNA read generation
        :param bool make_amplicons: should reads be generated from the ends of the target regions? Default False. \
        Make uniform random ends.
        """

        self.feature_file = feature_file
        self.ff_bedtool = pybedtools.BedTool(feature_file)
        self.ref = ref
        self.paired = paired
        self.rna = rna
        self.read_length = read_length
        self.frag_length = frag_length
        if self.frag_length < MIN_FRAG_LENGTH:
            raise InsufficientFragLengthException(
                "Desired fragment length will produce reads that are too small to map. "
                "Please choose fragment length >= %i bp." % MIN_FRAG_LENGTH)

        self.slop_length = slop_length
        self.make_amplicons = make_amplicons

    @abc.abstractmethod
    def make_reads(self, nreads=DEFAULT_NREADS):
        """Enforced method to override for making reads.

        :param int nreads: mean number of reads/pairs to generate
        """
        pass

    @staticmethod
    def remove_version_extension(id):
        """Removes the version number from a RefSeq ID.

        :param str id: Contig or RefSeq ID, e.g. NM_014453.3
        :return str: pared-down ID, e.g. NM_014453
        """

        return id.split(".")[0]

    def _make_random_ends(self, seq):
        """Generates reads from a random fragment within sequence.

        :param str seq: DNA or RNA sequence; method should be agnostic to type
        :return tuple: (str, str | None) R1 and R2 sequences with random ends
        """

        if self.make_amplicons:
            # Simple case where we generate reads from the termini of the target (e.g. tiles)
            frag = seq
            len_rand_frag = len(seq)
        else:
            # Typical case where we want uniform coverage over the feature(s)
            seq_len = len(seq)
            seq_indices = range(seq_len)

            # Define some initial bounds for a fragment length, then choose one
            random_len = 0
            while random_len < MIN_FRAG_LENGTH:
                # Determine a fragment length based on a normal distribution; make SD about a fourth of the mean
                random_len = int(numpy.random.normal(float(self.frag_length), self.frag_length * SD_FROM_MEAN_FACTOR))

            # Choose a random break point in the feature
            start_index = random.choice(seq_indices)
            rand_index_plus = start_index + random_len
            rand_index_minus = start_index - random_len

            # There is likely a better solution to this
            if rand_index_plus <= seq_len:
                frag = seq[start_index:rand_index_plus]
            elif rand_index_minus >= 0:
                frag = seq[rand_index_minus:start_index]
            else:
                # Last ditch effort to create a fragment in the region. Notice this does not use random_len;
                # required for small features that are smaller than the fragment length desired
                frag = seq[random.randint(0, self.slop_length):random.randint(seq_len - self.slop_length, seq_len)]

            len_rand_frag = len(frag)

        # Now that we have a fragment, make reads from it
        # Make R1 from the start to the end of the fragment
        r1 = frag[:self.read_length]
        r2 = None
        if self.paired:
            # Make R2 from the end to the start
            r2 = reverse_complement(frag[-len_rand_frag:])
            # Trim to read length if the fragment is larger
            r2 = r2[:self.read_length]
            # Alternate which strand is R1, R2
            if random.random() >= 0.5:
                r1, r2 = r2, r1
        elif random.random() >= 0.5:
            # If we have single end reads, reverse every 1 of 2 so they're not all on the same strand
            r1 = reverse_complement(r1)

        return r1, r2

    def _make_reads_from_fragment(self, seq, nreads, r1_fh, r2_fh, prefix=""):
        """Given a sequence and number of reads, generate random fragments within the sequence and write to FASTA.

        :param str seq: contig sequence
        :param int nreads: number of reads to generate
        :param file r1_fh: R1 FASTA filehandle
        :param file | None r2_fh: Optional R2 FASTA filehandle
        :param str prefix: optional prefix to include for read IDs, which will be random and unique
        """

        feature_reads = [self._make_random_ends(seq) for _ in range(nreads)]

        # Generate unique read IDs for each sequence
        feature_read_ids = make_unique_ids(nreads, prefix=prefix)

        # Output the reads to FASTA
        for (r1, r2), read_id in zip(feature_reads, feature_read_ids):

            seq_name = FASTA_HEADER_CHAR + read_id + FILE_NEWLINE
            r1_fh.write(seq_name)
            r1_fh.write(r1.upper() + FILE_NEWLINE)

            if self.paired:
                r2_fh.write(seq_name)
                r2_fh.write(r2.upper() + FILE_NEWLINE)

    def fastqs_from_features(self, output_dir=".", output_prefix=None, nreads=DEFAULT_NREADS,
                             add_snps=DEFAULT_SNPS, add_indels=DEFAULT_INDELS):
        r"""Makes realistic in silico sequencing reads from a feature file.

        :param str output_dir: Optional output directory to store generated FASTQs and BAM
        :param str output_prefix: Optional output prefix for the FASTQ(s); if None, will use same prefix as feature file
        :param int nreads: mean number of reads to generate; will draw counts from a normal distribution to simulate \
        sampling effects
        :param bool add_snps: add SNP errors to reads
        :param bool add_indels: add InDel errors to reads
        :return tuple: (str, str) path of the output FASTQ(s); if single end, the second element will be None
        """

        na_type = "rna" if self.rna else "dna"

        if output_prefix is None:
            out_prefix = replace_extension(self.feature_file, na_type)
        else:
            out_prefix = ".".join([output_prefix, na_type])

        out_prefix = os.path.join(output_dir, out_prefix)

        out_r2_fh = open(out_prefix + ".r2.fq", "w") if self.paired else None
        with open(out_prefix + ".r1.fq", "w") as out_r1_fh:

            # Generate the reads
            read_gen_args = [self.feature_file, self.ref, self.paired, self.read_length, self.frag_length,
                             self.slop_length, self.make_amplicons]

            read_gen = RnaReadGenerator(*read_gen_args) if self.rna else DnaReadGenerator(*read_gen_args)
            r1_fasta, r2_fasta = read_gen.make_reads(nreads)

            # Convert to FASTQ and add sequencing error
            fasta_to_fastq(
                in_fasta=r1_fasta, out_fastq=out_r1_fh.name, add_error_snps=add_snps, add_error_indels=add_indels)

            safe_remove((r1_fasta,))

            if self.paired:
                fasta_to_fastq(
                    in_fasta=r2_fasta, out_fastq=out_r2_fh.name, add_error_snps=add_snps, add_error_indels=add_indels)
                out_r2_fh.close()
                safe_remove((r2_fasta,))
                return out_r1_fh.name, out_r2_fh.name

            return out_r1_fh.name, None


class DnaReadGenerator(ReadGenerator):
    """Class for generation of in silico DNA sequencing reads."""

    def __init__(self, feature_file, ref, paired=True, read_length=DEFAULT_READ_LEN,
                 frag_length=DEFAULT_FRAG_LEN, slop_length=75, make_amplicons=False):
        """Ctor for DnaReadGenerator."""

        super(DnaReadGenerator, self).__init__(feature_file, ref, paired, rna=False, read_length=read_length,
                                               frag_length=frag_length, slop_length=slop_length,
                                               make_amplicons=make_amplicons)

        # Assuming the features are exons, slop the region to include intronic sequences, the defining feature of DNA
        # emulates hybrid capture based enrichment
        # Use larger slop to deal with logic of generate_random_ends(), which picks a random frag point in the feature
        with tempfile.NamedTemporaryFile(suffix=".genome_file.txt", mode="w", delete=False) as genome_file:
            genome_fn = genome_file.name
            get_genome_file(ref=self.ref, output_file=genome_file.name)

        slopped_features = slop_features(
            feature_file=self.feature_file, genome_file=genome_fn, bp_left=self.slop_length, bp_right=self.slop_length)

        self.slopped_bedtool = pybedtools.BedTool(slopped_features)

        self.feature_weights = self._weight_features()

        safe_remove((genome_fn,))

    def _weight_features(self, feature_weight_min_len=300):
        """Weights individual components to enable uniform coverage.

        :param int feature_weight_min_len: min feature len for which to upweight read count
        """

        weighting_dict = collections.defaultdict(float)

        for feature in self.ff_bedtool:
            feature_name = COORD_FORMAT.format(*list(map(str, [feature.chrom, feature.start + 1, feature.stop])))
            feature_len = len(feature)

            feature_weight = 1.0
            # Only upweight longer features
            if feature_len > feature_weight_min_len:
                feature_weight = feature_len / MEAN_EXON_LEN

            weighting_dict[feature_name] = feature_weight

        return weighting_dict

    def make_reads(self, nreads=DEFAULT_NREADS):
        r"""Make some realistic DNA reads from a feature file of targeted regions.

        :param int nreads: mean number of reads/pairs to generate; random noise is added and longer features are \
        upweighted- see DnaReadGenerator.weight_features()
        :return tuple: (str, str) path of the output FASTA(s); if single end, the second element will be None

        In the future, consider creating sample-specific errors (chemical changes like deamination, oxidation)
        """

        r2_is_temp = False if self.paired else True
        with tempfile.NamedTemporaryFile("w", suffix=".dna.r1.fa", delete=False) as out_r1_fh, \
                tempfile.NamedTemporaryFile("w", suffix=".dna.r2.fa", delete=r2_is_temp) as out_r2_fh:

            # Iterate over every feature to create reads;  iter() method of BedTool returns Interval objects
            for feature in self.slopped_bedtool:

                # coords = list(map(str, [feature.chrom, feature.start + 1, feature.stop]))
                # coord_str = COORD_FORMAT.format(*coords)
                seq = extract_seq(feature.chrom, feature.start + 1, feature.stop, self.ref)

                # Now generate reads for each feature
                feature_score = str(feature.score)
                if is_number(feature_score):
                    num_reads = int(feature_score)
                else:
                    coords = list(map(str, [feature.chrom, feature.start + 1 - self.slop_length,
                                            feature.stop + self.slop_length]))
                    coord_str = COORD_FORMAT.format(*coords)
                    weighted_reads = self.feature_weights[coord_str] * nreads
                    num_reads = int(numpy.random.normal(weighted_reads, weighted_reads * SD_FROM_MEAN_FACTOR))

                r2_fh = out_r2_fh if self.paired else None
                self._make_reads_from_fragment(seq=seq, nreads=num_reads, r1_fh=out_r1_fh, r2_fh=r2_fh, prefix="")

            if self.paired:
                return out_r1_fh.name, out_r2_fh.name

            return out_r1_fh.name, None


class RnaReadGenerator(ReadGenerator):
    r"""Class for generation of in silico RNA sequencing reads.

    Warning: features should exist in 5'-3' direction under the metafeature line, regardless of strand. This enables \
    not only generation of WT transcript RNA reads, but the creation of fusions, even when on different strands.
    All that is required is a proper ordering.

    Note BED format is not a viable format for transcript annotations if normal RNA reads (WT transcripts) are desired. \
    For that, use GFF/GTF format with exon features sorted by coordinate under a descriptive metafeature with GFF \
    score equal to the expression level. Without the score, the nreads value sets equal counts between metafeatures. \
    Because of the variability in GFF formats, group_features.features_to_omit allows the gathering of information \
    on metafeatures and may need to be modified according to your need. Standard naming is NCBI's gff3 convention.
    """

    def __init__(self, feature_file, ref, paired=True, read_length=DEFAULT_READ_LEN,
                 frag_length=DEFAULT_FRAG_LEN, slop_length=DEFAULT_BP_SLOP):
        """Ctor for RnaReadGenerator."""

        super(RnaReadGenerator, self).__init__(feature_file, ref, paired, rna=True, read_length=read_length,
                                               frag_length=frag_length, slop_length=slop_length)

    # Check score
    def _group_features(self, nreads):
        """Groups features by metafeature (if present); otherwise groups into an "Unknown" metafeature.

        :param int nreads: mean number of reads/pairs to generate; draws counts from a normal distr to simulate sampling effects
        :return tuple: (dict, dict): dicts of grouped features and feature scores
        """

        # We only want exon features, which should link to the gene or transcript
        # Keep this private, because we want tight control over what features to grab metadata from
        features_to_omit = re.compile("|".join(["gene", "transcript", "mRNA" , "cDNA_match"]))

        # Attempt to group exons by transcript; iter() method of BedTool returns Interval objects
        grouped_features = collections.defaultdict(list)
        metafeature_scores = {}

        for feature in self.ff_bedtool:

            # Assuming we can identify the transcript associated with each exon, group them
            trx_id = None
            gene_id = None
            if feature.file_type in {GFF_FILETYPE, GTF_FILETYPE}:
                trx_id = str(feature.attrs[GFF_ATTR_TRANSCRIPT_ID]) if GFF_ATTR_TRANSCRIPT_ID in feature.attrs else None
                gene_id = str(feature.attrs[GFF_ATTR_GENE_ID]) if GFF_ATTR_GENE_ID in feature.attrs else None

            # Filter the feature file of metafeatures (transcripts, genes, etc.) so we have just exons, but first save
            # the metafeature score
            if features_to_omit.search(feature.fields[GFF_FEATURE_TYPE_FIELD]):

                metafeature_score = nreads
                if is_number(feature.score):
                    metafeature_score = int(float(feature.score))

                # First keep track of metafeature expression scores; ensure we have a gene or transcript to group on
                if trx_id is not None and trx_id not in PYBEDTOOLS_NULL_CHARS or \
                        gene_id is not None and gene_id not in PYBEDTOOLS_NULL_CHARS:
                    metafeature_scores[(gene_id, trx_id)] = metafeature_score

                continue

            if trx_id is not None and trx_id not in PYBEDTOOLS_NULL_CHARS or \
                        gene_id is not None and gene_id not in PYBEDTOOLS_NULL_CHARS:
                grouped_features[(gene_id, trx_id)].append(feature)
            else:
                # Note if an exon is not found as a part of a larger metafeature then it will create adjoined sequences
                # with its adjacent feature. This may or may not be desired behavior. Safest case is to provide a
                # legit GTF with metafeatures, or provide a feature file with one transcript or gene.
                grouped_features[("Unknown", "Unknown")].append(feature)

        return grouped_features, metafeature_scores

    def make_reads(self, nreads=DEFAULT_NREADS):
        r"""Makes realistic RNA reads from a feature file of targeted regions.

        :param int nreads: mean number of reads/pairs to generate; draws counts from a normal distr to simulate sampling effects
        :return tuple: (str, str) path of the output FASTA(s); if single end, the second element will be None

        Warning! This method only produces reads for GFFs with correct ordering of features in a larger \
        metafeature (e.g. transcript).

        Noncanonical uses:
        Two-way fusions can be made (even on genes on different strands), BUT the metafeature order should have \
        the 5' partner first, then the 3' partner (regardless of strand).

        Expression imbalance within a transcript (e.g. 3' exons with elevated expression) is not supported, though \
        expression at the level of the whole transcript can be modulated with the GFF/GTF score field.
        """

        r2_is_temp = False if self.paired else True
        with tempfile.NamedTemporaryFile("w", suffix=".rna.r1.fa", delete=False) as out_r1_fh, \
                tempfile.NamedTemporaryFile("w", suffix=".rna.r2.fa", delete=r2_is_temp) as out_r2_fh:

            # First group the features/exons, assuming the features are exons
            grouped_exons, metafeature_scores = self._group_features(nreads)

            # Now that we have the exons grouped, concatenate their sequences in proper order
            for metafeature, exons in grouped_exons.items():

                metafeature_seq = self._concatenate_feature_seqs(exons)

                # We assembled the transcript or fusion, time to make reads
                # Note we support differential transcript expression (differential isoform expression can be done, but
                # requires some forthought on the part of the script caller- include all isoforms in GTF)
                if metafeature in metafeature_scores:
                    num_reads = metafeature_scores[metafeature]
                else:
                    num_reads = int(numpy.random.normal(float(nreads), nreads * SD_FROM_MEAN_FACTOR))

                r2_fh = out_r2_fh if self.paired else None
                gene_id, trx_id = metafeature
                name_prefix = "{}.{}".format(gene_id, trx_id)
                self._make_reads_from_fragment(metafeature_seq, num_reads, out_r1_fh, r2_fh, prefix=name_prefix)

            if self.paired:
                return out_r1_fh.name, out_r2_fh.name

            return out_r1_fh.name, None

    def _concatenate_feature_seqs(self, exons):
        """Concatenates a list of features. Assumes the exons are ordered.

        :param list exons: list of pybedtools.BedTool objects
        :return str: sequence of the metafeature, 5'-3'
        """

        # Add up the sequence as we iterate over the exons
        metafeature_seq = ""

        for exon in exons:

            # Normalize different contig naming conventions
            # e.g. in the hg19 RefSeq GFF, chromosomes can be annotated as NC_000001.10
            # disregard contig/RefSeq versions and assume the annotation file is up to dat
            no_version_contig = self.remove_version_extension(exon.chrom)
            contig = HumanContig(no_version_contig).value

            coords = list(map(str, [contig, exon.start + 1, exon.stop]))
            seq = extract_seq(*coords, ref=self.ref)

            # Normal transcript or fusion on the same strand
            if Strand(str(exon.strand)) == Strand.MINUS:
                seq = reverse_complement(seq)

            metafeature_seq += seq

        return metafeature_seq
