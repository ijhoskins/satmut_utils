#!/usr/bin/env python3
"""Objects for generating variants."""

import collections
import logging
import pysam
import random
import re
import tempfile

import analysis.coordinate_mapper as cm
import analysis.seq_utils as su
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

tempfile.tempdir = os.getenv("SCRATCH", "/tmp")
logger = logging.getLogger(__name__)

VAR_ID_TUPLE = collections.namedtuple("VAR_ID_TUPLE", "pos, ref, alt, aa_change, mut_sig_match")


class VariantGenerator(object):
    """Class for generating all variant permutations of each codon in a transcript's CDS."""

    DEFAULT_OUTDIR = "."
    DEFAULT_OUTFILE = None
    DEFAULT_TARGETS = None
    DEFAULT_RACE_LIKE = False
    DEFAULT_EXT = "codon.permuts.vcf"
    DEFAULT_VAR_TYPE = "total"
    DEFAULT_MNP_BASES = 3
    DEFAULT_HAPLO = False
    DEFAULT_HAPLO_LEN = 12
    DEFAULT_HAPLO_SEED = 9
    DEFAULT_INFO_DELIM = ","
    DEFAULT_SNP_WEIGHT = 0.5
    DEFAULT_MUTAGENESIS_PRIMER_LEN = 25

    def __init__(self, gff, ref, mut_sig=DEFAULT_MUT_SIG, haplotypes=DEFAULT_HAPLO, haplotype_len=DEFAULT_HAPLO_LEN,
                 outdir=DEFAULT_OUTDIR, random_seed=DEFAULT_HAPLO_SEED):
        r"""Constructor for VariantGenerator.

        :param str gff: transcript GFF; must have "transcript_id" metafeature and "exon", "CDS", "start_codon", \
        and "stop_codon" features
        :param str ref: reference FASTA with contigs matching those in the GFF seqname field
        :param str mut_sig: mutagenesis signature- one of {NNN, NNK, NNS}. Default NNN.
        :param bool haplotypes: should haplotypes be created with uniform number to codon variants? Default True.
        :param int haplotype_len: max length to create haplotypes. No longer than read length.
        :param str outdir: Output directory to write results
        :param int random_seed: integer seed for haplotype generation if haplotypes=True.
        """

        file_ext = fu.get_extension(gff)
        if not (re.match(ffu.GFF_FILETYPE, file_ext) or re.match(ffu.GTF_FILETYPE, file_ext)):
            raise NotImplementedError("Input file must be a GFF/GTF filetype.")

        self.gff = gff
        self.ref = ref
        self.mut_sig = mut_sig
        self.haplotypes = haplotypes
        self.haplotype_len = haplotype_len
        self.outdir = outdir
        self.random_seed = random_seed

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        # Create the necessary mapping object for getting CDS sequence
        self.aam = cm.AminoAcidMapper(gff=gff, ref=ref, outdir=outdir, mut_sig=mut_sig)
        self.contig_lookup = su.get_contig_lookup(self.ref)

    def _create_vcf_header(self, contigs):
        """Creates a VCF header for the candidate variants.

        :param tuple contigs: Contigs to include in the header
        :return pysam.VariantHeader: VCF header object
        """

        vcf_header = pysam.VariantHeader()

        class_module_path = ".".join(
            [PROJECT_LAB, PROJECT_AUTHOR, os.path.basename(PROJECT_ROOT), __name__, self.__class__.__name__])

        aa_mapper_module_path = ".".join(
            [PROJECT_LAB, PROJECT_AUTHOR, os.path.basename(PROJECT_ROOT), cm.__name__, cm.AminoAcidMapper.__name__])

        vcf_header.add_line(vu.VCF_METADATA_CHAR + "source=" + class_module_path)

        # Add the contig IDs
        for contig in contigs:
            vcf_header.add_line(vu.VCF_CONTIG_HEADER_FORMAT.format(contig))

        new_info_ids = (
            (vu.VCF_AAM_AA_CHANGE_ID, ".", "String", "%s comma-delimited amino acid changes." % aa_mapper_module_path),
            (vu.VCF_VARTYPE_ID, ".", "String", "Variant type."),
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
            (vu.VCF_MUT_SIG_MATCH, ".", "String", "Whether or not the variant matches the mutagenesis signature."),
        )

        vu.update_header(header=vcf_header, info_ids=new_info_ids, contigs=contigs)

        return vcf_header

    def _get_codon_variant_pos(self, codon_start, ref, alt):
        """Infers the variant given a codon REF and ALT.

        :param int codon_start: coordinate 1-based start
        :param str ref: REF codon
        :param str alt: ALT codon
        :return tuple: (transcriptomic_start, REF, ALT) with 0-based start
        """

        if len(ref) != len(alt):
            raise NotImplementedError("Incorrect REF and ALT args: (%s, %s)" % (ref, alt))

        pos = 0
        first_pos_seen = False
        mm_bps = []

        for i, (b1, b2) in enumerate(zip(ref, alt)):

            if b1 != b2:
                mm_bps.append(i)

                if not first_pos_seen:
                    pos = codon_start + i
                    first_pos_seen = True

        min_index = min(mm_bps)
        max_index = max(mm_bps) + 1
        new_ref = ref[min_index:max_index]
        new_alt = alt[min_index:max_index]

        return pos, new_ref, new_alt

    def _combine_variants(self, trx_id, first_var, second_var):
        """Combines two variants into a longer-range haplotype.

        :param str trx_id:
        :param collections.namedtuple first_var: pos, ref, alt of first variant
        :param collections.namedtuple second_var: pos, ref, alt of second variant
        :return collections.namedtuple: pos, ref, alt of new variant
        """

        min_var = first_var if first_var.pos <= second_var.pos else second_var
        max_var = second_var if first_var == min_var else first_var

        end_ref_read_pos = max_var.pos + len(max_var.alt)
        # Extraction uses 1-based pos
        base_seq = su.extract_seq(contig=trx_id, start=min_var.pos + 1, stop=end_ref_read_pos, ref=self.ref)

        # We have the span's sequence; now just update with the mismatched bases from each variant
        # Paste together the components
        var_seq = list(min_var.alt) + list(base_seq[len(min_var.alt):-len(max_var.alt)]) + list(max_var.alt)
        var_seq = "".join(var_seq)

        aa_change = self.DEFAULT_INFO_DELIM.join([min_var.aa_change, max_var.aa_change])

        var_id_tup = VAR_ID_TUPLE(pos=min_var.pos, ref=base_seq, alt=var_seq, aa_change=aa_change, mut_sig_match=None)

        return var_id_tup

    def _get_all_trx_variants(self, trx_id, outfile, var_type=DEFAULT_VAR_TYPE, mnp_bases=DEFAULT_MNP_BASES):
        """Generates a VCF of all permutation variants in  the coding region of a transcript.

        :param str trx_id: transcript ID to generate variants for; only one version may be available in the input GFF
        :param str outfile: output VCF filename
        :param str var_type: one of {"snp", "mnp", "total"}
        :param int mnp_bases: report for di- or tri-nt MNP? Must be either 2 or 3. Default 3.
        :return str: name of the output VCF
        """

        vcf_header = self._create_vcf_header(contigs=(self.contig_lookup[trx_id],))

        trx_cds_info = self.aam.cds_info[trx_id]
        trx_start_pos = trx_cds_info[self.aam.CDS_INFO_CDS_START_INDEX]
        trx_cds_seq = trx_cds_info[self.aam.CDS_INFO_CDS_SEQ_INDEX]

        cds_codons = [trx_cds_seq[i:i + 3] for i in range(0, len(trx_cds_seq), 3)]
        cp = CodonPermuts(var_type=var_type, mnp_bases=mnp_bases)

        # Get the full name of the contig
        full_trx_id = self.contig_lookup[trx_id]

        # Keep SNPs and MNPs apart for differential weighting
        snps = set()
        mnps = set()

        patch_vcf = tempfile.NamedTemporaryFile(suffix=".patch.vcf", delete=False).name
        with pysam.VariantFile(patch_vcf, "w", header=vcf_header) as out_vf:

            # Don't worry about parallelizing because the max number of elements is low
            for i, codon in enumerate(cds_codons):

                codon_start = trx_start_pos + (i * 3)
                ref_aa = cm.translate(codon)

                # Get the codon alts and some add'l metadata for annotating the VCF records
                var_permut_num, var_permut_alts = cp.codon_var_alts[codon]
                aa_permut_num, _ = cp.codon_aa_alts[codon]

                for permut in var_permut_alts:

                    pos, ref, alt = self._get_codon_variant_pos(codon_start=codon_start, ref=codon, alt=permut)
                    var_type = vu.get_variant_type(ref=ref, alt=alt, split_mnps=True)

                    alt_aa = cm.translate(permut)
                    aa_change = cm.HGVS_AA_FORMAT.format(ref_aa, i + 1, alt_aa)

                    mut_info_tuple = self.aam.get_codon_and_aa_changes(trx_id=trx_id, pos=pos + 1, ref=ref, alt=alt)

                    # Only generate variants that match the mutagenesis signature
                    if mut_info_tuple.matches_mut_sig == "False":
                        continue

                    info_dict = {vu.VCF_VARTYPE_ID: str(var_type.value),
                                 vu.VCF_AAM_LOCATION_ID: mut_info_tuple.location,
                                 vu.VCF_AAM_CODON_REF_ID: mut_info_tuple.wt_codons,
                                 vu.VCF_AAM_CODON_ALT_ID: mut_info_tuple.mut_codons,
                                 vu.VCF_AAM_AA_REF_ID: mut_info_tuple.wt_aas,
                                 vu.VCF_AAM_AA_ALT_ID: mut_info_tuple.mut_aas,
                                 vu.VCF_AAM_AA_CHANGE_ID: mut_info_tuple.aa_changes,
                                 vu.VCF_AAM_AA_POS_ID: mut_info_tuple.aa_positions,
                                 vu.VCF_MUT_SIG_MATCH: mut_info_tuple.matches_mut_sig}

                    new_variant_record = out_vf.new_record(
                        contig=full_trx_id, start=pos, alleles=(ref, alt), info=info_dict)

                    out_vf.write(new_variant_record)

                    var_tup = VAR_ID_TUPLE(pos, ref, alt, aa_change, None)
                    if var_type == vu.VariantType.SNP:
                        snps.add(var_tup)
                    else:
                        mnps.add(var_tup)

            if self.haplotypes:
                # Once we have determined all codon combinations for the transcript, pick some and combine until we meet
                # the same number of variants, to create a uniform mixture with hierarchical structure.
                # For now just create haplotypes between two component codon mutants. Note this does not imply the
                # number of component SNPs is 2. These haplotypes could have up to 6 bases changed.
                all_vars = snps | mnps
                number_variants = len(all_vars)
                haplotype_vars = 0
                haplotype_snp_vars = 0
                needed_snps = int(round(self.DEFAULT_SNP_WEIGHT * number_variants))
                random.seed(self.random_seed)
                variant_pop = snps

                # The combinatorics ensure we will complete this loop provided the number of variants is large
                # (the transcript does not consist of a single or a few codons). This is generally a safe assumption.
                while haplotype_vars < number_variants:

                    # Require a greater composition of SNPs to haplotype components, otherwise MNPs will dominate
                    # and not realistically model mutagenesis efficiencies
                    # Here we have a heuristic that first collects the haplotypes comprised of only SNPs,
                    # then accumulation of haplotypes of mixed composition (either SNP-MNP or MNP-MNP)
                    if haplotype_snp_vars == needed_snps:
                        variant_pop = all_vars

                    haplotype_components = random.sample(variant_pop, k=2)
                    first_var = haplotype_components[0]
                    second_var = haplotype_components[1]

                    # Enforce that the variant positions be within the length parameter (read length)
                    haplo_span = abs(first_var.pos - (second_var.pos + len(second_var.ref)))
                    if haplo_span > self.haplotype_len or haplo_span <= self.DEFAULT_MUTAGENESIS_PRIMER_LEN:
                        continue

                    first_var_type = vu.get_variant_type(ref=first_var.ref, alt=first_var.alt, split_mnps=True)
                    second_var_type = vu.get_variant_type(ref=second_var.ref, alt=second_var.alt, split_mnps=True)

                    # We did this above but this is required for the heuristic to work properly
                    if first_var_type == vu.VariantType.SNP and second_var_type == vu.VariantType.SNP:
                        haplotype_snp_vars += 1

                    haplotype_vars += 1

                    new_var = self._combine_variants(trx_id=full_trx_id, first_var=first_var, second_var=second_var)

                    info_dict = {vu.VCF_AAM_AA_CHANGE_ID: new_var.aa_change,
                                 vu.VCF_VARTYPE_ID: "{}:{}".format(
                                     str(first_var_type.value), str(second_var_type.value))}

                    new_variant_record = out_vf.new_record(
                        contig=full_trx_id, start=new_var.pos, alleles=(new_var.ref, new_var.alt), info=info_dict)

                    out_vf.write(new_variant_record)

        vu.remove_end_info_tag(in_vcf=patch_vcf, out_vcf=outfile)
        fu.safe_remove((patch_vcf,))
        return outfile

    def workflow(self, trx_id, targets=DEFAULT_TARGETS, outfile=DEFAULT_OUTFILE,
                 var_type=DEFAULT_VAR_TYPE, mnp_bases=DEFAULT_MNP_BASES):
        """Runs the VariantGenerator workflow.

        :param str trx_id: transcript ID to generate variants for; only one version may be available in the input GFF
        :param str | None targets: optional target feature file. Only variants intersecting the target will be generated.
        :param str | None outfile: optional output VCF filename
        :param str var_type: one of {"snp", "mnp", "total"}. Default total.
        :param int mnp_bases: report for di- or tri-nt MNP? Must be either 2 or 3. Default 3.
        :return str: name of the output VCF
        """

        if outfile is None:
            output_filepath = os.path.join(self.outdir, fu.add_extension(trx_id, self.DEFAULT_EXT))
        else:
            output_filepath = os.path.join(self.outdir, os.path.basename(outfile))

        # If there are no targets, output all codon permutation variants
        if targets is None:
            out_vcf = self._get_all_trx_variants(trx_id, output_filepath, var_type, mnp_bases)
            return out_vcf

        # Otherwise intersect the codon permutation variants with the targets
        with tempfile.NamedTemporaryFile(suffix=".codon.permuts.vcf", delete=False) as temp_vcf_fh:
            temp_vcf = self._get_all_trx_variants(trx_id, temp_vcf_fh.name, var_type, mnp_bases)
            out_vcf = ffu.intersect_features(ff1=temp_vcf, ff2=targets, outfile=output_filepath, f=1.0, header=True)
            fu.safe_remove((temp_vcf,))
            return out_vcf


class CodonPermuts(object):
    """Class for enumerating all permutations of all possible codons to another codon."""

    AA_CHANGE_FORMAT = "p.{}>{}"
    ALL_CODONS = tuple(["".join((e1, e2, e3)) for e1 in su.DNA_BASES for e2 in su.DNA_BASES for e3 in su.DNA_BASES])
    TABLE_HEADER = ["Codon", "AA_alternates", "AA_changes"]

    def __init__(self, codons=ALL_CODONS, var_type="total", mnp_bases=3, mut_sig=DEFAULT_MUT_SIG):
        """Constructor for CodonPermuts.

        :param tuple codons: tuple of codons to find permutations for
        :param str var_type: one of {"snp", "mnp", "total"}
        :param int mnp_bases: report di- or tri-nt MNPs? Must be either 2 or 3. Default 3.
        :param str mut_sig: mutagenesis signature- one of {NNN, NNK, NNS}. Default NNN.
        :raises NotImplementedError: if mut_sig or mnp_bps is invalid.
        """

        self.codons = codons
        self.var_type = var_type
        self.mnp_bases = mnp_bases
        self.mut_sig = mut_sig

        if self.var_type not in {"snp", "mnp", "total"}:
            raise NotImplementedError("Not a valid var_type: %s" % self.var_type)

        if self.mnp_bases not in {2, 3}:
            raise NotImplementedError("mnp_bases %i is not one of {2, 3}." % mnp_bases)

        if mut_sig not in VALID_MUT_SIGS:
            raise NotImplementedError("mut_sig %s is not one of {NNN, NNK, NNS}." % mut_sig)

        self.codon_var_alts = self._get_all_codon_var_changes()
        self.codon_aa_alts = self._get_all_codon_aa_changes()

    @property
    def all_codons(self):
        """Returns all the WT codon sequences.

        :return tuple: tuple of codon sequences
        """
        return self.ALL_CODONS

    @property
    def mut_sig_proportions(self):
        """Gets the proportion of REF and ALT codons that match the mutagenesis signature.

        :return tuple: (proportion of REF codons with mut_sig, proportion of ALT codons with mut_sig)
        """

        ref_prop = len([codon for codon in self.codons if codon[2] not in
                        cm.MUT_SIG_UNEXPECTED_WOBBLE_BPS[self.mut_sig]]) / len(self.codons)

        matching_codons = 0
        total_codons = 0
        for k, v in self.codon_var_alts.items():
            alt_codons = v[1]
            for alt_codon in alt_codons:
                total_codons += 1
                if alt_codon[2] not in cm.MUT_SIG_UNEXPECTED_WOBBLE_BPS[self.mut_sig]:
                    matching_codons += 1

        alt_prop = matching_codons / total_codons

        return ref_prop, alt_prop

    def _get_codon_snp_alts(self, codon):
        """Gets the set of codons that may derive from the input codon by a single SNP.

        :param str codon: input codon
        :return set: set of codons that may derive from the codon by SNPs
        """

        bps = tuple([e for e in codon])
        codon_permutations = set()

        for i, b in enumerate(bps):

            b_alts = {e for e in su.DNA_BASES if e != b}

            if i == 0:
                codon_permutations |= {"".join((ba, bps[1], bps[2])) for ba in b_alts
                                       if bps[2] not in cm.MUT_SIG_UNEXPECTED_WOBBLE_BPS[self.mut_sig]}
            elif i == 1:
                codon_permutations |= {"".join((bps[0], ba, bps[2])) for ba in b_alts
                                       if bps[2] not in cm.MUT_SIG_UNEXPECTED_WOBBLE_BPS[self.mut_sig]}
            elif i == 2:
                codon_permutations |= {"".join((bps[0], bps[1], ba)) for ba in b_alts
                                       if ba not in cm.MUT_SIG_UNEXPECTED_WOBBLE_BPS[self.mut_sig]}

        return codon_permutations

    @staticmethod
    def _get_permut_alts(codon, bi):
        """Gets the set of alternate bases for each base index pair in a codon.

        :param str codon: codon to interrogate
        :param tuple bi: indices for the codon (bases to get alts for)
        :return tuple: alternates for the 1st and 2nd index bases
        """

        bi_one_b_alts = {e for e in su.DNA_BASES if e != codon[bi[0]]}
        bi_two_b_alts = {e for e in su.DNA_BASES if e != codon[bi[1]]}

        return bi_one_b_alts, bi_two_b_alts

    def _get_codon_mnp_alts(self, codon):
        """Gets the set of codons that may derive from the input codon by a single MNP.

        :param str codon: input codon
        :return set: set of codons that may derive from the codon by MNPs
        """

        if self.mnp_bases == 3:
            codons_mismatch_mut_sig = {
                c for c in self.all_codons if c[2] in cm.MUT_SIG_UNEXPECTED_WOBBLE_BPS[self.mut_sig]}

            possible_codons = set(self.all_codons) - {codon} - codons_mismatch_mut_sig
            return possible_codons

        bps = tuple([e for e in codon])
        codon_permutations = set()

        base_indices = ((0, 1), (0, 2), (1, 2))
        for i, bi in enumerate(base_indices):

            bi1_alts, bi2_alts = self._get_permut_alts(codon, bi)

            if i == 0:
                codon_permutations |= {"".join((ba1, ba2, bps[2])) for ba1 in bi1_alts for ba2 in bi2_alts}
            elif i == 1:
                codon_permutations |= {"".join((ba1, bps[1], ba2)) for ba1 in bi1_alts for ba2 in bi2_alts
                                       if ba2 not in cm.MUT_SIG_UNEXPECTED_WOBBLE_BPS[self.mut_sig]}
            elif i == 2:
                codon_permutations |= {"".join((bps[0], ba1, ba2)) for ba1 in bi1_alts for ba2 in bi2_alts
                                       if ba2 not in cm.MUT_SIG_UNEXPECTED_WOBBLE_BPS[self.mut_sig]}

        return codon_permutations

    def _get_codon_total_alts(self, codon):
        """Gets the set of codons that may derive from the input codon by a single MNP.

        :param str codon: input codon
        :return set: set of codons that may derive from the codon by SNPs or MNPs
        """

        codon_permutations = set()
        codon_permutations |= self._get_codon_snp_alts(codon)
        codon_permutations |= self._get_codon_mnp_alts(codon)
        return codon_permutations

    def _get_codon_var_changes(self, codon):
        """Gets the codon alternate variants possible for the codon.

        :param str codon: input codon
        :param set: total set of codon alternates possible
        """

        codon_alts = set()

        if self.var_type == "total":
            codon_alts = self._get_codon_total_alts(codon)
        elif self.var_type == "snp":
            codon_alts = self._get_codon_snp_alts(codon)
        elif self.var_type == "mnp":
            codon_alts = self._get_codon_mnp_alts(codon)

        return codon_alts

    def _get_all_codon_var_changes(self):
        """Returns a dict of all the possible codon nt permutations for each codon.

        :return collections.defaultdict: {str: tuple} of {codon: (num_alts, set of ALT variants arising from codon)}}
        """

        codon_alts = collections.defaultdict(tuple)

        for codon in self.codons:
            alt_vars = self._get_codon_var_changes(codon=codon)
            codon_alts[codon] = (len(alt_vars), alt_vars)

        return codon_alts

    def get_codon_aa_changes(self, codon, pos=None):
        """Gets the AA changes possible for the codon.

        :param str codon: input codon
        :param int | None pos: optional position for the codon, for consistent output formatting
        :return set: set of AA changes arising from codon
        """

        codon_alts = self.codon_var_alts[codon][1]

        aa_changes = set()
        for codon_alt in codon_alts:

            wt_aa = cm.translate(codon)
            alt_aa = cm.translate(codon_alt)

            if pos is None:
                aa_changes.add(self.AA_CHANGE_FORMAT.format(wt_aa, alt_aa))
            else:
                aa_changes.add(cm.HGVS_AA_FORMAT.format(wt_aa, pos, alt_aa))

        return aa_changes

    def _get_all_codon_aa_changes(self):
        """Gets a dict of the total AA changes possible for all codons.

        :return collections.defaultdict: {str: tuple} of {codon: (num_alts, set of AA changes arising from codon)}}
        """

        codon_counts = collections.defaultdict(tuple)

        for codon in self.codons:
            alt_aas = self.get_codon_aa_changes(codon=codon, pos=None)
            codon_counts[codon] = (len(alt_aas), alt_aas)

        return codon_counts

    def dump_table(self, outdir="."):
        """Generates a table of the codon counts dictionary.

        :param str outdir: output directory to write codon tables
        """

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        table_name = "codon_lookup.var_type-{}.mnp_bases-{}.txt".format(self.var_type, self.mnp_bases)

        with open(os.path.join(outdir, table_name), "w") as outfile:

            outfile.write(fu.FILE_DELIM.join(self.TABLE_HEADER) + fu.FILE_NEWLINE)

            for k, v in self.codon_aa_alts.items():

                outfile.write(fu.FILE_DELIM.join([k, str(v[0]), ",".join(v[1])]) + fu.FILE_NEWLINE)

    # The following static methods are meant to be used without first instantiating the class- separate from class?
    @staticmethod
    def sum_var_changes(orf_seq, var_type="total", mnp_bases=3):
        """Gets the total number of codon variants possible given an ORF nt sequence.

        :param str orf_seq: open reading frame sequence.
        :param str var_type: "one of {"snp", "mnp", "total"}
        :param bool with_mnps: should MNPs be considered in addition to SNPs? Default False.
        :param int mnp_bases: report for di- or tri-nt MNP? Must be either 2 or 3. Default 3.
        :return int: number of variants possible by mutating orf_seq via SNPs and/or MNPs.
        """

        cp = CodonPermuts(var_type=var_type, mnp_bases=mnp_bases)

        orf_codons = [orf_seq[i:i + 3] for i in range(0, len(orf_seq), 3)]

        total_variants = 0
        for codon in orf_codons:

            if codon not in cp.codon_var_alts:
                raise NotImplementedError("Did not find codon %s" % codon)
            else:
                total_variants += cp.codon_var_alts[codon][0]

        return total_variants

    @staticmethod
    def sum_aa_changes(orf_seq, var_type="total", mnp_bases=3):
        """Gets the total number of AA changes possible given an ORF nt sequence.

        :param str orf_seq: open reading frame sequence
        :param str var_type: one of {"snp", "mnp", "total"}
        :param bool with_mnps: should MNPs be considered in addition to SNPs? Default False.
        :param int mnp_bases: report for di- or tri-nt MNP? Must be either 2 or 3. Default 3.
        :return int: number of AA changes possible by mutating aa_seq via SNPs and/or MNPs.
        """

        cp = CodonPermuts(var_type=var_type, mnp_bases=mnp_bases)

        orf_codons = [orf_seq[i:i + 3] for i in range(0, len(orf_seq), 3)]

        total_aa_changes = 0
        for codon in orf_codons:

            if codon not in cp.codon_aa_alts:
                raise NotImplementedError("Did not find codon %s" % codon)
            else:
                total_aa_changes += cp.codon_aa_alts[codon][0]

        return total_aa_changes


class AminoAcidTypes(CodonPermuts):
    """Class for enumerating all possible AA chemical type changes from codon permutations."""

    POSITIVE_CHARGED_AAS = frozenset(("R", "H", "K"))
    NEGATIVE_CHARGED_AAS = frozenset(("D", "E"))
    POLAR_UNCHARGED_AAS = frozenset(("S", "T", "N", "Q"))
    NONPOLAR_UNCHARGED_AAS = frozenset(("A", "V", "I", "L", "M", "F", "Y", "W"))
    OTHER_AAS = frozenset(("C", "G", "P"))
    AA_TYPE_CHANGE_FORMAT = "{}:{}"

    def __init__(self, codons=CodonPermuts.ALL_CODONS, var_type="total", mnp_bases=3):
        """Ctor for AminoAcidTypes.

        :param tuple codons: tuple of codons to find permutations for
        :param str var_type: one of {"snp", "mnp", "total"}
        :param int mnp_bps: report for di- or tri-nt MNP? Must be either 2 or 3. Default 3.
        :raises NotImplementedError: if the mnp_bps was not either 2 or 3.
        """

        super(AminoAcidTypes, self).__init__(codons, var_type, mnp_bases)

        self.aa_type_counts = self.get_aa_type_counts()
        self.aa_type_probs = self.get_aa_type_probabilities()

    def get_aa_type(self, aa):
        """Gets the type of an AA.

        :param str aa: AA in single letter format
        :return str: type of the AA
        """

        if aa in self.POSITIVE_CHARGED_AAS:
            return "Positive_charged"
        elif aa in self.NEGATIVE_CHARGED_AAS:
            return "Negative_charged"
        elif aa in self.POLAR_UNCHARGED_AAS:
            return "Polar_uncharged"
        elif aa in self.NONPOLAR_UNCHARGED_AAS:
            return "Nonpolar_uncharged"
        elif aa in self.OTHER_AAS:
            return "Other"
        elif aa == "*":
            return "Stop"
        else:
            raise RuntimeError("Unknown AA %s." % aa)

    def get_aa_type_counts(self):
        """Gets the counts for each AA type change."""

        aa_type_counts = collections.defaultdict(int)

        for codon, v in self.codon_aa_alts.items():

            ref_type = self.get_aa_type(cm.translate(codon))
            alt_types = v[1]

            for at in alt_types:

                alt_aa = at[-1]
                alt_type = self.get_aa_type(alt_aa)
                type_change = self.AA_TYPE_CHANGE_FORMAT.format(ref_type, alt_type)
                aa_type_counts[type_change] += 1

        return aa_type_counts

    def get_aa_type_probabilities(self):
        """Gets the probabilities for each AA type change."""

        aa_type_probs = collections.defaultdict(float)

        total_counts = sum(self.aa_type_counts.values())
        for k, v in self.aa_type_counts.items():
            aa_type_probs[k] = v / total_counts

        return aa_type_probs

    def dump_aa_type_probs(self, outdir="."):
        """Generates a table of the AA type probs dictionary.

        :param str outdir: output directory to write to
        """

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        table_name = "aa_type_probs.var_type-{}.mnp_bases-{}.txt".format(self.var_type, self.mnp_bases)

        with open(os.path.join(outdir, table_name), "w") as outfile:

            outfile.write(fu.FILE_DELIM.join(["AA_type_change", "Probability"]) + fu.FILE_NEWLINE)

            for k, v in self.aa_type_probs.items():

                outfile.write(fu.FILE_DELIM.join([k, str(v)]) + fu.FILE_NEWLINE)
