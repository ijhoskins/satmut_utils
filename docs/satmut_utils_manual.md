![](satmut_utils_logo.png)

# satmut\_utils manual

## Table of Contents

[Installation](#Installation)

[Reference files](#Reference-files)

[satmut_utils 'sim'](#satmut_utils-sim)

[satmut_utils 'call'](#satmut_utils-call)

[Code examples](#Code-examples)

['sim' code examples](#sim-code-examples)

['call' code examples](#call-code-examples)

[satmut_utils command-line interface](#satmut_utils-command-line-interface)

[Common options](#Common-options)

['sim' options](#sim-options)

['call' options](#call-options)

[Accessory scripts](#Accessory-scripts)

[Tests](#Tests)

satmut\_utils is a Python package for simulation and variant calling of saturation mutagenesis data. The two main subcommands are:
1. 'sim'
2. 'call'

satmut\_utils commands are designed to simulate and call variants in paired-end, targeted sequencing reads. Alignments to a mature mRNA reference (contiguous, spliced coding sequence with possible untranslated regions) are expected. Genome-wide and transcriptome-wide variant calling is not supported.


'sim' and 'call' support two different types of paired-end read chemistries. The first supported chemistry is a tiled amplicon PCR approach with interleaved sets of PCR amplicons. In this mode, R1 and R2 start at primer ends.

The second supported chemistry is a RACE-like approach (for example, Anchored Multiplex PCR; Zheng et al. 2014), where R1 starts at a variable fragment end and R2 starts at a primer. RACE-like PCR may be used to enrich entire coding regions in one PCR, in addition to targeting of RNA 5' and 3' ends.

satmut\_utils 'call' allows for unique molecular indices (UMIs) at the start of R1 in either library preparation chemistry. Currently, it does not support analysis of UMIs on both reads.

Furthermore, satmut\_utils does not support barcode sequencing, wherein sequencing is used to link a variant with a barcode, followed by readout of unique barcodes.

## Installation

Currently, only Linux and MacOSX operating systems are supported. To get started, follow these steps:

1. If conda is not installed, install miniconda for managing environments. See this [link](https://docs.conda.io/en/latest/miniconda.html) for installation on your particular architecture.

2. Clone the satmut\_utils repository:
```
git clone https://github.com/CenikLab/satmut_utils.git
SATMUT_ROOT="$PWD/satmut_utils"
```

3. Execute the provided shell script to generate the satmut\_utils environment, install the package, and optionally download curated reference files, which are required if using Ensembl identifiers ([see Reference files](#Reference-files)). Finally, activate the satmut\_utils environment.
```
REF_DIR="$HOME/satmut_utils_refs"
$SATMUT_ROOT/install_satmut_utils.sh -h
$SATMUT_ROOT/install_satmut_utils.sh -t -g -r "$REF_DIR" "$SATMUT_ROOT"
conda activate satmut_utils
```

You are now ready to call the command-line executable ```satmut_utils```

satmut\_utils is the primary command, with subcommands 'sim' and 'call'.

## Reference files

For convenience, a curated transcriptome of primary human transcripts from [APPRIS](https://apprisws.bioinfo.cnio.es/landing_page/) is provided, allowing the user to pass an Ensembl gene or transcript ID. However, if the requested Ensembl ID is not found in this set, custom reference files must be passed, which include:

A. Transcript reference (FASTA)

B. Transcript annotations (GFF)

C. GFF reference (FASTA)

Common transcript annotations available in GFF format map exon coordinates in the genome. For this case, file A should specify a transcript FASTA and file C should specify the genome FASTA. This is the default satmut\_utils configuration when using Ensembl IDs.

In typical saturation mutagenesis experiments a transgene is expressed from a vector. In this case, the user ideally sets file A and C to a custom composite (vector + target) reference FASTA, then makes a custom annotation GFF (B). Custom references are useful for mapping PCR tiles that span the vector-transgene junctions. Because local alignment is employed in satmut\_utils 'call', variants near these junctions may be clipped along with vector sequence unless custom reference files are provided.

Providing a genome-based GFF reference file (C) with a custom transcript reference and annotations (A, B) allows alignment to novel transcript isoforms, such as those with 5' or 3' extensions. This is possible so long as the transcript reference FASTA and corresponding exon annotations comprise a contiguous (spliced) coding sequence (possibly with noncoding 5' and 3' untranslated regions).

If you already have the GRCh38 FASTA, you may use it provided it uses Ensembl contig nomenclature (not the NCBI "NC_" or UCSC "chr" nomenclature). That is, the chromosome names should start with a single integer for autosomes, and X, Y, and MT for the sex chromosomes and mitochondrial chromosome, respectively. This is important for compatibility with the curated transcript annotations when using Ensembl identifiers. If the FASTA meets these requirements, it should be copied or linked in the REF\_DIR (optionally set by install\_satmut\_utils.sh -r option). The genome should have an index generated with samtools faidx, which is generated by default if running the installation script.

## satmut_utils sim

'sim' makes certain design decisions that may impact the performance of 'sim' for your particular application.

1. *'sim' edits variants in only those fragments with coverage from both mates of the pair*

For compatibility with the 'call' subcommand, which requires mate concordance for a variant call, 'sim' will generate the variants in only those fragments with supporting coverage from both mates. We term this concordant alternate observations (CAO). The denominator to the variant frequency calculation (DP) is also fragment-based. That is, depth is not counted twice at positions where the mates overlap. The resulting concordant allele frequency, CAF=CAO/DP, is a more conservative measure of variant abundance.

Note that satmut\_utils 'sim' truth frequencies may differ from variants quantified by other callers due to this design choice. Thus, validation of 'sim' results is only recommended with satmut\_utils 'call', as workflows were designed to follow the same filtering logic, enabling high-accuracy variant editing and subsequent calling.

2. *'sim' relies on a heuristic to select fragments for editing*

To enable editing of multiple variants at a single position, as well as prohibit unintentional phasing of variants, 'sim' employs a heuristic to enforce several rules. The designs assumes that the collection of edited variant frequencies may not exceed 1. Thus, 'sim' assumes that **all variants are configured for editing at low frequencies (<1:1000)**. If the user seeks to edit variants at higher frequencies, the pool of fragments that are amenable for editing rapidly dwindles. This is because 'sim' prohibits phasing of edited variants (each fragment/mate pair may be edited only once).

This constraint is particularly problematic if the user seeks to edit even a few variants at germline-like frequencies (0.5 - 1). 'sim' will raise a InvalidVariantConfig exception if the sum of variant frequencies across all variants in the input VCF exceeds 1. To ignore this exception and edit as many variants as possible, pass --force\_edit. The output truth VCF will indicate all *configured* variants, whether or not they were edited.

The heuristic for selecting reads enforces the following rules:

1) Any fragment selected for editing has not previously been selected for editing of another variant. This requirement ensures variants edited at nearby positions are never phased, which may lead to simultaneous false negative and false positive calls. Note this applies to variants that are different records in the VCF. To create phased variants (e.g. MNPs), provide single VCF records with REF and ALT fields that span multiple bases.

2) Any fragment selected for editing must have error-free coverage in both mates within a symmetric buffer spanning the variant position. This prohibits editing a true variant adjacent to an error, which would convert the variant to higher order (for example, SNP to di-nucleotide MNP).

3) If a primer BED file has been provided, 'sim' enforces that no part of the variant span (POS + REF field length) intersects a read segment arising from a synthetic primer. Passing a primer BED file to 'sim' is highly recommended because it enables the avoidance of editing at read termini. Otherwise, under subsequent local alignment of edited reads (in satmut\_utils 'call'), variants at read termini may be clipped, leading to their loss from the alignment and causing false negatives.

Note that rule 3 does *not* mean variants cannot be edited within or overlapping primer regions, as long as the position has read-through coverage (typically provided by an adjacent PCR amplicon/tile). 

Collectively, while effectively constraining the number and frequency of variants that can be edited at once, these rules ensure high fidelity of variant editing at ultra-low frequencies.

## satmut\_utils call

'call' makes certain design decisions that may impact the performance of variant calling for your particular application.

1. *'call' enforces mate pair concordance for candidate calls*

To reduce the impact of single-stranded sequencing errors, satmut\_utils requires the same base call is made in both mates of the pair for enumeration of a variant call. Note this does not handle polymerase errors made during PCR enrichment.

2. *'call' applies quality filters prior to verifying mate concordance*

Before pairs are considered for mate concordance, reads or base calls may be filtered by read edit distance (--max\_nm) and base quality (--min\_bq) filters.

3. *'call' uses a unique algorithm to resolve SNPs and MNPs contained in the same fragment*

The general algorithm for variant calling is as follows. See satmut\_utils.src.analysis.variant\_caller.VariantCaller for implementation.

A. Apply quality filters on read edit distance and base quality.

B. Find all concordant mismatched base calls between read mates of a pair.

C. For D-F, ensure no mismatch is part of more than one MNP call.

D. At a mismatch, determine if there is a second mismatch downstream within the window (--max\_mnp\_window).

E. If D), determine if there is third mismatch downstream within the window. If so, call a tri-nucleotide MNP.

F. If D), and the next two mismatches are not adjacent, call a di-nucleotide MNP.

G. After all MNPs within pair have been called, call remaining mismatches as SNPs.

This algorithm leads to specific calling expectations for consecutive mismatch runs/tracts with nearby isolated mismatches. In these cases, consecutive mismatch runs are called as compact MNPs, as opposed to including an isolated mismatch as part of the MNP call.

### Mate pair concordance

Enforcing concordance may lead to under-quantification of variant frequencies in cases where a significant proportion of pairs do not have appreciable overlap between R1 and R2. Insufficient overlap between mates may be due to inadequate primer design and/or overzealous 3’ base quality trimming. Additionally, for RACE-like libraries, a high proportion of fragment lengths greater than the read length may also limit sensitivity.

### satmut\_utils frequency calculation

The denominator for the variant frequency calculation- fragment coverage depth (DP)- is determined after filtering on read edit distance (NM tag) and mate pair overlap (concordance). If a primer BED file is provided to satmut\_utils, read segments originating from synthetic primer sequence do not contribute to fragment coverage depth.
When the MNP span (variant reference coordinates) covers both a synthetic primer position and an adjacent position in the amplicon, the minimum depth of coverage in the MNP span is used as the denominator for the concordant allele frequency: CAF=CAO/DP.

### satmut\_utils 'call' read preprocessing steps

satmut\_utils 'call' supports multiple methods for moderating false positives calls. Application of these methods is optional but improves specificity.

1. Synthetic primer base quality masking

To prohibit false positive calls arising from primer synthesis errors, satmut\_utils 'call' provides a masking strategy to demarcate read segments originating from synthetic primer sequence. In masking, base qualities for these segments are set to 0 and omitted from variant calls. This step requires the user to provide a primer BED file specifying the primers used in target enrichment. 

The read subsequences to mask depend on the library preparation chemistry and mate read identity and orientation. For amplicon libraries, both R1 and R2 may be masked at both ends. For RACE-like data, only the 3’ end of R1 and the 5’ end of R2 are masked. This assumes the presence of a unique molecular index (UMI) at the start of R1, which is extracted through consensus deduplication (see #2 below).

Alignments are intersected with primers with ‘bedtools intersect -bed -wa -wb’. The resulting BED file is grouped with ‘bedtools groupby -o collapse’ to group the intersecting primers for each read, then custom satmut\_utils masks primer subsequences in intersecting reads. Any primers which start within 15 nt of the read termini are selected for masking, which captures reads that have undergone trimming or clipping in the primer subsequence. This design choice also assumes primers are not designed back-to-back (<15 nt start coordinate offset). After masking reads that intersect with primers, reads that do not intersect are concatenated to recover all input reads. See satmut\_utils.src.analysis.read\_preprocessor.ReadMasker for implementation.  

2. Unique-molecular-identifier (UMI)-based read consensus deduplication

UMIs tag unique molecules prior to target enrichment with PCR, allowing for a more accurate estimation of molecular counts. UMI-based consensus deduplication improves the specificity of calls and variant frequency accuracy, as PCR jackpots are deduplicated/collapsed and some PCR and sequencing errors may be removed in the process of generating the consensus read.

This option assumes UMIs are at the start of R1, and utilizes UMI-tools (Smith, Heger, Sudbery, Genome Res. 2017) to group read duplicates based on R1 UMI-POS, where POS is the R1 aligned start position. In grouping pairs, the template length is ignored by passing --ignore-tlen. satmut\_utils then employs a consensus deduplication workflow to correct errors within duplicate groups and generate a consensus R1 and R2. The consensus base call at each aligned position of the duplicate group is determined by majority vote. In ties (two duplicates in the UMI group), if one read matches the reference, the reference base is chosen. Otherwise, the base call with the higher base quality is selected. If both bases have the same quality, the consensus base is chosen at random.

Consensus deduplication is supported for RACE-like libraries. In one RACE-like library preparation method (Anchored Multiplex PCR, Zheng et al. Nat Med 2014), the R2 starts at a primer and R1 contains a UMI-adapter which ligates to an internal fragmentation site or the 5’ or 3’ end of the molecule.

For RACE-like libraries, as a result of passing umi\_tools the --ignore-tlen option, R2s that share the same R1 UMI-POS but do not share the same R2 start coordinate may be merged into a R2 consensus contig. This logic ensures accurate fragment depth of coverage reporting for RACE-like data, but may create R2 contigs with read lengths greater than the sequencing read length. It may also lead to a small proportion of R2 contigs with internal unknown base calls (N), which arise when merged R2s do not overlap one another.

To disable such merging of RACE-like R2s into larger contigs during consensus deduplication, the user can provide --primer\_fasta. When a primer FASTA is provided, the primer from which R2 originates is appended to the UMI prior to grouping. This ensures R2s emanating from different primers are assigned to separate groups, despite sharing the same R1 UMI-POS. Note that fragment depth of coverage will not be as accurate with this option. However, this option prohibits deletion artifacts arising from merging of non-overlapping R2s. See also the --contig\_del\_threshold option.

If --primer_fasta is provided, read names (QNAME) are modified with a R2 primer barcode (last 16 nt of primer) by using fuzzy matching of all primers to each R2. This step helps moderate artifactual consensus contigs arising from reads across separate amplicons in RACE-like chemistries (--race\_like). 

To recap, when provided RACE-like data, satmut\_utils assembles R2 contigs from multiple R2s sharing a common R1 UMI-position. While not currently supported, dual UMIs on R1 and R2 may be possible by modification of the source code that calls UMI-tools.


## Code examples

Parameter help:
```
satmut_utils -h
satmut_utils sim -h
satmut_utils call -h
```

Common arguments to both 'sim' and 'call' subcommands should be provided first, then the subcommand, and then the subcommand-specific arguments.

It is recommended that a new output directory is created for each job. Default is to output to the current directory.

```OUTPUT_DIR="/tmp/satmut_utils_test"```

### 'sim' code examples

Run the 'sim' workflow by providing a BAM containing paired-end, single-contig alignments, and a VCF file specifying variants and their desired frequencies. 

```
TEST_DIR="satmut_utils/src/tests/test_data"
OUTPUT_DIR="/tmp/satmut_utils_test"
satmut_utils -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed sim -f -a $TEST_DIR/CBS_sim.bam -v $TEST_DIR/CBS_sim.vcf
```

The BAM file should be generated with satmut\_align to ensure required alignment tags (e.g. MD tag) are present.

See satmut\_utils/src/tests/test\_data/cbs\_sim.vcf for an example input VCF with configured AF INFO tags. 

Typically the alignments are from a non-mutagenized negative control library and the variants are subsamples from a VCF containing all codon permutations matching the mutagenesis signature across the target region. See satmut\_utils/src/prototype/run\_variant_generatory.py to generate variants in a target *de novo*.

See 'sim' requirement 4 in the section above for a description of primer masking in the context of simulation. The primer BED must have a strand field. See satmut\_utils/tests/test\_data/CBS\_insilico_primers.bed for an example BED file.

Specify that the alignments are from a RACE-like library preparation chemistry, and use a new random seed to select different read pairs for editing:
```
TEST_DIR="satmut_utils/src/tests/test_data"
OUTPUT_DIR="/tmp/satmut_utils_test"
satmut_utils -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed --race_like sim -f -a $TEST_DIR/CBS_sim.bam -v $TEST_DIR/CBS_sim.vcf -s 99
```

The test data is not from RACE-like chemistry and the above code snippet is for illustrative purposes only.

### 'sim' outputs

The 'sim' workflow outputs paired FASTQs, a realigned BAM file (bowtie2 global alignment mode), and a truth VCF containing expected variants and their expected concordant counts (CAO) and frequencies (CAF). The user should be aware that edited reads that do not realign will not give rise to variant calls in downstream analysis.

### 'call' code examples

Run 'call' by specifying an Ensembl transcript or gene ID, the directory containing curated reference files, and adapters to trim (or -v, --omit\_trim):
```
# NOT RUN
satmut_utils -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz --r1_fiveprime_adapters TACACGACGCTCTTCCGATCT --r1_threeprime_adapters AGATCGGAAGAGCACACGTCT --r2_fiveprime_adapters AGACGTGTGCTCTTCCGATCT --r2_threeprime_adapters AGATCGGAAGAGCGTCGTGTA
```

More than one 5' adapter and more than one 3' adapter are often needed to additionally trim vector sequences (e.g. attB sites) from reads of terminal PCR tiles that span the vector-transgene junctions. In this case, provide multiple comma-delimited adapters:

```
--r1_fiveprime_adapters TACACGACGCTCTTCCGATCT,CAAGTTTGTACAAAAAAGTTGGC
```

For an Ensembl gene ID with more than one transcript isoform, satmut\_utils will select a primary transcript in the curated APPRIS transcriptome. If the Ensembl ID is not in the curated set of primary transcripts, the user must provide custom reference files:
```
TEST_DIR="satmut_utils/src/tests/test_data"
OUTPUT_DIR="/tmp/satmut_utils_test"
satmut_utils -r $TEST_DIR/CBS.fa -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed call -1 $TEST_DIR/CBS_sim.R1.fq.gz -2 $TEST_DIR/CBS_sim.R2.fq.gz -v -m 1 -g $TEST_DIR/CBS.gff -k $REF_DIR/GRCh38.fa.gz
```
Additional files may be passed, such as a primer and target BED file.

The primer BED file must have a strand field. See satmut\_utils/src/tests/test\_data/CBS\_insilico_primers.bed for an example BED file.

Passing the target BED file only impacts reporting of variants and does not speed up analysis by intersection of alignments with the target prior to variant calling. This is due to the strict requirement that paired reads are input to variant calling; intersection of reads against the target prior to calling often leads to dropout of reads during alignment to vector-transgene references.

### 'call' outputs

The 'call' workflow produces a VCF of candidate variant calls, a tab-delimited summary file, and a bedgraph file reporting fragment coverage across the reference.

The output VCF and its corresponding summary.txt file contain records for each mismatched base in an MNP, so that quality information for the mismatches can be used for machine learning-based error correction. See [output fields](#satmut\_utils-'call'-output-fields) or the output VCF header for column/feature descriptions. 

### satmut\_utils 'call' output fields

The output tab-delimited summary.txt file contains the standard VCF fields and INFO tag-value pairs split into unique columns. The VCF INFO fields are described below.

POS\_NT: Coordinate position of component nucleotide.

REF\_NT: Component reference nucleotide.

ALT\_NT: Component alternate nucleotide.

UP\_REF\_NT: -1 upstream reference nucleotide.

DOWN\_REF\_NT: +1 downstream reference nucleotide.

DP: Fragment-based depth of coverage after quality filters (pair overlap, edit distance).

CAO: Concordant alternate observations- alternate found in both mates.

NORM\_CAO: Mate-concordant observations per 1000000 pairs.

CAF: Concordant allele frequency in range (0,1). Calculated as CAO/DP.

R1\_PLUS\_AO: Read 1 alternate observations on (+) strand.

R1\_MINUS\_AO: Read 1 alternate observations on (-) strand.

R2\_PLUS\_AO: Read 2 alternate observations on (+) strand.

R2\_MINUS\_AO: Read 2 alternate observations on (-) strand.

R1\_PLUS\_MED\_RP: Read 1 (+) strand median read position supporting call.

R1\_MINUS\_MED\_RP: Read 1 (-) strand median read position supporting call.

R2\_PLUS\_MED\_RP: Read 2 (+) strand median read position supporting call.

R2\_MINUS\_MED\_RP: Read 2 (-) strand median read position supporting call.

R1\_PLUS\_MED\_BQ: Read 1 (+) strand median Phred base quality supporting call.

R1\_MINUS\_MED\_BQ: Read 1 (-) strand median Phred base quality supporting call.

R2\_PLUS\_MED\_BQ: Read 2 (+) strand median Phred base quality supporting call.

R2\_MINUS\_MED\_BQ: Read 2 (-) strand median Phred base quality supporting call.

R1\_PLUS\_MED\_NM: Read 1 (+) strand median edit distance supporting call.

R1\_MINUS\_MED\_NM: Read 1 (-) strand median edit distance supporting call.

R2\_PLUS\_MED\_NM: Read 2 (+) strand median edit distance supporting call.

R2\_MINUS\_MED\_NM: Read 2 (-) strand median edit distance supporting call.

LOCATION: Location of the variant in the transcript. One of {CDS, 5\_UTR, 3\_UTR, intergenic, untranslated}.

REF\_CODON: Comma-delimited reference codon(s). NA if the variant is out of CDS bounds.

ALT\_CODON: Comma-delimited alternate codon(s). NA if the variant is out of CDS bounds.

REF\_AA: Comma-delimited reference amino acid(s). NA if the variant is out of CDS bounds.

ALT\_AA: Comma-delimited alternate amino acid(s). NA if the variant is out of CDS bounds.

AA\_CHANGE: Comma-delimited amino acid change(s). NA if the variant is out of CDS bounds.

AA\_POS: Comma-delimited amino acid position(s). NA if the variant is out of CDS bounds.

MATCHES\_MUT\_SIG: Whether or not the variant matches the mutagenesis signature.

## satmut\_utils command line interface

satmut\_utils provides the 'sim' and 'call' workflow as subcommands, which have common and unique options.

### Common options

1. -i, --ensembl_id

An Ensembl gene (ENSG) or transcript (ENST) identifier containing the minor version number. For example, ENST00000398165.7.

2. -r, --reference

A custom reference FASTA for a single transcript or target transgene. Will be indexed by samtools and bowtie2 if not already. This option is mutually exclusive with --ensembl\_id.

3. -x, --reference\_dir

If using --ensembl\_id, a reference directory containing curated a transcriptome FASTA, transcript annotations (GFF), and the genome FASTA.

4. -z, --race\_like

Flag to indicate input reads or alignments are from a rapid-amplification-of-cDNA-ends (RACE)-like target enrichment strategy (e.g. Anchored Multiplex PCR). Default is tiled PCR enrichment strategy with R1 and R2 both starting at the 5' end of opposing primers.

5. -p, --primers

Primer BED file with six fields, containing the contig, start, stop, name, score, and strand fields. BED files have a 0-based start position.

6. -o, --output\_dir

Optional output directory to write results to. Default is current working directory.

7. -j, --nthreads

Number of additional alignment threads (bowtie2), and threads to use for BAM sorting (samtools). Multiprocessing is not currently supported for satmut\_utils 'sim' or 'call' Python code.

8. -e, --max\_nm

Maximum edit distance for either mate of a pair to be considered for simulation and variant calling. Default 10.

### 'sim' options

1. -a, --alignments

This is a BAM file containing paired-end alignments to a single transcript or contig. If a BAM file is not generated, use satmut\_align for bowtie2 alignment. If reads have been aligned with another aligner, there is no guarantee 'sim' will complete without error, as certain alignment tags (e.g. MD tag) are required for 'sim'. 

2. -v, --vcf

Variant Call Format file that contains an in-line INFO tag (AF) specifying the desired variant frequency. See satmut\_utils/src/tests/test\_data/CBS\_sim.vcf for an example.

3. -b, --edit_buffer

When selecting reads to edit, verify the absence of errors within this buffer about the edit span (range of reference coordinates having mismatches). The read segment must match the reference within the span (POS - edit\_buffer, POS + edit_buffer + REF length). Increasing this value prevents variant conversion and/or unexpected clipping of the variant from the termini of the read under local alignment. Decreasing this value may be needed in cases where variants are being edited near read termini and there are no adjacent PCR amplicons to provide coverage at the position.

4. -f, --force\_edit

By default, 'sim' will raise a InvalidVariantConfig Exception if the sum of variant frequencies across all variants in the input VCF exceed 1. Invalid configurations indicate that all variants may not be edited due to insufficient read depth of coverage. In this case, to edit as many variants as possible, provide this flag.

InvalidVariantConfig exceptions are guaranteed to result in unedited variants if alignments from a single PCR tile are being edited. However, if variants are to be edited into alignments spanning multiple PCR tiles, all variants *may* be edited despite an InvalidVariantConfig exception. This flag is particularly useful for maximizing the number of variants edited for multi-tile alignments.

5. -y, --random\_seed

Integer seed to use for pseudorandom qname (read name) sampling. This may be used to select different reads for editing of variants.

### 'call' options

1. -1, --fastq1

R1 FASTQ of the pair. If using simulated data that has previously been trimmed, pass the -v/--omit\_trim option.

If -p (--primers) is provided, the read names must match either Illumina format or consist of a unique integer to facilitate primer base quality masking.

2. -2, --fastq2

R2 FASTQ of the pair.

3. -v, --omit\_trim

Useful for input where adapters have been previously trimmed (e.g. simulated data). Trimming when no adapters are present may degrade the data quality by nonspecific trimming of the reads. For this reason, it is best to directly re-align these inputs.

4. --r1\_fiveprime\_adapters

Comma-delimited R1 5' adapters, or None (default) if no R1 5' adapters exist.

5. --r1\_threeprime\_adapters

Comma-delimited R1 3' adapters, or None (default) if no R1 3' adapters exist.

6. --r2\_fiveprime\_adapters

Comma-delimited R2 5' adapters, or None (default) if no R2 5' adapters exist. 

7. --r2\_threeprime\_adapters

Comma-delimited R2 3' adapters, or None (default) if no R2 3' adapters exist.

8. -g, --transcript\_gff

Transcript GFF where **features are ordered from 5' to 3', regardless of strand**. For examples, see src/tests/test\_data/gencode.v29.annotation.gtf for a standard genome-based GFF example, or src/tests/test\_data/CBS_pEZY3.gff for a custom, composite vector GFF example.

IMPORTANT: the GFF seqname field should match a contig name in the GFF reference FASTA (-k/--gff\_reference). Furthermore, each GFF record should have a "transcript\_id" attribute in the attribute field that matches the contig name in the reference FASTA (-r/--reference).

The minimum expected GFF records have a feature field (3rd column of the GFF) specifying one of the following:
```
exon, CDS, stop_codon
```

Additional allowable features in the feature field (3rd column) of the GFF include:
```
gene, transcript, start_codon
```

While these latter features/records are not necessary, they are often helpful for thorough annotation of the transcript. While a start\_codon feature is contained within the most 5' CDS feature, a stop\_codon feature is *not* contained in the most 3' CDS feature, and must be provided if the user seeks to call variants altering the stop codon (nonstop variants). 

Otherwise, one is free to specify exon and CDS features that annotate the coding and noncoding portions of each mature transcript, provided the exon features sum up to the full reference FASTA sequence. Exon features describe both untranslated regions and coding regions, while CDS features only annotate coding regions. Protein annotations are determined only based on CDS and stop codon features.

9. -k, --gff\_reference

Reference FASTA that features in the GFF map to. Typically, GFFs map exons to genomic coordinates. However, this may also be a custom vector-transgene composite reference FASTA if a custom GFF was generated.

10. -t, --targets

Target BED file specifying target regions of the transcript to report variant calls in. Supplying this option only alters final reporting of variants, and does not speed up processing. This is because of the requirement for perfectly paired reads, which may be compromised by intersection of the alignments with the target region prior to variant calling.

11. -d, --consensus\_deduplicate

Flag to turn on consensus deduplication. Use with -u/--umi\_regex to specify a regular expression to match the UMI and anchoring adapter sequence. The UMI will be moved to the read names and any anchoring adapter sequence is discarded (anchoring is recommended but not required). 

UMI extraction is performed by UMI-tools extract prior to adapter trimming with cutadapt.

12. -u, --umi\_regex

Python regex package regular expression for matching the UMI within a desired edit distance and for matching and discarding anchoring adapter sequence.

13. -s, --mutagenesis\_signature

Mutagenesis signature which matches one of the IUPAC DNA codes NNN, NNK, NNS. Candidate variant calls will be tagged with a boolean to annotate a match. No filtering on the signature is performed.

14. -q, --min\_bq

Minimum base quality for either mate of a pair to be considered for variant calling. Default 30. 

15. -m, --min\_supporting

Minimum number of fragments for a candidate variant call. Default 2 (discard singletons).

16. -w, --max\_mnp\_window
Integer window span to search for phased SNPs and call MNPs. Must be between 1 and 3 (default 3). satmut\_utils does not support long-range haplotype calling, which is challenged by exponentially increasing false positive calls with a wider window span.

17. -n, --ntrimmed

cutadapt option (-n) for number of adapters to be trimmed from each read. Default 3. 

Internal PCR tiles normally have two possible adapters whereas terminal PCR tiles may have three. This is because the read emanating from the insert towards the vector in a terminal PCR tile should have a 5' adapter (sequencing adapter) and potentially two 3' adapters (adjacent vector sequence and sequencing adapter).

18. -l, --overlap\_length

cutadapt option (-m) for the min length of matched adapter required for trimming.

This moderates the compromise made by --ntrimmed where three adapters are provided by default, which may cause over-zealous read trimming. As the length increases, adapter trimming becomes more specific but less sensitive. satmut_utils default local alignment should help clip adapters from aligned segments in cases where the adapter is not recognized with a lower min length value.

19. -b, --trim\_bq

cutadapt option (-q) for the length of adapter match required for trimming.

20. --ncores

Number CPU cores to use for cutadapt. Default 0, autodetect.

21. -c, --contig\_del\_threshold

If -z/--race\_like and -cd/--consensus\_deduplicate are provided, convert deletions spanning wider than this threshold to runs of the unknown base N. Required as some R2s may share the same R1 [UMI x POS] but align to non-overlapping coordinates. In other words, consensus deduplication of RACE-like data may generate an unknown segment in the R2 consensus. This allows more accurate reporting of fragment coverage. To avoid this behavior and omit R2 merging from separate amplicons, provide -f/--primer\_fasta, which will annotate read pairs with a unique amplicon/tile.

22. -f, --primer\_fasta
If -z/--race\_like and -cd/--consensus\_deduplicate are provided, reads can be annotated with an originating R2 primer, which prohibits merging of R2s in consensus deduplication. With this option, fragment coverage (DP) may be over-reported in certain regions because of the multi-amplicon coverage of RACE-like data.

23. -a, --primer\_nm\_allowance

If -f/--primer\_fasta, allow up to this number of edit operations for matching primers in the start of R2. Default 3. The last sixteen 3' nucleotides of the matched primer will be appended to the read names to avoid UMI grouping and consensus deduplication. R2s that do not match any primer will be reassigned the unknown primer regex X{16}.

24. --keep\_intermediates

Option to write intermediate files to the output directory. These include preprocessed FASTQ files (trimmed and/or UMI-extracted), alignment files, and log files for preprocessing steps.

## Tests

To run unit tests, execute the following from the satmut_utils repository:

```nose2 -q```

## Accessory scripts

Two command-line interfaces are provided to enable pre-processing of reads prior to  satmut_utils 'sim':

1. satmut\_trim
2. satmut\_align
 
satmut\_trim is a wrapper around cutadapt, and satmut\_align a wrapper around bowtie2.  satmut\_align should be used to generate the BAM file accepted by 'sim'. If reads have been aligned with some other method, there is no guarantee 'sim' will complete without error, as alignment tags output by bowtie2 are required for 'sim' (MD, NM).
