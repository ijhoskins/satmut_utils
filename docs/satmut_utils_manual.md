![](satmut_utils_logo.png)

[Installation](#Installation)

[Reference files](#Reference-files)

[satmut_utils sim](#satmut_utils-sim)

[satmut_utils call](#satmut_utils-call)

[Code examples](#Code-examples)

[sim code examples](#sim-code-examples)

[call code examples](#call-code-examples)

[satmut_utils command-line interface](#satmut_utils-command\-line-interface)

[Common options](#Common-options)

[sim options](#sim-options)

[call options](#call-options)

[Accessory scripts](Accessory-scripts)

[Tests](#Tests)

satmut_utils is a Python package for simulation and variant calling of saturation mutagenesis data. The two main subcommands are:
1. sim
2. call

satmut_utils is designed to simulate and call low-frequency variants (<1:1000) in paired-end, targeted RNA-sequencing reads. Alignments to a single transcript are expected; whole transcriptome or whole genome datasets are currently not supported.

sim and call support two different types of paired-end read chemistries. The first supported chemistry is a tiled amplicon PCR approach with interleaved sets of PCR amplicons. In this mode, R1 and R2 start at primer ends.

The second supported chemistry is a RACE-like approach (for example, Anchored Multiplex PCR), where R1 starts at a variable fragment end and R2 starts at a primer. RACE-like PCR may be used to sequence entire coding regions efficiently, in addition to targeting of RNA 5' and 3' ends.

satmut_utils call allows for unique molecular indices (UMIs) at the start of R1 in either library preparation chemistry. Currently, it does not support analysis of UMIs on both reads.

Furthermore, satmut_utils does not support barcode sequencing, wherein sequencing is used to link a variant with a barcode, followed by readout of unique barcodes.

## Installation

Currently, only Linux and MacOSX operating systems are supported. To get started, follow these steps:

1. If conda is not installed, install miniconda for managing environments. See this [link](https://docs.conda.io/en/latest/miniconda.html) for installation on your particular architecture.

2. Clone the satmut\_utils repository:
```
git clone https://github.com/ijhoskins/satmut_utils.git
```

3. Execute the bash script to generate the satmut\_utils environment, install the package, and optionally download curated reference files, which are required if using Ensembl identifiers ([see Reference files](#Reference-files)).
```
REF_DIR="~/satmut_utils_refs"
satmut_utils/install_satmut_utils.sh -t -g -r $REF_DIR
```

If you already have the GRCh38 FASTA, you may use it provided it uses Ensembl contig nomenclature (not the NCBI "NC_" or UCSC "chr" nomenclature). That is, the chromosome names should start with a single integer for autosomes, and X, Y, and MT for the sex chromosomes and mitochondrial chromosome, respectively. This is important for compatibility with the curated transcript annotations when using Ensembl identifiers. If the FASTA meets these requirements, it should be copied or linked in the REF\_DIR (optionally set by install\_satmut\_utils.sh -r option). The genome should be unzipped and have an index.

```
samtools faidx $REF_DIR/GRCh38.fa
```

You are now ready to call the command-line executable:
```satmut_utils```

satmut\_utils is the primary command, and has the subcommands sim and call.

## Reference files

For convenience, a curated transcriptome of primary human transcripts from [APPRIS](https://apprisws.bioinfo.cnio.es/landing_page/) is provided, allowing the user to pass an Ensembl gene or transcript ID to satmut_utils. However, if the requested Ensembl ID is not found in this set, custom reference files must be passed, which include:

A. Transcript reference (FASTA)

B. Transcript annotations (GFF)

C. GFF reference (FASTA)

Common transcript annotations available in GFF format map exon coordinates in the genome. For this case, file A should specify a transcript FASTA and file C should specify the genome FASTA. This is the default satmut_utils configuration when using Ensembl IDs.

In typical saturation mutagenesis experiments a transgene is expressed from a vector. In this case, the user ideally sets A and C to a custom composite (vector + target) reference FASTA, then makes a custom annotation GFF (B). Custom references are useful for mapping PCR tiles that span the vector-transgene junctions (terminal tiles). Because local alignment is employed in satmut\_utils call, variants near these junctions may be clipped along with vector sequence unless custom reference files are provided and adapters are fully trimmed.

Otherwise, genome-based reference files allow full flexibility for alignment to novel isoforms such as splicing reporters that contain retained introns or 5'/3' extensions that deviate from canonical isoforms, so long as these inclusions to the coding sequence generate a contiguous coding sequence within the transcript reference FASTA (-r).

## satmut_utils sim

sim makes certain design decisions that may impact the performance of sim for your particular application.

1. *sim edits variants in only those fragments with coverage from both mates of the pair*

For compatibility with the call subcommand, which requires mate concordance for a variant call, sim will generate the variants in only those fragments with supporting coverage from both mates. We term this concordant alternate observations (CAO). The denominator to the variant frequency calculation (DP) is also fragment-based. That is, depth is not counted twice at positions where the mates overlap. The resulting frequency, CAF=CAO/DP, is a more conservative measure of variant abundance.

Note that satmut\_utils sim truth frequencies may differ from variants quantification by other callers due to this design choice. Thus, validation of sim results is only recommended with satmut\_utils call, as sim and call were designed to follow the same filtering logic, enabling high-accuracy variant editing.

2. *sim relies on heuristics to select fragments for editing*

To enable editing of multiple variants at a single position, as well as prohibit unintentional phasing of variants, sim employs a heuristic to enforce several rules. The designs assumes that the collection of edited variant frequencies may not exceed 1. Thus, sim assumes that **all variants are configured for editing at low frequencies (<1:1000)**. If the user seeks to edit variants at higher frequencies, the pool of fragments that are amenable for editing rapidly dwindles. This is again because sim prohibits phasing of edited variants.

This constraint is particularly problematic if the user seeks to edit even a few variants at germline-like frequencies (0.5 - 1). sim will raise a InvalidVariantConfig exception if the sum of variant frequencies across all variants in the input VCF exceeds 1. To ignore this exception and edit as many variants as possible, pass --force_edit.

The heuristic for selecting reads enforces the following rules:

1) Any fragment selected for editing has not previously been selected for editing of another variant. This requirement ensures variants edited at nearby positions are never phased, which may lead to simultaneous false negative and false positive calls. Note this applies to variants that are different records in the VCF. To create phased variants (e.g. MNPs), input single variant records with REF and ALT fields that span multiple bases.

2) Any fragment selected for editing must have error-free coverage in both mates within a symmetric buffer spanning the variant position. This prohibits editing a true variant adjacent to an error, which would convert the variant to higher order (for example, SNP to di-nucleotide MNP).

3) If a primer BED file has been provided, sim enforces that no part of the variant span (POS + REF field length) intersects a read segment arising from a synthetic primer. Passing a primer BED file to sim is highly recommended because it enables the avoidance of editing at read termini. Otherwise, under subsequent local alignment of edited reads, variants at read termini may be clipped, leading to their loss from the alignment and causing false negatives.

Note that rule 3 does *not* mean variants cannot be edited within or overlapping primer regions, as long as the position has read-through coverage (typically provided by an adjacent PCR tile). 

Collectively, these rules ensure high fidelity of variant editing at ultra-low frequencies.

## satmut_utils call

call makes certain design decisions that may impact the performance of variant calling for your particular application.

1. *call enforces mate pair concordance for candidate calls*

To reduce the impact of single-stranded sequencing errors, satmut_utils requires the same base call is made in both mates of the pair for enumeration of a variant call. Note this does not handle polymerase errors made during PCR enrichment.

2. *call applies quality filters prior to verifying mate concordance*

Before pairs are considered for mate concordance, reads or base calls may be filtered by read edit distance (--man_nm) and base quality (--min-bq) filters.

3. *call uses a unique algorithm to resolve SNPs and MNPs contained in the same fragment*

The algorithm for variant calling is as follows:

A. Apply quality filters on read edit distance and base quality.

B. Find all concordant mismatched base calls between read mates of a pair.

C. For D-F, ensure no mismatch is part of more than one MNP call.

D. At a mismatch, determine if there is a second mismatch downstream within the window (--max\_mnp\_window).

E. If D), determine if there is third mismatch downstream within the window. If so, call a tri-nucleotide MNP.

F. If D), and the next two mismatches are not adjacent, call a di-nucleotide MNP.

G. After all MNPs within pair have been called, call remaining mismatches as SNPs.

This algorithm leads to specific calling expectations for consecutive mismatch runs/tracts with nearby isolated mismatches. In these cases, consecutive mismatch runs are called as compact MNPs, as opposed to including the isolated mismatch as part of the MNP call.

### Mate pair concordance

Enforcing concordance may lead to under-quantification of variant frequencies in cases where a significant proportion of pairs do not have appreciable overlap between R1 and R2. Insufficient overlap between mates may be due to inadequate primer design and/or overzealous 3’ base quality trimming. Additionally, for RACE-like libraries, a high proportion of fragment lengths greater than the read length may also limit sensitivity.

### satmut_utils frequency calculation

The denominator for the variant frequency calculation- i.e. fragment coverage depth- is determined after filtering on read edit distance (NM tag) and mate pair overlap (concordance). If a primer BED file is provided to satmut_utils, read segments originating from synthetic primer sequence do not contribute to fragment depth (DP).
When the MNP span (variant reference coordinates) covers both a synthetic primer position and an adjacent position in the amplicon, the minimum depth of coverage in the MNP span is used as the denominator for the frequency: CAF=CAO/DP.

### satmut\_utils call error correction preprocessing steps

satmut\_utils call supports multiple methods for moderating false positives calls. Application of these methods is optional but improves specificity.

1. Synthetic primer base quality masking

To prohibit false positive calls arising from primer synthesis errors, satmut_utils call provides a masking strategy to demarcate read segments originating from synthetic primer sequence. In masking, base qualities for these segments are set to 0 and omitted from variant calls. This step requires the user to provide a primer BED file specifying the primers used in target enrichment.

2. Unique-molecular-identifier (UMI)-based read consensus deduplication

UMIs tag unique molecules prior to target enrichment with PCR, allowing for a more accurate estimation of molecular counts. UMI-based consensus deduplication improves the specificity of calls and variant frequency accuracy, as PCR jackpots are deduplicated/collapsed and some PCR and sequencing errors may be removed in the process of generating the consensus read.

This option assumes UMIs are at the start of R1, and utilizes UMI-tools (Smith et al. 2017) to group read duplicates based on R1 UMI-POS, where POS is the R1 aligned start position. In grouping pairs, the template length is ignored by passing --ignore-tlen. satmut\_utils then employs a consensus deduplication workflow to correct errors within duplicate groups and generate a consensus R1 and R2. The consensus base call at each aligned position of the duplicate group is determined by majority vote. In ties (two duplicates in the UMI group), if one read matches the reference, the reference base is chosen. Otherwise, the base call with the higher base quality is selected. If both bases have the same quality, the consensus base is chosen at random.

Unlike other consensus deduplication tools (Clement et al. 2018; Chen et al. 2019), consensus deduplication is supported for RACE-like libraries. In one RACE-like library preparation method (Anchored Multiplex PCR, Zheng et al. 2014), the R2 starts at a primer and R1 contains a UMI-adapter which ligates to an internal fragmentation site or the 5’ or 3’ end of the molecule.

For RACE-like libraries, as a result of passing umi\_tools the --ignore-tlen option, R2s that share the same R1 UMI-POS but do not share the same R2 start coordinate may be merged into a R2 consensus contig. This logic ensures accurate fragment depth of coverage reporting for RACE-like data, but may create R2 contigs with read lengths greater than the sequencing read length. It may also lead to a small proportion of R2 contigs with internal unknown base calls (N), which arise when merged R2s do not overlap one another.

To disable such merging of RACE-like R2s during consensus deduplication, the user can provide --primer_fasta. When a primer FASTA is provided, the primer from which R2 originates is appended to the UMI prior to grouping. This ensures R2s emanating from different primers are assigned to separate groups, despite sharing the same R1 UMI-POS. Note that fragment depth of coverage will not be as accurate with this option. However, this option prohibits deletion artifacts arising from merging of non-overlapping R2s. See also the --contig\_del\_threshold option.

To recap, when provided RACE-like data, satmut\_utils assembles R2 contigs from multiple R2s sharing a common R1 UMI-position. While not supported, dual UMIs on R1 and R2 may be possible by modification of the source code that calls umi\_tools. 

## Code examples

Parameter help:
```
satmut_utils -h
satmut_utils sim -h
satmut_utils call -h
```

Common arguments to both sim and call subcommands should be provided first, then the subcommand, and then the subcommand-specific arguments.

It is recommended that a new output directory is created for each job. Default is to output to the current directory.

```OUTPUT_DIR="/tmp/satmut_utils_test"```

### sim code examples

Run the sim workflow by providing a BAM containing paired-end, single-contig alignments and a VCF file specifying variants and their desired frequencies. 

```
satmut_utils -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR sim -a negative_control.bam -v sim_variants.vcf
```

The BAM file should be generated with satmut\_align to ensure required alignment tags (e.g. MD tag) are present.

Typically the alignments are from a non-mutagenized negative control library and the variants are samples from a VCF containing all codon permutations matching the mutagenesis signature across the target region. See run\_variant_generatory.py to generate test VCFs.

See src/tests/test\_data/cbs\_sim.vcf for an example VCF.

Provide a primer BED file for primer-aware read editing.
```
satmut_utils -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p primers.bed sim -a negative_control.bam -v sim_variants.vcf
```

See sim requirement 4 in the section above for a description of primer masking in the context of simulation. The primer BED must have a strand field. See tests/test\_data/CBS\_insilico_primers.bed for an example BED.

Specify that the alignments are from a RACE-like library preparation chemistry, and use a new random seed to select different read pairs for editing:
```
satmut_utils -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p primers.bed --race_like sim -a negative_control.bam -v sim_variants.vcf -s 99
```

### sim outputs

The sim workflow outputs paired FASTQs, a realigned BAM file, and a truth VCF containing expected variants and their quantification. Edited reads that do not realign will not give rise to variant calls in downstream analysis.

### call code examples

Run call by specifying an Ensembl transcript ID and the directory containing curated reference files:
```
satmut_utils -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5 TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT
```

-5 and -3 specify adapter sequences for trimming read pairs *relative to the R1*. Equivalently, for R2, its 5' adapter is provided to -3 and its 3' adapter to -5.

More than one 5' adapter and more than one 3' adapter are often needed to additionally trim vector sequences (e.g. attB sites) from reads of terminal PCR tiles that span the vector-transgene junctions. In this case, provide multiple comma-delimited adapters:

```
-5 TACACGACGCTCTTCCGATCT,CAAGTTTGTACAAAAAAGTTGGC -3 AGATCGGAAGAGCACACGTCT,CCAACTTTCTTGTACAAAGTGGT
```

The Ensembl ID may also specify a gene:
```
satmut_utils -i ENSG00000160200.17 -x $REF_DIR -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5 TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT
```

For Ensembl gene ID with more than one transcript isoform, satmut\_utils will select the primary transcript in the curated APPRIS transcriptome.

If the Ensembl ID is not in the curated set of primary transcripts, the user must provide custom reference files:
```
satmut_utils -r CBS.fa -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5 TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT -g CBS.gtf -k GRCh38.fa
```

Additional files may be passed, such as a primer and target BED file:
```
satmut_utils -i ENST00000398165.7 -x $REF_DIR -p primers.bed -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5 TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT -t target.bed
```

The primer BED must have a strand field. See src/tests/test\_data/CBS\_insilico_primers.bed for an example BED.

The target BED only impacts reporting of variants and does not speed up analysis by intersection of alignments with the target.

### call outputs

The call workflow produces a VCF of candidate variant calls, a tab-delimited summary file, and a bedgraph file reporting fragment coverage across the reference.

The output VCF and its corresponding summary.txt file contain records for each mismatched base in an MNP, so that quality information for the mismatches can be used for machine learning-based error correction. See the VCF header for column/feature descriptions.

A number of useful R functions exist in prototype.summarization_utils.r for parsing and summarizing the VCF summary file.

## satmut_utils command line interface

satmut_utils provides the sim and call workflow as subcommands, which have common and unique options.

### Common options

1. -i, --ensembl_id

An Ensembl gene (ENSG) or transcript (ENST) identifier containing the minor version number. For example, ENST00000398165.7.

2. -r, --reference

A custom reference FASTA for a single transcript or target transgene. Will be indexed by samtools and bowtie2 if not already. This option is mutually exclusive with --ensembl_id.

3. -x, --reference\_dir

If using --ensembl_id, a reference directory containing curated a transcriptome FASTA, transcript annotations (GFF), and the genome FASTA.

4. -z, --race\_like

Flag to indicate input reads or alignments are from a rapid-amplification-of-cDNA-ends (RACE)-like target enrichment strategy (e.g. Anchored Multiplex PCR). Default is tiled PCR enrichment strategy with R1 and R2 both starting at the 5' end of opposing primers.

5. -p, --primers

Primer BED file with six fields, containing the contig, start, stop, name, score, and strand fields. BED files have a 0-based start position.

6. -o, --output_dir

Optional output directory to write results to. Default is current working directory.

7. -j, --nthreads

Number of CPU cores (cutadapt), additional alignment threads (bowtie2), and threads to use for BAM sorting (samtools). Multiprocessing is not currently supported for satmut\_utils sim or call Python code.

8. -e, --max_nm

Maximum edit distance for either mate of a pair to be considered for simulation and variant calling. Default 10.

### sim options

1. -a, --alignments

This is a BAM file containing paired-end alignments to a single transcript or contig. If a BAM file is not generated, use satmut\_align for bowtie2 alignment. If reads have been aligned with another aligner, there is no guarantee sim will complete without error, as certain alignment tags (e.g. MD tag) are required for sim. 

2. -v, --vcf

Variant Call Format file that contains an in-line INFO tag (AF) specifying the desired variant frequency. See src/tests/test\_data/CBS\_sim.vcf for an example.

3. -b, --edit_buffer

When selecting reads to edit, verify the absence of errors within this buffer about the edit span (range of reference coordinates having mismatches). The read segment must match the reference within the span (POS - edit_buffer, POS + REF length + edit_buffer). Increasing this value prevents variant conversion and/or unexpected clipping of the variant from the termini of the read under local alignment. Decreasing this value may be needed in cases where variants are being edited near read termini and there are no adjacent PCR tiles to provide coverage.

4. -f, --force_edit

By default, sim will raise a InvalidVariantConfig Exception if the sum of variant frequencies across all variants in the input VCF exceed 1. Invalid configurations indicate that all variants may not be edited due to insufficient read depth of coverage. In this case, to edit as many variants as possible, provide this flag.

InvalidVariantConfig exceptions are guaranteed to result in unedited variants if alignments from a single PCR tile are being edited. However, if variants are to be edited into alignments spanning multiple PCR tiles, all variants *may* be edited despite an InvalidVariantConfig exception. This flag is particularly useful for maximizing the number of variants edited for multi-tile alignments.

5. -y, --random_seed

Integer seed to use for pseudorandom qname (read name) sampling. This may be used to select different reads for editing of variants.

### call options

1. -1, --fastq1

R1 FASTQ of the pair. If using simulated data that has previously been trimmed, pass the -v/--omit trim option.

If -p (--primers) is provided, the read names must match either Illumina format or consist of a unique integer to facilitate primer BQ masking.

2. -2, --fastq2

R2 FASTQ of the pair.

3. -5 --r1\_fiveprime_adapters

Comma-delimited R1 5' adapters, or None if no R1 5' adapters exist. See also -v option.

4. -3 --r1\_threeprime_adapters

Comma-delimited R1 3' adapters, or None if no R1 3' adapters exist. See also -v option.

5. -5 --r2\_fiveprime_adapters

Comma-delimited R2 5' adapters, or None if no R2 5' adapters exist. See also -v option.

6. -3 --r2\_threeprime_adapters

Comma-delimited R2 3' adapters, or None if no R2 3' adapters exist. See also -v option.

7. -g, --transcript_gff

Transcript GFF where **features are ordered from 5' to 3', regardless of strand**. For examples, see src/tests/test\_data/gencode.v29.annotation.gtf for a standard genome-based GFF example, or src/tests/test\_data/CBS_pEZY3.gff for a custom, composite vector GFF example.

IMPORTANT: the GFF seqname field and a "transcript_id" attribute in the attribute field should match the contig name in the reference FASTA (-r/--reference).

The minimum expected GFF records have a feature field (3rd column of the GFF) specifying one of the following:
```
exon, CDS, stop_codon
```

Additional allowable features in the feature field (3rd column) of the GFF include:
```
gene, transcript, start_codon
```

While these latter features/records are not necessary, they are often helpful for thorough annotation of the transcript. While a start\_codon feature is contained within the most 5' CDS feature, a stop_codon feature is *not* contained in the most 3' CDS feature, and must be provided if the user seeks to call variants altering the stop codon (nonstop variants). 

Otherwise, one is free to specify exon and CDS features that annotate the coding and noncoding portions of each transcript, provided the exon features sum up to the full reference FASTA sequence. Exon features describe both untranslated regions and coding regions, while CDS features only annotate coding regions. Protein annotations are determined only based on CDS and stop codon features.

8. -k, --gff_reference

Reference FASTA that features in the GFF map to. Typically, GFFs map exons to genomic coordinates. However, this may also be a custom vector-transgene composite reference FASTA if a custom GFF was generated.

9. -t, --targets

Target BED file specifying target regions of the transcript to report variant calls in. Supplying this option only alters final reporting of variants, and does not speed up processing. This is because of the requirement for perfectly paired reads, which may be compromised by intersection of the alignments with the target region prior to variant calling.

10. -d, --consensus_deduplicate

Flag to turn on consensus deduplication. Use with -u/--umi_regex to specify a regular expression to match the UMI and anchoring adapter sequence. The UMI will be moved to the read names and any anchoring adapter sequence is discarded (anchoring is recommended but not required). 

UMI extraction is performed by umi_tools extract prior to adapter trimming with cutadapt.

11. -u, --umi_regex

Python regex package regular expression for matching the UMI within a desired edit distance and for matching and discarding anchoring adapter sequence.

12. -s --mutagenesis_signature

Mutagenesis signature which matches one of the IUPAC DNA codes NNN, NNK, NNS. Candidate variant calls will be tagged with a boolean to annotate a match. No filtering on the signature is performed. Instead, variants are annotated for signature match.

13. -q, --min_bq

Minimum base quality for either mate of a pair to be considered for variant calling. Default 30. 

14. -m, --min_supporting

Minimum number of fragments for a candidate variant call. Default 2 (discard singletons).

15. -w, --max\_mnp\_window
Integer window span to search for phased SNPs and call MNPs. Must be 2 or 3 (default 3). satmut_utils does not support long-range haplotype calling, which is challenged by exponentially increasing false positive calls with a wider window span.

16. -n, --ntrimmed

cutadapt option (-n) for number of adapters to be trimmed from each read. Default 3. 

Internal PCR tiles normally have two possible adapters whereas terminal PCR tiles may have three. This is because the read emanating from the insert towards the vector in a terminal PCR tile should have a 5' adapter (sequencing adapter) and potentially two 3' adapters (adjacent vector sequence and sequencing adapter).

17. -l, --overlap_length

cutadapt option (-m) for the min length of matched adapter required for trimming.

This moderates the compromise made by --ntrimmed where three adapters are provided by default, which may cause over-zealous read trimming. As the length increases, adapter trimming becomes more specific but less sensitive. satmut_utils default local alignment should help clip adapters from aligned segments in cases where the adapter is not recognized with a lower min length value.

18. -b, --trim_bq

cutadapt option (-q) for the length of adapter match required for trimming.

19. -v, --omit_trim

Useful for input where adapters have been previously trimmed. Note that trimming when no adapters are present may degrade the data quality by nonspecific trimming of segments internal and terminal to the reads. For this reason, it is best to directly re-align these inputs.

20. -c, --contig\_del_threshold

If -z/--race\_like and -cd/--consensus\_deduplicate are provided, convert deletions spanning wider than this threshold to runs of the unknown base N. Required as some R2s may share the same R1 [UMI x POS] but align to non-overlapping coordinates. In other words, consensus deduplication of RACE-like data may generate an unknown segment in the R2 consensus. This allows more accurate reporting of fragment coverage. To avoid this behavior and omit R2 merging from separate amplicons, provide -f/--primer_fasta, which will annotate read pairs with a unique amplicon/tile.

21. -f, --primer_fasta
If -z/--race\_like and -cd/--consensus\_deduplicate are provided, reads can be annotated with an originating R2 primer, which prohibits merging of R2s in consensus deduplication. With this option, fragment coverage (DP) may be over-reported in certain regions because of the multi-amplicon coverage of RACE-like data.

22. -a, --primer\_nm_allowance

If -f/--primer_fasta, allow up to this number of edit operations for matching primers in the start of R2. Default 3. The last sixteen 3' nucleotides of the matched primer will be appended to the read names to avoid UMI grouping and consensus deduplication. R2s that do not match any primer will be reassigned the unknown primer regex X{16}.

23. --keep\_intermediates

Option to write intermediate files to the output directory. These include preprocessed FASTQ files (trimmed and/or UMI-extracted), alignment files, and log files for preprocessing steps.


## Tests

To run unit tests, execute the following from the satmut_utils repository:

```nose2 -q```

## Accessory command-line interfaces and scripts

Two command-line interfaces are provided to enable pre-processing of reads prior to  satmut_utils sim:

1. satmut\_trim
2. satmut\_align
 
satmut\_trim is a wrapper around cutadapt, and satmut\_align a wrapper around bowtie2.  satmut\_align should be used to generate the BAM file accepted by sim. If reads have been aligned with some other method, there is no guarantee sim will complete without error, as alignment tags output by bowtie2 are required for sim (MD, NM).

To facilitate simulation of reads and variants in a desired transcript *de novo*, various accessory scripts are provided in the src/scripts directory. Code here is not fully tested and is only provided for convenience to facilitate error modeling.

1. run\_bowtie2\_aligner.py.
This script may be called directly, but for convenience satmut\_utils provides the CLI satmut\_align.

```
satmut_align -f1 satmut_utils/src/tests/test_data/CBS_sim.R1.fq.gz -f2 satmut_utils/src/tests/test_data/CBS_sim.R2.fq.gz -r satmut_utils/src/tests/test_data/CBS.fa -d /tmp/align_example
```

2. run\_variant\_generator.py
This script may be used to generate a VCF of all SNP and MNP permutations matching a given mutagenesis signature in desired transcript coding regions. Together with run\_read\_generator.py, this provides *de novo* inputs for satmut\_utils sim. Provide --trx\_id exactly as it exists in the reference FASTA, target BED, and annotation GFF. Often these files can be easily obtained by first running satmut\_utils call on a negative control library and using the outputted reference files.

```
python -m scripts.run_variant_generator -i CBS_pEZY3 -s NNK -r satmut_utils/src/tests/test_data/CBS_pEZY3.fa -g satmut_utils/src/tests/test_data/CBS_pEZY3.gff -t satmut_utils/src/tests/test_data/CBS_pEZY3_targets.bed -d /tmp/var_gen_example
```

3. run\_vcf\_subsampler.py
This script can be used to subsample variants or mismatched bases from the VCF produced by run\_variant\_generator.py. Balancing true and false positive variants is highly recommended for training error models.

```
python -m scripts.run_vcf_subsampler -v satmut_utils/src/tests/test_data/CBS_sim.vcf -n 2 -d /tmp/vcf_subsamp_example
```

4. run\_ec\_data\_generator.py
This script generates simulated datasets modeled after true data. As input, it requires satmut\_utils call output summary.txt files for a true mutagenized library and a non-mutagenized negative control library. It then generates codon-permutation variants, subsamples them to balance true and false positives, and configures variant frequencies by estimating parameters for SNPs and MNPs using variants the mutagenized summary.txt file (optionally for those variants only matching the mutagenesis signature). It finally invokes satmut\_utils sim to generate the simulated dataset.

satmut\_utils call should be ran thereafter (ideally with loose quality parameters: -m 1 -q 1 -e 150) to extract quality features and complete the validation dataset. R utilities for training machine learning models on the resulting calls are provided in src/prototype/modeling_utils.R.

5. run\_read\_generator.py.

This may be used to simulate paired-end targeted RNA-sequencing reads given a reference FASTA. Useful for generating input reads for sim. However, one of the many Next-Generation Sequencing read simulators that construct error models are recommended to generate more realistic test reads.

To generate error-free, tiled RNA read alignments, provide a transcript reference FASTA and target BED file. Split the transcript target region into chunks with span roughly [read length - 2\*(average primer length)] if the target region requires multiple tiles. Then run with --make\_amplicons but *not* --rna. (Use of --rna is for simulation of random RNA fragments and when starting with a genomic reference FASTA).

--make\_amplicons is intended for direct simulation of reads from a transcript FASTA. For general DNA or RNA read generation using a genome FASTA and standard GFF annotations, omit --make\_amplicons. In this mode, generated fragments start and end at random coordinates in the sequence space informed by the --frag\_length argument.

As an example, for 2 x 150 bp chemistry, create ~100 bp interleaved target chunks in BED format and configure the number of reads to generate for each amplicon by the BED score field. Finally, pass --make\_amplicons and --slop_length 0 so that reads start at the termini of the targets.

