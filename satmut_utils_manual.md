# satmut_utils:

satmut_utils is a Python package for simulation and variant calling of saturation mutagenesis data. The two main subcommands are:
1. sim
2. call

satmut_utils is designed to simulate and call low-frequency variants (<1:1000) in paired-end, targeted RNA-sequencing datasets. Alignments to a single transcript are expected; whole transcriptome or whole genome datasets are currently not supported.

sim and call support two different types of paired-end read chemistries. The first supported chemistry is a tiled PCR approach with interleaved sets of PCR amplicons. In this mode, R1 and R2 start at primer ends.

The second supported chemistry is a RACE-like PCR approach (for example, Anchored Multiplex PCR), where R1 starts at a variable fragment end and R2 starts at a primer.

call allows for unique molecular indices (UMIs) at the start of R1 in either library preparation chemistry. Currently, it does not support analysis of UMIs on both reads.

Furthermore, satmut_utils does not support "barcode sequencing", wherein sequencing is  used to link a variant with a barcode, followed by readout of the linked barcodes.


## Reference files

For convenience, a curated transcriptome of primary transcripts from APPRIS is provided, allowing the user to pass an Ensembl gene or transcript identifier to satmut_utils sim and call. However, if the requested Ensembl ID is not found in this set, the user must pass their own reference files. These include:

A. Transcript reference (FASTA)

B. Transcript annotations (GFF)

C. GFF reference (FASTA)

Common transcript annotations available in GFF format map exon coordinates in the genome. For this case, file A should specify a transcript FASTA and file C should specify the genome FASTA. This is the default satmut_utils configuration when using Ensembl IDs.

However, in typical saturation mutagenesis experiments a transgene is expressed from a vector.Therefore, I recommend that the user set A and C to the same composite (vector + target) reference FASTA, then make a custom GFF annotation (B) that denotes the coding regions within the composite reference. This is primarily useful for mapping of PCR tiles that span the vector-transgene junctions. Because local alignment is used, variants near these junctions may be clipped along with vector sequence unless custom reference files are provided.

Finally, note that genome-based reference files for B and C allows full flexibility for alignment to novel isoforms such as splicing reporters that contain  retained introns or 5'/3' extensions that deviate from canonical isoforms.


## Setting up the running environment

To get started, you must clone the satmut_utils repository, create a conda environment, and download or generate reference files. This process is as follows:

1. Acquire the code base by cloning the satmut_utils repository:

```
git clone https://github.com/ijhoskins/satmut_utils.git
```

2. Install miniconda for managing environments and packages. See this [link](https://docs.conda.io/en/latest/miniconda.html) for installation.

3. Create the conda environment:
```
cd satmut_utils && conda env create -f satmut_utils_env.yaml
```

4. Activate the environment
```
conda activate satmut_utils
```

5. Download reference files (if using Ensembl identifiers).

Create a reference directory:
```
REF_DIR="/path_to_refs"
mkdir $REF_DIR
```

Download a curated set of primary transcripts and unzip the files:
```
git clone https://github.com/ijhoskins/satmut_utils_refs.git
cp satmut_utils_refs/* $REF_DIR
gunzip $REF_DIR/*gz
```

Download the [human genome FASTA](https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz) and move it to the reference directory. Then index it with samtools:
```
samtools faidx $REF_DIR/GRCh38.fa
```

Other genome reference files for other organisms are available [here](http://daehwankimlab.github.io/hisat2/download/)

6. Navigate to the calling directory within the repository:
```
cd satmut_utils
```


## Code help

Parameter help:
```
python satmut_utils.py -h
python satmut_utils.py sim -h
python satmut_utils.py call -h
```

It is recommended that the user create a new output directory for each job.
```OUTPUT_DIR="/path_to_output_folder"```

Common arguments to both sim and call subcommands are provided first, then the subcommand, and then subcommand-specific arguments.


## satmut_utils sim

sim makes certain design decisions that may impact the performance of sim for your particular application.

1. *sim edits variants in only those fragments with coverage from both mates of the pair*

For compatibility with the call subcommand, which requires mate concordance for a variant call, sim will generate the variants in only those fragments with supporting coverage from both mates. sim will not filter out read pairs that do not overlap at the variant position. Thus, when using satmut\_utils sim with another variant caller, observed frequencies may differ from configuration depending on the logic of variant calling. As such, for validation I recommend the use of satmut\_utils call on sim outputs. 

2. *sim relies on heuristics to select fragments for editing*

To enable editing of multiple variants at a single position and at nearby positions, sim employs heuristics to enforce several rules. A side-effect is that **variants must be edited at low frequencies (<1:1000)**.

If the user seeks to edit variants at higher frequencies, the pool of fragments that are amenable for editing rapidly dwindles. This is particularly problematic if the user seeks to edit a handful of variants at germline-like frequencies. sim will issue a warning if the sum of configured frequencies across all variants and positions exceeds 1 - error rate.

1) Any fragment selected for editing has not previously been selected for editing of another variant. This requirement ensures variants edited at nearby positions are never phased, which may lead to simultaneous false negative and false positive calls. Note this applies to variants that are different records in the VCF.

2) Any fragment selected for editing must have error-free coverage in both mates at the variant position. This ensures native errors in the dataset are preserved.

3) Any fragment selected for editing must have error-free coverage in both mates within a window spanning the variant position. This prohibits editing a true variant adjacent to an error, which would convert the variant to higher order (for example, SNP to di-nucleotide MNP).

4) If a primer feature file has been provided, sim enforces that no part of the variant span intersects a synthetic primer sequence. This does *not* mean variants cannot be edited within or overlapping primer regions. It instead requires read-through coverage for such editing, which is typically provided by an adjacent, interleaved PCR amplicon.

### sim code examples

Run the sim workflow by providing a BAM containing paired-end, single-contig alignments and a VCF file specifying variants and their desired frequencies. 

```
python satmut_utils.py -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR sim -a negative_control.bam -v sim_variants.vcf
```

Typically the alignments are from a non-mutagenized negative control library and the variants are samples from a VCF containing all codon permutations matching the mutagenesis signature across the target region. See run\_variant_generatory.py to generate test VCFs.

See cbs\_sim.vcf in the tests/test_data directory for an example VCF.

Provide a primer BED file for primer-aware read editing:
```
python satmut_utils.py -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p primers.bed sim -a negative_control.bam -v sim_variants.vcf
```
See sim requirement 4 in the section above for a description of primer masking in the context of simulation. The primer BED must have a strand field. See tests/test\_data/CBS\_insilico_primers.bed for an example BED.

Specify that the alignments are from a RACE-like library preparation chemistry, and use a random seed to select different reads for editing:
```
python satmut_utils.py -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p primers.bed --race_like sim -a negative_control.bam -v sim_variants.vcf -s 99
```

### sim outputs

The sim workflow outputs paired FASTQs, a realigned BAM file, and a truth VCF containing expected variants and their frequencies. Edited reads that do not realign will not give rise to variant calls in downstream analysis.


## satmut_utils call

Generation of custom reference files are recommended for typical saturation mutagenesis experiments where a transgene is expressed from a vector. This is particularly important if target regions span vector-transgene junctions. However, to facilitate ease-of-use, curated primary transcript sequences are available as an alignment reference.

call makes certain design decisions that may impact the performance of variant calling for your particular application.

1. *call enforces mate pair concordance for candidate calls*

To reduce the impact of single-stranded sequencing errors, satmut_utils requires the same base call is made in both mates of the pair for enumeration of a variant call. Note this does not handle polymerase errors made during PCR enrichment.

2. *call applies quality filters prior to verifying mate concordance*

Before pairs are considered for mate concordance, reads or base calls may be filtered by read edit distance (--man_nm) and base quality (--min-bq) filters.

3. *call uses a unique algorithm to resolve SNPs and MNPs contained in the same fragment*

The algorithm for variant calling is as follows:

A. Apply quality filters.

B. Find all concordant base calls between read mates of a pair.

C. For D-F, ensure no mismatch is part of more than one MNP call.

D. At a mismatch, determine if there is a second mismatch downstream within the window (--max_mnp_window).

E. If D), determine if there is third mismatch downstream within the window. If so, call a tri-nucleotide MNP.

F. If D), and the next two mismatches are not adjacent, call a di-nucleotide MNP.

G. After all MNPs within pair have been called, remaining mismatches are called as SNPs.


This algorithm leads to specific calling expectations for consecutive mismatch runs/tracts with nearby isolated mismatches. In these cases, consecutive mismatch runs are called as compact MNPs as opposed to including the isolated mismatch as part of the MNP call.


### call code examples

Run call by specifying an Ensembl transcript ID and the directory containing curated reference files:
```
python satmut_utils.py -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5 TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT
```

-5 and -3 specify adapter sequences for trimming read pairs *relative to the R1*. Equivalently, for R2, its 5' adapter is provided to -3 and its 3' adapter to -5.

More than one 5' adapter and more than one 3' adapter are often needed to additionally trim vector sequences (e.g. attB sites) from reads of terminal PCR tiles that span the vector-transgene junctions. In this case, provide multiple comma-delimited adapters:

```
-5 TACACGACGCTCTTCCGATCT,CAAGTTTGTACAAAAAAGTTGGC -3 AGATCGGAAGAGCACACGTCT,CCAACTTTCTTGTACAAAGTGGT
```

The Ensembl ID may also specify a gene:
```
python satmut_utils.py -i ENSG00000160200.17 -x $REF_DIR -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5 TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT
```

For Ensembl gene ID with more than one transcript isoform, satmut_utils will select the primary transcript in the curated APPRIS transcriptome.

If the Ensembl ID is not in the curated set of primary transcripts, the user must provide custom reference files:
```
python satmut_utils.py -r CBS.fa -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5 TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT -g CBS.gtf -k GRCh38.fa
```

Additional files may be passed, such as a primer and target BED file:
```
python satmut_utils.py -i ENST00000398165.7 -x $REF_DIR -p primers.bed -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5 TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT -t target.bed
```

The primer BED must have a strand field. See tests/test\_data/CBS\_insilico_primers.bed for an example BED.


### call outputs

The call workflow produces a VCF of candidate variant calls as well as a BED file reporting fragment coverage across the reference.

The output VCF and its corresponding tab-delimited summary.txt file contain records for each mismatched base in an MNP, so that quality information for the mismatches can be used for learning-based error correction.

A number of useful R functions exist in prototype.summarization_utils.r for parsing and summarizing the VCF summary file. See the VCF header for column/feature descriptions.


## satmut_utils command-line interface

### Common options

1. -i, --ensembl_id

An Ensembl gene (ENSG) or transcript (ENST) identifier containing the minor version number. For example, ENST00000398165.7.

2. -r, --reference

A custom reference FASTA for a single transcript or target transgene. Will be indexed by samtools and bowtie2 if not already. This option is mutually exclusive with --ensembl_id.

3. -x, --reference\_dir

If using --ensembl_id, a reference directory containing curated a transcriptome FASTA, transcript annotations (GFF), and the genome FASTA.

4. -z, --race\_like

Flag to indicate input reads or alignments are from a rapid-amplification-of-cDNA-ends based target enrichment strategy (e.g. Anchored Multiplex PCR). Default is tiled PCR enrichment strategy with R1 and R2 both starting at the 5' end of opposing primers.

5. -p, --primers

Primer BED file with six fields, containing the contig, start, stop, name, score, and strand fields. BED files have a 0-based start position.

6. -o, --outdir

Optional output directory to write results to. Default is current working directory.

7. -j, --nthreads

Number of CPU cores (cutadapt), additional alignment threads (bowtie2), and threads to use for BAM sorting (samtools). Multiprocessing is not currently supported for satmut_utils.


### sim options

1. -a, --alignments

This is a BAM file containing paired-end alignments to a single contig. If a BAM file is not generated, use scripts/run\_bowtie2\_aligner.py for bowtie2 alignment. Alternatively, intermediate paired BAM files from a satmut_utils call may be used.

2. -v, --vcf

Variant Call Format file that contains an in-line INFO tag (AF) specifying the desired variant frequency.

3. -s, --random_seed

Integer seed to use for pseudorandom qname (read name) sampling. This may be used to select different reads for editing of variants. Default 9.


### call options

1. -1, --fastq1

R1 FASTQ of the pair. If using simulated data that has previously been trimmed, pass the -v/--omit trim option.

If -p (--primers) is provided, the read names must match either Illumina format or consist of a unique integer to facilitate primer BQ masking.

2. -2, --fastq2

R2 FASTQ of the pair.

3. -5 --r1\_fiveprime_adapters

Comma-delimited 5' adapters relative to the R1. The reverse complement of these adapters will be trimmed from the 3' end of R2.

4. -3 --r1\_threeprime_adapters

Comma-delimited 3' adapters relative to the R1. The reverse complement of these adapters will be trimmed from the 5' end of R2.

5. -g, --transcript_gff

Transcript GFF where **features are ordered from 5' to 3', regardless of strand**. For examples, see tests/test\_data/gencode.v29.annotation.gtf for a standard genome-based GFF, or tests/test\_data/CBS_pEZY3.gff for a custom composite vector GFF.

The GFF seqname field and a "transcript_id" attribute in the attribute field should specify the contig name of the custom reference (-r/--reference).

The expected features in the feature field (3rd column) of the GFF include:
```
transcript, start_codon, exon, CDS, stop_codon
```

While the start\_codon feature intersects the first CDS feature, stop_codon is not contained in the last CDS feature, so it should be provided if the user seeks to call variants changing the stop codon. Otherwise, the user is free to specify intervening exon and CDS features that annotate the coding and noncoding portions of each transcript. Protein annotations are determined based on CDS and stop codon features. For both canonical isoforms and novel isoforms (e.g. with retained noncoding exons), exon features describe both untranslated regions (UTRs) and coding exons, and should add up to the  reference FASTA.

6. -k, --gff_reference

Reference FASTA that features in the GFF map to. Typically, GFFs map exons to genomic coordinates. However, this may also be a custom vector-transgene composite reference FASTA.

7. -t, --targets

Target BED file specifying target regions of the transcript to report variant calls in. Supplying this option only alters reporting of variants, and does not speed up processing. This is because 

8. -d, --consensus_deduplicate

Flag to turn on consensus deduplication. Use with -u/--umi_regex to specify a regular expression to match the UMI and anchoring adapter sequence. The UMI will be moved to the read names and any anchoring adapter sequence may be discarded (recommended). 

UMI extraction and anchoring adapter trimming are performed by umi_tools extract prior to processing the FASTQs with cutadapt.

9. -u, --umi_regex

Python regex package regular expression for matching the UMI within a desired edit distance and for matching and discarding anchoring adapter sequence.

10. -s --mutagenesis_signature

Mutagenesis signature which matches one of the IUPAC DNA codes NNN, NNK, NNS. Candidate variant calls will be tagged with a boolean to annotate a match. Specificity may be improved by requiring a match; however, by default not filtering on the signature is performed.

11. -q, --min_bq

Minimum base quality for either mate of a pair to be considered for variant calling. Default 30. 

12. -e, --max_nm

Maximum edit distance for either mate of a pair to be considered for variant calling. Default 7.

13. -m, --min_supporting

Minimum number of fragments for a candidate variant call. Default 2 (discard singletons).

14. -w, --max\_mnp\_window
Integer window span to search for phased SNPs and call MNPs. Must be 2 or 3 (default 3). satmut_utils does not support long-range haplotype calling, which is challenged by exponentially increasing false positive calls with a wider window span.

15. -n, --ntrimmed

cutadapt option (-n) for number of adapters to be trimmed from each read. Default 3. 

Internal PCR tiles normally have two possible adapters whereas terminal PCR tiles may have three. This is because the read emanating from the insert towards the vector in a terminal PCR tile should have a 5' adapter (sequencing adapter) and possibly two 3' adapters (adjacent vector sequence and sequencing adapter).

16. -l, --overlap_length

cutadapt option (-m) for the min length of matched adapter required for trimming.

This moderates the compromise made by --ntrimmed where three adapters are provided by default, which may cause over-zealous read trimming. As the length increases, adapter trimming becomes more specific but less sensitive. satmut_utils default local alignment should help clip adapters from aligned segments in cases where the adapter is not recognized with a lower min length value.

17. -b, --trim_bq

cutadapt option (-q) for the length of adapter match required for trimming.

18. -v, --omit_trim

Useful for satmut_utils sim data where adapters have been previously trimmed by the call workflow. Note that trimming when no adapters are present may degrade the data quality by nonspecific trimming of insert sequences contained in the trimmed reads. For this reason, 
it is best to directly realign simulated datasets.

19. -c, --contig\_del_threshold

If -z/--race\_like and -cd/--consensus\_deduplicate are provided, convert deletions spanning wider than this threshold to runs of the unknown base N. Required as some R2s may share the same R1 UMI-position but align to non-overlapping coordinates. In other words, consensus deduplication of RACE-like data may generate an unknown segment in the R2 consensus read. This allows more accurate reporting of fragment coverage, which reports more accurate variant frequencies.

To avoid this behavior and omit R2 merging from non-overlapping coordinates, provide -f/--primer_fasta, which will annotate read pairs with a unique R2 primer.

20. -f, --primer_fasta
If -z/--race\_like and -cd/--consensus\_deduplicate are provided, reads can be annotated with an originating R2 primer, which prohibits merging of R2s in consensus deduplication. With this option, fragment coverage (the denominator to variant frequency calculations) may be over-reported in certain regions because of multi-amplicon reporting in RACE-like data.

21. -a, --primer\_nm_allowance

If -f/--primer_fasta, allow up to this number of edit operations for matching primers in the start of R2. Default 3. The last sixteen 3' nucleotides of the matched primer will be appended to the read names to avoid UMI grouping and consensus deduplication. R2s that do not match any primer will be reassigned the unknown primer regex X{16}  


## Accessory scripts

To facilitate simulation of reads and variants in a desired transcript *de novo*, accessory scripts are provided in the scripts directory. Code here is not fully tested and is only provided for convenience to facilitate error modeling.

1. run\_read_generator.py.
This may be used to simulate paired-end reads with optional random addition of noise. However, I recommend one of the many NGS read simulators that construct error models from real data to generate test reads.

To simulate error-free, tiled RNA sequencing alignments in full-length coding regions, I provide a transcript reference FASTA and target BED file in transcriptomic space. I split the transcript target region into smaller targets (read length minus the sum of the average lengths of two primers, ~ 50 bp) then run with --make\_amplicons but *not* --rna. The use of --rna is limited to simulating RNA reads starting from a genomic reference as opposed to a transcript reference. Mock primers can be generated that are flush with the termini of the segmented target regions.

As an example, for 2 x 150 bp reads I create ~100 bp interleaved target regions in BED format and configure the number of reads to generate for each amplicon by the BED score field. Finally, I pass --make\_amplicons so that reads start from the terminus of the target segments.

--make\_amplicons is intended for direct simulation of reads from a transcript FASTA and target features. For general DNA or RNA read generation using a genome FASTA and standard GFF annotations, omit --make\_amplicons. In this mode, generated fragments start and end at random coordinates in the sequence space informed by the --frag\_length and --read\_length arguments.

2. run\_variant_generator.py
This script may be used to generate a VCF of all SNP and MNP codon permutations in a desired transcript coding region that match a mutagenesis signature.

3. run\_vcf\_subsampler.py
This script can be used to subsample variants from the VCF produced by run\_variant_generator.py. Balancing true and false positive variants is highly recommended for training models.

4. run\_ec\_data\_generator.py
As input, this script requires satmut\_utils call output summary.txt files for a true mutagenized library and a non-mutagenized negative control library. It then generates variants and configures variant frequencies by estimating parameters from the true mutagenized summary.txt file. It finally invokes satmut_utils sim and call to generate training data for modeling.

Secondary analysis functions are provided in R for parsing and summarizing of satmut\_utils call summary.txt files (summarization\_utils.R) and for training several models (modeling_utils.R).


## Tests

To run unit tests, execute the following from the satmut_utils repository:

```nose2 -v```
