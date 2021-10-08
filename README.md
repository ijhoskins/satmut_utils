# satmut_utils:

satmut_utils is a Python package for simulating and calling SNPs and MNPs in saturation mutagenesis datasets. The two main subcommands are:
1. sim
2. call

satmut_utils is designed to call variants in paired-end, targeted RNA-seq datasets. That is, alignments to a single transcript are expected. For germline and somatic variant calling in whole genomes or transcriptomes, other tools are recommended.


## Setting up the running environment

To get started, you must clone the satmut_utils repository and create a conda environment with the required dependencies. This process is outlined as follows:

1. Clone the repository:

Change to directory of your choice and clone the repository:

```git clone https://github.com/ijhoskins/satmut_utils.git```


2. Download curated reference files if using Ensembl identifiers.
TBD
```REF_DIR="/path_to_refs"```


3. Install miniconda for managing environments and packages.

For convenience, a miniconda installer script is provided for Linux distributions under the environment directory. Otherwise, see the following link for installation instructures for your particular hardware architecture.
https://docs.conda.io/en/latest/miniconda.html

4. Create the conda environment:
```cd satmut_utils && conda env create -f satmut_utils_env.yaml```

5. Activate the environment
```conda activate satmut_utils```

6. Navigate to the calling directory within the repository:
```cd satmut_utils; python 


## Reference files

For convenience, a curated transcriptome of primary transcripts from APPRIS is provided, allowing the user to pass an Ensembl gene or transcript ID to satmut_utils. However, if the requested Ensembl ID is not found in this set, the user must pass their own reference files. These include:

A. Transcript reference (FASTA)

B. Transcript annotations (GFF)

C. GFF reference (FASTA)

Common transcript annotations in GFF format map coordinates in the genome. For this typical case, file A should specify a transcript FASTA and file C should specify the genome FASTA.

In typical saturation mutagenesis datasets, an intron-less coding sequence, often lacking endogenous untranslated (UTR) regions, is expressed from a vector. With the three reference files noted above, the user can set A and C to the composite (vector and target CDS) reference FASTA, then make a custom GFF annotation (file B). See documentation for more information on creating reference files.


## Code examples

Parameter help:
```
python satmut_utils -h
python satmut_utils sim -h
python satmut_utils call -h
```

It is recommended that the user create a new running directory/folder for each job.
```OUTPUT_DIR="/path_to_output_folder"```

Common arguments to both sim and call subcommands are provided first, then the subcommand, and then subcommand-specific arguments.


Run the sim workflow:

IH TODO


Run the call workflow by specifying an Ensembl transcript ID and the directory containing curated reference files:
```
python satmut_utils.py -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5  TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT
```

The Ensembl ID may also specify a gene:
```
python satmut_utils.py -i ENSG00000160200.17 -x $REF_DIR -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5  TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT
```

If an Ensembl gene ID has more than one transcript isoform, satmut_utils will select the first transcript in the gencode.v29.annotation.gtf. IH TODO: select the only transcript available in the transcriptome GTF.

If the Ensembl ID is not in the curated set of primary transcripts, the user must provide their own reference files:
```
python satmut_utils.py -r CBS.fa -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5  TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT -g CBS.gtf -k GRCh38.fa
```

Additional files may be passed, such as a primer and target BED file:
```
python satmut_utils.py -i ENST00000398165.7 -x $REF_DIR -p primers.bed -o $OUTPUT_DIR call -1 R1.fq.gz -2 R2.fq.gz -5 TACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCT -t target.bed
```

More than one 5' adapter and more than one 3' adapter are often needed to additionally trim vector sequences (e.g. attB sites) from reads of terminal PCR tiles that span the vector-CDS junctions. In this case, provide multiple comma-delimited adapters:
```-5 TACACGACGCTCTTCCGATCT,CAAGTTTGTACAAAAAAGTTGGC -3 AGATCGGAAGAGCACACGTCT,CCAACTTTCTTGTACAAAGTGGT
```

## call outputs

The call workflow produces a VCF of candidate variant calls as well as a BED file reporting fragment coverage across the reference.

The output VCF and its corresponding tab-delimited summary.txt file contain records for each mismatched base in an MNP, so that quality information for the mismatches can be used for learning-based error correction.

A number of useful R functions exist in prototype.summarization_utils.r for parsing and summarizing the VCF summary file. See the VCF header for column/feature descriptions.


## Tests

To run unit tests, execute the following from the satmut_utils repository:

```nose2 -v```


## Accessory scripts

To facilitate simulation of reads and variants in a desired transcript de novo, additional command-line scripts are provided in the scripts directory. Code here is not fully tested and is only provided for convenience.

1. run_read_generator.py.
This may be used to simulate paired-end RNA reads with random addition of noise. However, we recommend one of the many NGS read simulators that construct error models to generate test reads.

2. run_variant_generator.py
This script may be used to generate a VCF of all SNP and MNP codon permutations in a desired transcript coding region.

3. run_vcf_subsampler.py
This script can be used to subsample variants from the VCF produced by run_variant_generator.py.

4. run_ec_data_generator.py
As input, this script requires satmut_utils call output summary.txt files for a true mutagenized library and a non-mutagenized negative control library. It then generates variants and configures variant frequencies by estimating parameters from the true mutagenized summary.txt file. It finally invokes satmut_utils sim and call to generate training data for modeling.

Secondary analysis functions are provided in R for parsing and summarizing of satmut_utils call summary.txt files (summarization_utils.R) and for training several models (modeling_utils.R).

