# satmut_utils:

satmut_utils is a Python package for simulation and variant calling of saturation mutagenesis data. The two main subcommands are:
1. sim
2. call

satmut_utils commands are designed to simulate and call variants in paired-end, targeted RNA-seq datasets. That is, alignments to a single transcript are expected. 

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

5. Download reference files if using Ensembl identifiers.

Create a reference directory:
```
REF_DIR="./path_to_refs"
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

6. Navigate to the package parent directory within the repository:
```
cd satmut_utils
```

## Reference files

For convenience, a curated transcriptome of primary human transcripts from APPRIS is provided, allowing the user to pass an Ensembl gene or transcript ID to satmut_utils. However, if the requested Ensembl ID is not found in this set, the user must pass their own reference files. These include:

A. Transcript reference (FASTA)

B. Transcript annotations (GFF)

C. GFF reference (FASTA)

Common transcript annotations in GFF format map coordinates in the genome. For this typical case, file A should specify a transcript FASTA and file C should specify the genome FASTA.

In typical saturation mutagenesis datasets, an intron-less coding sequence, often lacking endogenous untranslated (UTR) regions, is expressed from a vector. With the three reference files above, the user can set A and C to the composite (vector + target CDS) reference FASTA, then make a custom GFF annotation (file B). See documentation for more information on creating reference files.


## Code examples

Parameter help:
```
python satmut_utils.py -h
python satmut_utils.py sim -h
python satmut_utils.py call -h
```

It is recommended that the user create a new output directory for each job.
```OUTPUT_DIR="./satmut_utils_test"```

Common arguments to both sim and call subcommands are provided first, then the subcommand, and then subcommand-specific arguments.


### Run sim

Run sim on *in silico* alignments to generate SNPs, MNPs, and InDels:

```
TEST_DIR="tests/test_data"

python satmut_utils.py -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed sim -f -a $TEST_DIR/CBS_sim.bam -v $TEST_DIR/CBS_sim.vcf
```

### sim outputs

The sim workflow outputs paired FASTQs and a truth VCF containing expected variants and their frequencies. Optionally, edited reads may be filtered and realigned for further investigation.


### Run call

Run call on our simulated data by specifying an Ensembl transcript or gene ID and the directory containing curated reference files:
```
TEST_DIR="tests/test_data"

python satmut_utils.py -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed call -1 $TEST_DIR/CBS_sim.R1.fq.gz -2 $TEST_DIR/CBS_sim.R2.fq.gz -v
```

Here, we call variants on a simulated dataset with no adapter sequences, so we pass -v. However, in typical cases, the user should provide 5' and 3' adapters for trimming.


If the Ensembl ID is not in the curated set of primary transcripts, or if the user wishes to align to a custom reference, several reference files most be provided (see Reference Files):

```
python satmut_utils.py -r $TEST_DIR/CBS.fa -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed call -1 $TEST_DIR/CBS_sim.R1.fq.gz -2 $TEST_DIR/CBS_sim.R2.fq.gz -v -g $TEST_DIR/CBS.gff -k $REF_DIR/GRCh38.fa
```

### call outputs

The call workflow produces a VCF of candidate variant calls as well as a BED file reporting fragment coverage across the reference.

The output VCF and its corresponding tab-delimited summary.txt file contain records for each mismatched base in an MNP, so that quality information for the mismatches can be used for learning-based error correction.

A number of useful R functions exist in prototype.summarization_utils.r for parsing and summarizing the VCF summary file. See the VCF header for column/feature descriptions.


## Tests

To run unit tests, execute the following from the satmut_utils repository:

```nose2 -q```


## Accessory scripts

To facilitate simulation of reads and variants in a desired transcript de novo, additional command-line scripts are provided in the scripts directory. Code here is not fully tested and is only provided for convenience.

To call these scripts, make sure you are in the satmut_utils directory and execute as modules:

```
python -m scripts.run_bowtie2_aligner -h
``` 

1. run_read_generator.py.
This may be used to simulate paired-end RNA reads with random addition of noise. However, we recommend one of the many NGS read simulators that construct error models to generate test reads.

2. run_variant_generator.py
This script may be used to generate a VCF of all SNP and MNP codon permutations in a desired transcript coding region.

3. run_vcf_subsampler.py
This script can be used to subsample variants from the VCF produced by run_variant_generator.py.

4. run_ec_data_generator.py
As input, this script requires satmut_utils call output summary.txt files for a true mutagenized library and a non-mutagenized negative control library. It then generates variants and configures variant frequencies by estimating parameters from the true mutagenized summary.txt file. It finally invokes satmut_utils sim and call to generate training data for modeling.

Secondary analysis functions are provided in R for parsing and summarizing of satmut_utils call summary.txt files (summarization_utils.R) and for training several models (modeling_utils.R).

