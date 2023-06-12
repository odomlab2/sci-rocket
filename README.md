# Introduction

**snakemake-sciseq** is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow which performs processing of sci-RNA-Seq3 data.

## Pipeline summary

### **Pre-processing**

---

1. Imports and sanity-checks the supplied sample sheet.
2. Demultiplexes (paired-end) sequencing (.fastq) using the supplied sample-specific RT and ligation design.
    - Generates sample-specific fastq file(s) with correct read-name for R2.
3. Performs adapter and low-quality base-trimming (**fastp**).
4. Performs quality control (**fastqc**).
5. Aligns reads to the supplied reference genome (**STAR**).
6. Marks duplicates (**sambamba**).
7. Generates QC:
   - Alignment stats (**sambamba**)
   - Insert size distribution (**picard**)
   - Sequencing coverage (**mosdepth**)

> **Note**  
> See [Methodology](#methodology) for more information on various key aspects of the methodology.  
> Where possible, parallization is performed per sample.

## Installation and pre-requirements

### Conda environment

This workflow makes use of a conda environment to handle dependencies:

```bash
# Creates conda environment (default name: snakemake-sciseq)
micromamba create -f environment.yml

# Active environment.
micromamba activate snakemake-sciseq
```

### Configuration

The workflow is configured using a `config.yaml` file. See `workflow/config.yaml.example` for an example configuration file.

### Sample sheet

The workflow requires a sample sheet (.tsv) with at least the following required columns:

- **sample**: Name of the sequencing sample (e.g. AS-123456)
- **sample_name**: Descriptive name of the sample

> **Note**  
> **sample** is used to retrieve multiple lane-specific runs of the same sample which are merged together after alignment and duplicate marking.

## Usage

To run this workflow on the DKFZ LSF cluster, first [set-up the proper LSF profile](https://github.com/Snakemake-Profiles/lsf) and run the following command for either the WGS or WTS workflow:

`snakemake --profile lsf_dkfz -n`

> **Note**  
> Remove `-n` to disable dry-run.

## Output

The major output files are the following:

- **Sample-specific fastq file(s)**:
  - `demultiplex_fastq/untrimmed/{sample}_R1.fq.gz`
  - `demultiplex_fastq/untrimmed/{sample}_R2.fq.gz`
  - `demultiplex_fastq/fastp/{sample}_R1.fq.gz`
  - `demultiplex_fastq/fastp/{sample}_R2.fq.gz`
- **Sample-specific BAM file(s)**:
  - `alignment/{reference}/{sample}_sortedByCoord.bam`

> **Note**  
> {sample} is based on the sample name in the sample sheet.

## Methodology

### Demultiplexing scheme

For demultiplexing, the following steps are performed:

1. Extracts RT and ligation barcode from read 1 (R1).
2. (Optional): If no match, correct RT and ligation barcode to nearest match (with edit distance <= 1).
3. Set the RT, ligation barcode and UMI as read-name in read 2 (R2).

> **Note**  
> Reads with unmatched RT or ligation barcodes are discarded.
