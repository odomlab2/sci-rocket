# snakemake-sciseq

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Github issues](https://img.shields.io/github/issues/odomlab2/snakemake-sciseq)](https://img.shields.io/github/issues/odomlab2/snakemake-sciseq)

> **snakemake-sciseq** is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow which performs processing of sci-RNA-Seq3 data.

## Pipeline summary

```text
  Input:  Paired-end sequencing data from sci-RNA-Seq3 protocol (BCL files)
```

### **Demultiplexes and processes sci-RNA-Seq3 data from Illumina NovaSeq sequencer**

---

- Check for sanity of provided barcodes and sample-sheet.
- Converts BCL files to R1/R2.fastq.gz files with p5 and p7 indexes in header (**bcl2fastq**).
- Splits R1/R2 fastq.gz files into smaller (evenly-sized) chunks for parallization (**fastqsplitter**).
- Demultiplexes (paired-end) sequencing using the supplied sample-specific RT and ligation barcodes.
  - Finds exact or nearest match for p5, p7, ligation and/or RT barcode (<=1 Levenshtein distance with only single match).
  - Generates sample-specific .fastq.gz file(s) with correct read-name for R2.
  - Read-pairs with no p5, p7, ligation and/or RT barcode match are discarded into separate .fastq.gz files.
    - Log file contains information on discarded read-pairs detailing which barcodes are matching.
- Generates QC:
  - Demultiplexing stats

### **Processing of sci-RNA-Seq3 data**

---

- Performs adapter and low-quality base-trimming (**fastp**).
- Aligns reads to the supplied reference genome (**STAR**).
- Marks duplicates (**sambamba**).
- Generates QC:
  - Alignment stats (**sambamba**)
  - Sequencing coverage (**mosdepth**)

> **Note**  
> See [Methodology](#methodology) for more information on various key aspects of the methodology.  
> Where possible, parallization is performed per sequencing run.

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

- **sequencing_name**: Sequencing sample (e.g. AS-123456), used to determine respective .fastq file(s) and from which sample-specific BAM file(s) are generated.
- **barcode_rt**: RT barcode (e.g. P01-A01) used to identify the sample during demultiplexing.
- **sample_name**: Name of the RT-barcoded sample.
- **species**: Reference species (e.g. mouse or human).

> **Note**
>
> - **species** should be present in the `config.yaml` file with their respective genome sequences (.fa) and gene-annotations (.gtf) used to generate mapping indexes.

### Barcode design

The workflow requires a file (.tsv) containing the barcodes used in the experiment with at least the following required columns:

- **type**: Type of barcode (`ligation`, `p5`, `p7` or `rt`).
- **barcode**: Name of the barcode (e.g. A01).
- **sequence**: Nucleotide sequence of the barcode.

## Usage

To run this workflow on the DKFZ LSF cluster, first [set-up the proper LSF profile](https://github.com/Snakemake-Profiles/lsf) and run the following command for either the WGS or WTS workflow:

`snakemake --profile lsf_dkfz -n`

> **Note**  
> Remove `-n` to disable dry-run.

## Output

The major output files are the following:

- **Sequence and sample-specific fastq file(s)**:
  - `fastq/{sequencing_name}/demux/{sample_name}_R1.fastq.gz`
  - `fastq/{sequencing_name}/demux/{sample_name}_R2.fastq.gz`
  - `fastq/{sequencing_name}/demux/{sequencing_name}_R1_discarded.fastq.gz`
  - `fastq/{sequencing_name}/demux/{sequencing_name}_R2_discarded.fastq.gz`
- **Sample-specific R2 file used for alignment (fastp-trimmed)**:
  - `fastp/{sequencing_name}/{sample_name}_R2_fastq.gz`
- **Sample-specific BAM file(s)**:
  - `alignment/{reference}/{sample}_sortedByCoord.bam`

> **Note**  
> {sample} is based on **sample** in the sample sheet.

## Methodology

### Demultiplexing scheme

See [here](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq3.html) for more information on the library generation of the sci-RNA-Seq3 protocol.

#### **Sample demultiplexing**

For sample-demultiplexing, the following steps are performed:

1. Extracts ligation and RT barcodes from read 1 (R1).

   - The following scheme is used for paired-end Illumina NovaSeq:

   ```text
     Example R1:
     ACTTGATTGTGAGAGCTCCGTGAAAGGTTAGCAT

     First 9 or 10nt:  Ligation barcode
     Next 8nt:    UMI
     Next 6nt:    Primer
     Last 10nt:   RT Barcode (sample-specific)

     Anatomy of R1:
     |ACTTGATTGT| |GAGAGCTC| |CGTGAA| |AGGTTAGCAT|
     |-LIGATION-| |---UMI--| |Primer| |----RT----|

   ```

2. If no match, correct RT barcode to nearest match (with max. 1nt difference). If multiple close matches, discard read-pair.
3. Set the RT, ligation barcode and UMI as read-name in read 2 (R2).

> **Note**  
> Read-pairs with unmatched RT barcodes are discarded into separate R1/R2 fastq.gz files.
