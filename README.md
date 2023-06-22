# snakemake-scirocket

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Github issues](https://img.shields.io/github/issues/odomlab2/snakemake-sciseq)](https://img.shields.io/github/issues/odomlab2/snakemake-sciseq)

> **snakemake-scirocket** is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow which performs processing of sci-RNA-Seq3 data.

## Pipeline summary

```text
  Input:  Paired-end sequencing data from sci-RNA-Seq3 protocol (BCL files)
```

### **Demultiplexes and processes sci-RNA-Seq3 data from Illumina NovaSeq sequencer**

---

- Check for sanity of provided barcodes and sample-sheet.
- Converts BCL files to R1/R2.fastq.gz files with p5 and p7 indexes in header (**bcl2fastq**).
- Splits R1/R2 fastq.gz files into smaller (evenly-sized) chunks for parallization (**fastqsplitter**).
- Demultiplexes (paired-end) sequencing using the supplied sample-specific barcodes.
  - Finds exact or nearest match for p5, p7, ligation and/or RT barcode (<=1 hamming distance with only single match).
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

- **path_bcl**: Path to BCL files.
- **sequencing_name**: Sequencing name (e.g., run123), used to store sequencing-specific files.
- **experiment_name**: Experiment name (e.g, experiment_x).
- **P5**: PCR (p5) index (e.g. A01:H01) used to identify the sample during demultiplexing.
- **P7**: PCR (p7) index (e.g. G01:G12) used to identify the sample during demultiplexing.
- **barcode_rt**: RT barcode (e.g. P01-A01) used to identify the sample during demultiplexing.
- **sample_name**: Name of the demultiplexed sample.
- **species**: Reference species (e.g. mouse or human).

> **Note**
>
> - **P5** and **P7** are used to denote the PCR indexes belonging to a particular sample. The indexes are translated to all relevant combinations.
>   - To specify one or multiple p5/p7 strips, use the following format:
>     - p5 (1 strips): `A01:H01`
>     - p5 (1.5 strips): `A01:H01,A02:D02`
>     - p5 (2 strips): `A01:H01,A02:H02`
>     - p7 (1 strip): `G01:G12`
>     - p7 (1.5 strips): `G01:G12,H01:H06`
>     - p7 (2 strips): `G01:G12,H01:H12`
> - **species** should be present in the `config.yaml` file with their respective genome sequences (.fa) and gene-annotations (.gtf) used to generate mapping indexes.

### Barcode design

The workflow requires a file (.tsv) containing the barcodes used in the experiment with at least the following required columns:

- **type**: Type of barcode (`ligation`, `p5`, `p7` or `rt`).
- **barcode**: Name of the barcode (e.g. A01).
- **sequence**: Nucleotide sequence of the barcode.

## Usage

To run this workflow on the DKFZ LSF cluster, first [set-up the proper LSF profile](https://github.com/Snakemake-Profiles/lsf) and run the following command for either the WGS or WTS workflow:

`snakemake --profile lsf_dkfz -n --cluster-config cluster.json --configfile <config>`

> **Note**  
> Remove `-n` to disable dry-run.
> `--cluster-config cluster.json` is optional, but recommended for better load-sharing on cluster systems.

## Output

The major output files are the following:

- **Sequence and sample-specific fastq file(s)**:
  - `{sequencing_name}/demux_reads/{sample_name}_R1.fastq.gz`
  - `{sequencing_name}/demux_reads/{sample_name}_R1.fastq.gz`
  - `{sequencing_name}/demux_reads/{sequencing_name}_R1_discarded.fastq.gz`
  - `{sequencing_name}/demux_reads/{sequencing_name}_R1_discarded.fastq.gz`
- **Sample-specific R2 file used for alignment (fastp-trimmed)**:
  - `{sequencing_name}/demux_reads/{sample_name}_fastp_R2.fastq.gz`

## Methodology

### Demultiplexing scheme

See [here](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq3.html) for more information on the library generation of the sci-RNA-Seq3 protocol.

#### **Sample demultiplexing**

For sample-demultiplexing, the following steps are performed:

1. Extracts p5, p7 indexes from the read-name of R1 and ligation and RT barcodes from sequence of read 1 (R1).

   - The following scheme is used for paired-end Illumina NovaSeq:

   ```text
     Example R1:
     @READNAME 1:N:0:<p7>+<p5>
     ACTTGATTGTGAGAGCTCCGTGAAAGGTTAGCAT

     First 9 or 10nt:  Ligation barcode
     Next 8nt:    UMI
     Next 6nt:    Primer
     Last 10nt:   RT Barcode

     Anatomy of R1:
     |ACTTGATTGT| |GAGAGCTC| |CGTGAA| |AGGTTAGCAT|
     |-LIGATION-| |---UMI--| |Primer| |----RT----|

   ```

2. If no match, corrects p5, p7, ligation and/or RT barcode to nearest match (with max. 1nt difference). If multiple close matches, discard read-pair.
3. Set the p5, p7, RT, ligation barcode and UMI as read-name in read 2 (R2): `@READNAME|P5-<p5>-P7-<p7>|<ligation>|<rt>_<UMI>`

> **Note**  
> Read-pairs with unmatched p5, p7, ligation or RT barcodes are discarded into separate R1/R2 fastq.gz files.
> Log files are generated which denote the mismatching barcodes per read (R1).
