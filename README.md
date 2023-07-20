# snakemake-scirocket

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Github issues](https://img.shields.io/github/issues/odomlab2/snakemake-sciseq)](https://img.shields.io/github/issues/odomlab2/snakemake-sciseq)

> **snakemake-scirocket** is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow which performs processing of sci-RNA-Seq3 data.

## Pipeline summary

```text
  Input:  Paired-end sequencing data from sci-RNA-Seq3 protocol (BCL files)
```

### **Demultiplexes sci-RNA-Seq3 samples (Illumina NovaSeq)**

---

- Check for sanity of provided barcodes and sample-sheet.
- Converts BCL files to R1/R2.fastq.gz files with p5 and p7 indexes in header (**bcl2fastq**).
- Splits R1/R2 fastq.gz files into smaller (evenly-sized) chunks for parallelization (**fastqsplitter**).
- Demultiplexes (paired-end) sequencing using the supplied sample-specific barcodes (**sci-rocket**).
  - Finds exact or nearest match for p5, p7, ligation and/or RT barcode (<=1 hamming distance with only single match).
  - Generates sample-specific .fastq.gz file(s) with correct read-name for R2.
  - Read-pairs with no p5, p7, ligation and/or RT barcode match are discarded into separate .fastq.gz files.
    - Log file contains information on discarded read-pairs detailing which barcodes are (non-)matching.

### **Processing of sci-RNA-Seq3 data**

---

- Performs adapter and low-quality base-trimming (**fastp**).
- Aligns reads to the supplied reference genome and perform cell-barcode/UMI counting (**STARSolo**).
- Marks duplicates (**sambamba**).
- Generates QC:
  - Alignment stats (**sambamba**)
  - Generate demultiplexing/alignment overview. (**sci-dash**)

> **Note**  
> See [Methodology](#methodology) for more information on various key aspects of the methodology.  
> Where possible, parallization is performed per sequencing run and over samples.

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
- **p5**: PCR (p5) index (e.g. A01:H01) used to identify the sample during demultiplexing.
- **p7**: PCR (p7) index (e.g. G01:G12) used to identify the sample during demultiplexing.
- **barcode_rt**: RT barcode (e.g. P01-A01) used to identify the sample during demultiplexing.
- **sample_name**: Name of the demultiplexed sample.
- **species**: Reference species (e.g. mouse or human).

> **Note**
>
> - **p5** and **p7** are used to denote the PCR indexes belonging to a particular sample. The indexes are translated to all relevant combinations.
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

```shell
cd workflow/
snakemake --profile <profile_name> --configfile <path_config>
```

> **Useful parameters**
>
> - `-n`: Perform dry-run (generate commands without executing).
> - `-p`: Print shell commands.
> - `--notemp`: Do not remove files flagged as temporary.
> - `--rerun-incomplete`: Rerun all jobs with missing output files.

## Output

The major output files are the following:

- **Sequence and sample-specific fastq file(s)**:
  - `{sequencing_name}/demux_reads/{sample_name}_R1.fastq.gz`
  - `{sequencing_name}/demux_reads/{sample_name}_R2.fastq.gz`
  - `{sequencing_name}/demux_reads/{sequencing_name}_R1_discarded.fastq.gz`
  - `{sequencing_name}/demux_reads/{sequencing_name}_R2_discarded.fastq.gz`
  - `{sequencing_name}/demux_reads/log_{sequencing_name}_discarded_reads.tsv.gz`
- **Alignment files**:
  - `{sequencing_name}/alignment/{sample_name}_{species}_Aligned_sortedByCoord_markDup.bam/bai`
  - `{sequencing_name}/alignment/{sample_name}_{species}_Solo.out/`
- **Demultiplexing/alignment overview**:
  - `{sequencing_name}/sci-dash/`

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

#### **LSF profile**

These commands were used to set up the LSF profile on the DKFZ LSF cluster:

```bash
# Set up snakemake profile
micromamba activate sci-rocket
module rm python
pip install --user cookiecutter

# create configuration directory that snakemake searches for profiles
profile_dir="/home/<username>/.config/snakemake"
mkdir -p "$profile_dir"
# use cookiecutter to create the profile in the config directory
template="gh:Snakemake-Profiles/lsf"
cookiecutter --output-dir "$profile_dir" "$template"

# parameters to set
LSF_UNIT_FOR_LIMITS=MB
UNKWN_behaviour=wait
ZOMBI_behaviour=ignore
use_conda=False
use-singularity=False
latency-wait=30
printshellcmds=True
restart-times=2
jobs=50
max-jobs-per-second=10
max-status-checks-per-second=10
profile=lsf_dkfz
```
