# Introduction

**sci-rocket** is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow which performs processing of sci-RNA-seq3 sequencing, including barcode demultiplexing and downstream alignment / UMI-counting using [STARSolo](https://github.com/alexdobin/STAR).

Please see the set-up instructions below for more information on how to install and run the workflow.

## Pre-requirements

1. A conda system, e.g., [conda](https://docs.conda.io/en/latest/), [mamba](https://mamba.readthedocs.io/en/latest/) or [micromamba](https://micromamba.readthedocs.io/en/latest/)
2. Snakemake (v7.32.4) and a cluster-specific Snakemake configuration for batch-job submission (see instructions below).
      * E.g., [LSF](https://github.com/Snakemake-Profiles/lsf) or [SLURM](https://github.com/Snakemake-Profiles/slurm)

We make use of pre-defined environment(s) which houses all software dependencies (`workflow/envs/`). These are installed automatically by Snakemake when running the workflow (`--use-conda`).

## Set-up

1. Clone the repository:

      ```bash
      git clone https://github.com/odomlab2/sci-rocket
      ```

2. Download and install snakemake (e.g. using conda or micromamba):

      ```bash
      # This will install snakemake (7.32.4) into a new conda environment called 'snakemake'
      micromamba create -c conda-forge -c bioconda -n snakemake snakemake==7.32.4 mamba
      # Switch to the 'snakemake' environment
      micromamba activate snakemake
      ```

3. Run the workflow:

      ```bash
      cd workflow/
      snakemake --use-conda --profile <profile_name> --configfile <path_config>
      ```

**Useful Snakemake parameters**:

> * `-n`: Perform dry-run (generate commands without executing).
> * `-p`: Print shell commands.
> * `--notemp`: Do not remove files flagged as temporary.
> * `--rerun-incomplete`: Rerun all jobs with missing output files.

## Configuration

The workflow requires a configuration file (`config.yaml`) which can be copied from the [example configuration file](https://github.com/odomlab2/sci-rocket/blob/main/workflow/examples/example_config.yaml) and adjusted to your needs.

Within the configuration file, the [sample-sheet](overview_files.md) (`path_samples`) needs to be specified. This file contains the sample names and paths to the raw sequencing data (BCL or FASTQ).
