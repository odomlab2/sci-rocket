# Introduction

**sci-rocket** is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow which performs processing of sci-RNA-Seq3 data, including sample demultiplexing and downstream alignment and UMI-counting using [STARSolo](https://github.com/alexdobin/STAR).

Please see the set-up instructions below for more information on how to install and run the workflow.

## Pre-requirements

1. A conda system, e.g., [conda](https://docs.conda.io/en/latest/), [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [micromamba](https://micromamba.readthedocs.io/en/latest/)
2. Cluster-specific Snakemake configuration for batch-job submission
      * E.g., [LSF](https://github.com/Snakemake-Profiles/lsf) or [SLURM](https://github.com/Snakemake-Profiles/slurm)

We make use of a pre-defined [conda](https://docs.conda.io/en/latest/) environment in which all software dependencies (incl. Snakemake) are installed (`workflow/envs/env.yaml`).

## Set-up

1. Clone the repository:

      ```bash
      git clone https://github.com/odomlab2/sci-rocket
      ```

2. Create the conda environment (default name: `sci-rocket`):

      ```bash
      cd sci-rocket
      micromamba create -f workflow/envs/env.yaml
      ```

3. Activate the conda environment:

      ```bash
      micromamba activate sci-rocket
      ```

4. Run the workflow:

      ```bash
      cd workflow/
      snakemake --profile <profile_name> --configfile <path_config>
      ```

**Useful Snakemake parameters**:

> * `-n`: Perform dry-run (generate commands without executing).
> * `-p`: Print shell commands.
> * `--notemp`: Do not remove files flagged as temporary.
> * `--rerun-incomplete`: Rerun all jobs with missing output files.

## Configuration

The workflow requires a configuration file (`config.yaml`) which can be copied from the [example configuration file](https://github.com/odomlab2/sci-rocket/blob/main/workflow/examples/example_config.yaml) and adjusted to your needs.

Within the configuration file, the [sample-sheet](overview_files.md) (`path_samples`) needs to be specified. This file contains the sample names and paths to the raw sequencing data.
