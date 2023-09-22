# Sample-sheet and barcodes

**sci-rocket** requires a sample sheet (.tsv) with at least the following required columns:

* **path_bcl**: Path to folder containing the BCL files.
* **sequencing_name**: Sequencing name (e.g., run123), used to store sequencing-specific files.
* **p5**: PCR (p5) index (e.g. A01:H01) used to identify the sample during demultiplexing.
* **p7**: PCR (p7) index (e.g. G01:G12) used to identify the sample during demultiplexing.
* **barcode_rt**: RT barcode (e.g. P01-A01) used to identify the sample during demultiplexing.
* **sample_name**: Name of the demultiplexed sample, used to generate sample-specific files.
* **species**: Reference species (e.g. mouse or human).
* **n_expected_cells**: Number of expected cells in the sample (used during UMI filtering).

See [example sample-sheet](https://github.com/odomlab2/sci-rocket/blob/main/workflow/examples/example_samplesheet.tsv).

> * **p5** and **p7** are used to denote the PCR indexes belonging to a particular sample. The indexes are translated to all relevant combinations within the sequencing-run.
>   * To specify one or multiple p5/p7 strips, use the following format:
>     * p5 (1 strips): `A01:H01`
>     * p5 (1.5 strips): `A01:H01,A02:D02`
>     * p5 (2 strips): `A01:H01,A02:H02`
>     * p7 (1 strip): `G01:G12`
>     * p7 (1.5 strips): `G01:G12,H01:H06`
>     * p7 (2 strips): `G01:G12,H01:H12`
> * **species** should be present in the `config.yaml` file with their respective genome sequences (.fa) and gene-annotations (.gtf) used to generate mapping indexes.

## Barcode design

The workflow requires a file (.tsv) containing the barcodes used in the experiment with at least the following required columns:

* **type**: Type of barcode (`ligation`, `p5`, `p7` or `rt`).
* **barcode**: Name of the barcode (e.g. A01).
* **sequence**: Nucleotide sequence of the barcode.
