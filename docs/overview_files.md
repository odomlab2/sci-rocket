# Sample-sheet and barcodes

See [example sample-sheet](https://github.com/odomlab2/sci-rocket/blob/main/workflow/examples/example_samplesheet.tsv).

**sci-rocket** requires a sample sheet (.tsv) with at least the following required columns.

**Either one or both of the following columns:**

* **path_bcl**: Path to folder containing the BCL files.
* **path_fastq**: Path to folder containing the `Undetermined_S0_R1_001.fastq.gz` and `Undetermined_S0_R2_001.fastq.gz` files (if `bcl2fastq` has already been run).

**Additional required columns:**

* **experiment_name**: Experiment name (e.g., experimentXYZ), used to associate all downstream files and underlying samples.
* **p5**: PCR (p5) index(es) (e.g. A01:H01, or **column(s)** of a 96-well index plate) used to identify the sample during demultiplexing.
* **p7**: PCR (p7) index(es) (e.g. G01:G12, or **rows(s)** of a 96-well index plate) used to identify the sample during demultiplexing.
* **rt**: RT barcode(s) (e.g. P01-A01:P01-A12) used to identify the sample during demultiplexing.
* **sample_name**: Name of the demultiplexed sample, used to generate sample-specific files.
* **species**: Reference species (e.g. mouse or human).
* **n_expected_cells**: Number of expected cells in the (demultiplexed) sample (used during UMI filtering).

> * **p5** and **p7** are used to denote the PCR indexes belonging to a particular sample / cell. The indexes are translated to all relevant combinations within the sequencing-run using the 96-well index plate layout of 8x rows (A:H) and 12 columns (01:12). To specify one or multiple p5/p7 ranges, use the following format:
>   * p5 (1 column): `A01:H01`
>   * p5 (1.5 columns): `A01:H01,A02:D02`
>   * p5 (2 columns): `A01:H01,A02:H02`
>   * p7 (1 row): `G01:G12`
>   * p7 (1.5 rows): `G01:G12,H01:H06`
>   * p7 (2 rows): `G01:G12,H01:H12`
> * The **rt** is used to denote the RT barcode belonging to a particular sample / cell. The indexes are translated to all relevant combinations within the sequencing-run. To specify one or multiple RT strips, use the following format:
>   * One RT: `P01-A01`
>   * Multiple RT (1 row): `P01-A01:P01-A12`
>   * Multiple RT (1 column): `P01-A01:P01-H01`
>   * Multiple RT (2 columns): `P01-A01:P01-H02`
>   * Multiple RT (rectangular region): `P01-B02:P01-E04`; This will include all RT barcodes from the rectangle with P01-B02 at the top left and P01-E04 at the bottom right:
>
> | P01  | 1 | 2 | 3 | 4 | 5 |
> | ---  | - | - | - | - | - |
> | A    | . | . | . | . | . |
> | B    | . | X | X | X | . |
> | C    | . | X | X | X | . |
> | D    | . | X | X | X | . |
> | E    | . | X | X | X | . |
> | F    | . | . | . | . | . |
>
>   * Multiple RT (multiple plates): `P01-B02:P01-E04,P02-A01:P02-A12`; This will include the same RT barcodes from P01 as the previous example, plus row A from P02.
>
> * **species** should be present in the `config.yaml` file with their respective genome sequences (.fa) and gene-annotations (.gtf) used to generate mapping indexes.

## Barcode design

The workflow requires a file (.tsv) containing the barcodes used in the experiment with at least the following required columns:

* **type**: Type of barcode (`ligation`, `p5`, `p7` or `rt`).
* **barcode**: Name of the barcode (e.g. A01).
* **sequence**: Nucleotide sequence of the barcode.

## Hashing sheet

The hashing workflow requires a separate file (.tsv) containing the hashing schematics used in the sample with at least the following required columns:

* hash_name: Name of the hashing experiment (e.g. hash_exp1).
* barcode: Sequence of the respective hashing barcode (e.g. GGTTGGCGAC).

To specify which samples are to be hashed (and using which hashing-sheet), add an additional column (`hashing`) in the sample-sheet to each each sample with the respective path to the hashing sheet.

## Haplotyping (optional; _Mus musculus_ cross-experiments only)

As optional procedure, **sci-rocket** can be used to further haplotype the sex-chromosome X of the demultiplexed samples, e.g. in the case of mouse F1 cross-hybrids. For this, the following columns can be added to the sample-sheet:

* **strain1**: Name of the first strain (e.g. B6 (C57BL/6J)).
* **strain2**: Name of the second strain (e.g. CAST/EiJ).

> These strains should be present in the [Mouse Genome Project (MGP) database](https://www.sanger.ac.uk/science/data/mouse-genomes-project).
> For C57BL/6J (wt), use `B6` as strain name.

This will add haplotype-specific read tags (HP) to the STARSolo BAM files and will output an additional `haplotyping` folder which contains cell-based read-counts per gene per haplotype (H1, H2, UA) for chromosome X.
