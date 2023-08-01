# Methodology

See [here](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq3.html) for more information on the library generation of the sci-RNA-Seq3 protocol.

## Major steps

1. Check for sanity of provided barcodes and sample-sheet.
2. Converts BCL files to paired-end .fq.gz files with PCR indexes in header (**bcl2fastq**).
3. Splits paired-end .fq.fz files into smaller (evenly-sized) chunks for parallelization (**fastqsplitter**).
4. Demultiplexing using the supplied sample-specific barcodes (**sci-rocket**).
      * Finds exact or nearest match for PCR Index #1 (p5), PCR Index #1 (p7), ligation and/or RT barcode (single match with ≤1 hamming distance).
      * Generates sample-specific .fastq.gz files with corrected R1 sequence (48nt) and added read-names in R2.
      * Read-pairs without all four matching barcodes are discarded into separate .fastq.gz files with logs detailing which barcode(s) are (non-)matching.
5. Performs adapter and low-quality base-trimming (**fastp**).
      * Read-pairs with a mate ≤10nt after trimming are discarded.
6. Aligns reads to the supplied reference genome and perform cell-barcode/UMI counting (**STARSolo**).
      * STAR index can be generated based on supplied genome sequences and annotations.
      * Per gene and cellular barcode, intronic, exonics and UTR-overlapping reads (UMI) are counted and multi-mapping reads are distributed using the `EM` method.
7. Generate demultiplexing/alignment overview. (**sci-dash**)
        * Generates a HTML report with demultiplexing and alignment statistics.

> Parallization is performed per sequencing run and split chunk.

## Sample demultiplexing

**Example of R1 sequence:**

```text
  @READNAME 1:N:0:CCGTATGATT+AGATGCAACT
  ACTTGATTGTGAGAGCTCCGTGAAAGGTTAGCAT

  First 9 or 10nt:  Ligation barcode
  Next 8nt:    UMI
  Next 6nt:    Primer
  Last 10nt:   RT Barcode

  Anatomy of R1:
  |ACTTGATTGT| |GAGAGCTC| |CGTGAA| |AGGTTAGCAT|
  |-LIGATION-| |---UMI--| |Primer| |----RT----|

  Corrected R1 sequence (48nt):
  |CCGTATGATT| |CCGTATGATT| |ACTTGATTGT| |AGGTTAGCAT| |GAGAGCTC|
  |----p7----| |----p5----| |-LIGATION-| |----RT----| |---UMI--|
```

For sample-demultiplexing, the following steps are performed:

1. Extracts p5, p7 PCR indexes from the read-name of R1 and ligation, RT and UMI barcodes from sequence of read 1 (R1).
2. If no match, corrects p5, p7, ligation and/or RT barcode to nearest match (with max. 1nt difference). If multiple close matches, discard read-pair.
    * For ligation barcodes of 9nt in length, an extra G is added to the ligation sequence as padding to ensure 48nt R1 sequence.
3. Add the barcodes to the read-name of read 2 (R2): `@READNAME|P5-<p5>-P7-<p7>|<ligation>|<rt>_<UMI>`
4. Generate sample-specific paired-end fq.gz files with corrected R1 sequence (48nt) and R2 sequence.

## Output

The major output files are the following:

1. **Sequence and sample-specific fastq file(s)**:
      * `{sequencing_name}/demux_reads/{sample_name}_R1.fastq.gz`
      * `{sequencing_name}/demux_reads/{sample_name}_R2.fastq.gz`
      * `{sequencing_name}/demux_reads/{sequencing_name}_R1_discarded.fastq.gz`
      * `{sequencing_name}/demux_reads/{sequencing_name}_R2_discarded.fastq.gz`
      * `{sequencing_name}/demux_reads/log_{sequencing_name}_discarded_reads.tsv.gz`
2. **Alignment files**:
      * `{sequencing_name}/alignment/{sample_name}_{species}_Aligned_sortedByCoord_markDup.bam/bai`
      * `{sequencing_name}/alignment/{sample_name}_{species}_Solo.out/`
3. **Demultiplexing/alignment overview**:
      * `{sequencing_name}/sci-dash/`
