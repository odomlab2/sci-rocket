# Methodology

See [here](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq3.html) for more information on the library generation of the sci-RNA-Seq3 protocol.

## Major steps

1. Check for sanity of provided barcodes and sample-sheet.
2. Converts BCL files to paired-end .fq.gz files with PCR indexes in header (**bcl2fastq**).
   - Merges multiple sequencing runs (`path_bcl`) into one experiment-based file (`experiment_name`).
3. Splits paired-end .fq.fz files into smaller (evenly-sized) chunks for parallelization (**fastqsplitter**).
4. Demultiplexing using the supplied sample-specific barcodes (**sci-rocket**).
   - Finds exact or nearest match for PCR Index #1 (p5), PCR Index #1 (p7), ligation and/or RT barcode (single match with ≤1 hamming distance).
   - Generates sample-specific .fastq.gz files with corrected R1 sequence (48nt) and added read-names in R2.
   - Read-pairs without all four matching barcodes are discarded into separate .fastq.gz files with logs detailing which barcode(s) are (non-)matching.
   - For samples with a specified hashing sheet, additional hashing procedures are performed.
5. Performs adapter and low-quality base-trimming (**fastp**).
   - Read-pairs with a mate ≤10nt after trimming are discarded.
6. Aligns reads to the supplied reference genome and perform cell-barcode/UMI counting (**STARSolo**).
   - STAR index can be generated based on supplied genome sequences and annotations.
   - Per gene and cellular barcode, intronic, exonics and UTR-overlapping reads (UMI) are counted and multi-mapping reads are distributed using the `EM` method.
7. Generate demultiplexing/alignment overview. (**sci-dash**)
   - Generates a HTML report with demultiplexing and alignment statistics.

> Parallization is performed per experiment_name and split chunk.

### Optional steps

1. (_Mus musculus_-only) Haplotype demultiplexing.
   - Adds haplotype-specific read tags (HP) to the STARSolo BAM files using known haplotype-specific SNPs (**MGP** + **haplotag**).
   - Generate haplotype-specific read-counts per gene per cell (H1, H2, UA) (**umi_tools**).

## Downstream analysis

For downstream analysis, we also maintain an R package to analyze results produced by **sci-rocket** called [**scir**](https://github.com/odomlab2/scir).

## Sample demultiplexing (without hashing)

**Example of R1 sequence:**

```text
      @READNAME 1:N:0:CCGTATGATT+AGATGCAACT
                        |----p7---|+|----p5----|: p5 is reverse-complemented during demuxxing.
      ACTTGATTGTCAGAGCTTTGGTATCCTACCAGTT

      The R1 sequence should adhere to the following scheme:
      First 9 or 10nt:  Ligation barcode
      Next 6nt:    Primer
      Next 8nt:    UMI
      Last 10nt:   RT Barcode (sample-specific)

      Anatomy of R1 (ligation of 10nt):
      |ACTTGATTGT| |CAGAGC| |TTTGGTAT| |CCTACCAGTT|
      |-LIGATION-| |Primer| |---UMI--| |----RT----|

      Anatomy of R1 (ligation of 9nt):
      |CTCGTTGAT| |CAGAGC| |TTTGGTAT| |CCTACCAGTT| |T|
      |-LIGATION| |Primer| |---UMI--| |----RT----| |.| <- Extra base.

      Corrected R1 sequence (48nt):
      |CCGTATGATT| |AGTTGCATCT| |CTCGTTGAT| |CCTACCAGTT| |TTTGGTAT|
      |----p7----| |----p5----| |-LIGATION-| |----RT----| |---UMI--|
```

For sample-demultiplexing, the following steps are performed:

1. Extracts p5, p7 PCR indexes from the read-name of R1 and ligation, RT and UMI barcodes from sequence of read 1 (R1).
2. If no match, corrects p5, p7, ligation and/or RT barcode to nearest match (with max. 1nt difference). If multiple close matches, discard read-pair.
   - For ligation barcodes of 9nt in length, an extra G is added to the ligation sequence as padding to ensure 48nt R1 sequence.
3. Add the barcodes to the read-names of read 1 and 2:  
   `@READNAME|P5-<p5>-P7-<p7>|<ligation>|<rt>_<UMI>`
4. Generate sample-specific paired-end fq.gz files with corrected R1 sequence (48nt) and R2 sequence.

## Hashing

Reads (R2) containing both a polyA signal (AAAA) and a hashing barcode are used to flag reads as hashing-reads. These reads are used for collecting hashing metrics (with their respective R1) and subsequently removed from the analysis.

To flag reads as hashing-reads, we first check for the presence of the polyA signal (AAAA) in R2 (first occurence). If this signal is present, we check for the presence of the hashing barcode in R2 prior to this poly-A signal. It is assumed that the hashing barcodes are 10nt and are (directly) prior to the poly-A signal (5' - 1nt spacer).

If no match is found using the first 10nt (5' poly-A - 1nt spacer) ; we try again against the closest match (hamming distance=1). If no rescued match is found, we search for the presence of any hashing barcode in the entire R2 sequence prior to the poly-A signal.

The following metrics are generated from hashing reads, per cellular barcode / hash barcode combination:

```text
sequencing_name   hash_barcode    cell_barcode    count    n_umi
test            AGGTAGAGCT      F07_D09_LIG98_P01-C08      100      10
test            ACGTTGAATG      F07_D09_LIG98_P01-C08      200      15
```

These metrics are used to determine the hashing efficiency and to correct for UMI bias in downstream analysis:

- count: Total number of hashing reads for that specific cell-barcode / hash-barcode combination.
- n_umi: Number of unique UMIs for that specific cell-barcode / hash-barcode combination.

## Haplotyping (optional; _Mus musculus_ cross-experiments only)

As optional procedure, **sci-rocket** can be used to further haplotype the sex-chromosome X of the demultiplexed samples, e.g. in the case of mouse F1 cross-hybrids, see [here](overview_files.md#haplotyping-optional-mus-musculus-cross-experiments-only) for more information.
This will download (or symlink) the [MGP](http://www.sanger.ac.uk/science/data/mouse-genomes-project) database and perform haplotype-specific read-counting using **whatshap** on F1-informative heterozygous SNPs.

## Output

The major output files are the following:

1. **Sequence and sample-specific fastq file(s)**:
   - `{sequencing_name}/demux_reads/{sample_name}_R1.fastq.gz`
   - `{sequencing_name}/demux_reads/{sample_name}_R2.fastq.gz`
   - `{sequencing_name}/demux_reads/{sequencing_name}_R1_discarded.fastq.gz`
   - `{sequencing_name}/demux_reads/{sequencing_name}_R2_discarded.fastq.gz`
   - `{sequencing_name}/demux_reads/log_{sequencing_name}_discarded_reads.tsv.gz`
2. **Alignment files**:
   - `{sequencing_name}/alignment/{sample_name}_{species}_Aligned.sortedByCoord.out.bam/bai`
   - `{sequencing_name}/alignment/{sample_name}_{species}_Solo.out/`
3. **Demultiplexing/alignment overview**:
   - `{sequencing_name}/sci-dash/`
