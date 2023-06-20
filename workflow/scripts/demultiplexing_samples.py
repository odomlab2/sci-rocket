import argparse
import sys
import os
import logging
from rich.console import Console
from rich.logging import RichHandler

import pandas as pd
import pysam
import gzip
from Levenshtein import distance as hamming
import pickle

from sanity_checks import retrieve_p5_barcodes, retrieve_p7_barcodes


def sciseq_sample_demultiplexing(log: logging.Logger, sequencing_name: str, samples: pd.DataFrame, barcodes: pd.DataFrame, path_r1: str, path_r2: str, path_out: str):
    """
    Performs demultiplexing of the raw fastq files based on the RT and ligation barcodes to produce sample-specific R1 and R2 files.
    The RT barcode is used to identify the sample and generate the sample-specific R1 and R2 files.
    The read names for R2 are modified to include the sample name, the RT and ligation barcodes.
    The ligation barcode can be either 9nt or 10nt long and this can affect the location of the UMI and RT barcodes.

    Read-pairs are discarded if:
        - The R1 read is empty.
        - The R1 read is shorter than 34nt.
        - The R1 read is longer than 34nt.
        - The RT barcode is not found in the sample sheet (even with 1nt correction).

    Args:
        log (logging.Logger): Logger.
        sequencing_name (str): Sequencing sample.
        samples (pd.DataFrame): Sample sheet of the samples in the sequencing run.
        barcodes (pd.DataFrame): Barcode sheet.
        path_r1 (str): Path to R1 fastq file.
        path_r2 (str): Path to R2 fastq file.
        path_out (str): Path to output directory.

    Returns:
        None
    """

    log.info("Starting sample-based demultiplexing of %s:\n(R1) %s\n(R2) %s", sequencing_name, path_r1, path_r2)

    # region Open file handlers --------------------------------------------------------------------------------------------------------------------------------

    # Open file handlers for input R1 and R2 files.
    try:
        fh_r1 = pysam.FastxFile(path_r1)
        fh_r2 = pysam.FastxFile(path_r2)
    except OSError:
        log.error("Could not find the R1 and R2 .fastq files, please check the paths:\n(R1) %s\n(R2) %s", path_r1, path_r2)
        sys.exit(1)

    # Open file handlers for discarded R1 and R2 files.
    path_r1_discarded = os.path.join(path_out, sequencing_name + "_R1_discarded.fastq.gz")
    path_r2_discarded = os.path.join(path_out, sequencing_name + "_R2_discarded.fastq.gz")

    try:
        fh_discarded_r1 = gzip.open(path_r1_discarded, "wt")
        fh_discarded_r2 = gzip.open(path_r2_discarded, "wt")
    except OSError:
        log.error(
            "Could not generate the discarded output files, please check the paths:\n(R1 discarded) %s\n(R2 discarded) %s",
            path_r1_discarded,
            path_r2_discarded,
        )
        sys.exit(1)

    # Get the unique RT barcode and samples contained in the sequencing run.
    # Generate a dictionary with the RT barcodes as keys and the sample names as values.
    dict_rt_barcodes = dict(zip(samples["barcode_rt"], samples["sample_name"]))

    # Open file handlers for all sample-specific output files.
    dict_fh_out = {k: {} for k in set(samples.sample_name)}

    for sample in dict_fh_out:
        # Generate the output paths.
        path_r1_out = os.path.join(path_out, sample + "_R1.fastq.gz")
        path_r2_out = os.path.join(path_out, sample + "_R2.fastq.gz")

        try:
            # Open file handlers for output R1 and R2 files.
            dict_fh_out[sample]["R1"] = gzip.open(path_r1_out, "wt")
            dict_fh_out[sample]["R2"] = gzip.open(path_r2_out, "wt")

        except OSError:
            log.error("Could not generate the sample-specific demultiplexed output files, please check the paths:\n(R1) %s\n(R2) %s", path_r1_out, path_r2_out)
            sys.exit(1)

    # Open a second file logger to write the info for the discarded reads.
    path_log_discarded = os.path.join(path_out, sequencing_name + "_discarded_reads.log")
    log_discarded = logging.getLogger("log_discarded")
    log_discarded.setLevel(logging.INFO)
    fh_log_discarded = logging.FileHandler(path_log_discarded)
    log_discarded.addHandler(fh_log_discarded)
    log_discarded.propagate = False

    # Header of discard log.
    log_discarded.info("read_name\tp5\tp7\tligation\trt\tumi")

    # endregion --------------------------------------------------------------------------------------------------------------------------------

    # region Generate barcode lookup tables --------------------------------------------------------------------------------------------------------------------------------

    # Retrieve the ligation barcodes.
    barcodes_ligation = barcodes.query("type == 'ligation'")

    # Get the possible list of unique RT barcodes contained in the sequencing run.
    rt_barcodes_sequencing = samples.query("sequencing_name == @sequencing_name")["barcode_rt"].unique()
    barcodes_rt = barcodes.query("type == 'rt' & barcode in @rt_barcodes_sequencing")

    # Retrieve the p5 and p7 barcodes used in this sequencing run.
    indexes_p5 = retrieve_p5_barcodes(log, samples["p5"].unique(), barcodes.query("type == 'p5'")["barcode"].unique())
    barcodes_p5 = barcodes.query("type == 'p5' & barcode in @indexes_p5")

    indexes_p7 = retrieve_p7_barcodes(log, samples["p7"].unique(), barcodes.query("type == 'p7'")["barcode"].unique())
    barcodes_p7 = barcodes.query("type == 'p7' & barcode in @indexes_p7")

    # Generate dicts for fast lookup.
    dict_rt = dict(zip(barcodes_rt["sequence"], barcodes_rt["barcode"]))
    dict_p5 = dict(zip(barcodes_p5["sequence"], barcodes_p5["barcode"]))
    dict_p7 = dict(zip(barcodes_p7["sequence"], barcodes_p7["barcode"]))
    dict_ligation = dict(zip(barcodes_ligation["sequence"], barcodes_ligation["barcode"]))

    # The p5 index needs to be reverse complemented.
    dict_p5 = {k[::-1].translate(str.maketrans("ATCG", "TAGC")): v for k, v in dict_p5.items()}

    # endregion --------------------------------------------------------------------------------------------------------------------------------

    def find_closest_match(sequence, comparison):
        # Calculate the hamming distance between sequence and all keys in dict.
        distances = {k: hamming(sequence, k, score_cutoff=1) for k in comparison.keys()}
        distances = {k: v for k, v in distances.items() if v == 1}

        # If there is only one key with a distance of 1, return the name and sequence.
        if len(distances) == 1:
            sequence = next(iter(distances))
            return comparison[sequence]
        else:
            return None

    # QC metrics.
    qc = {}
    qc["sequencing_name"] = sequencing_name
    qc["n_pairs"] = 0  # Total number of initial read-pairs.
    qc["n_pairs_success"] = 0  # Total number of read-pairs with correct RT, p5, p7 and ligation barcodes.
    qc["n_pairs_failure"] = 0  # Total number of discarded read-pairs due to various reason.

    qc["n_corrected_p5"] = 0  # Total number of read-pairs with 1bp mismatch in p5.
    qc["n_corrected_p7"] = 0  # Total number of read-pairs with 1bp mismatch in p7.
    qc["n_corrected_ligation"] = 0  # Total number of read-pairs with 1bp mismatch in ligation.
    qc["n_corrected_rt"] = 0  # Total number of read-pairs with 1bp mismatch in RT.

    qc["n_uncorrectable_p5"] = 0  # Total number of read-pairs with >1bp mismatch in p5.
    qc["n_uncorrectable_p7"] = 0  # Total number of read-pairs with >1bp mismatch in p7.
    qc["n_uncorrectable_ligation"] = 0  # Total number of read-pairs with >1bp mismatch in ligation.
    qc["n_uncorrectable_rt"] = 0  # Total number of read-pairs with >1bp mismatch in RT.

    # Dictionary to store the number of (succesfull) read-pairs, no. and type of identified RT barcodes, total UMI and total unique UMI per sample.
    samples_dict = {k: {"n_pairs_success": 0, "rt": {}} for k in samples.sample_name}

    # Iterate over the read-pairs and search for the indexes within R1.
    # If not any index is not found, try to correct that index with 1bp mismatch and search again in respective dictionary.
    # If any of the indexes (p5, p7, LIG or RT) are not found, discard the read-pair and report status.
    for read1, read2 in zip(fh_r1, fh_r2):
        # Increase the number of read-pairs read.
        qc["n_pairs"] += 1
        name_p5, name_p7, name_ligation, name_rt = None, None, None, None

        if not read1.comment:
            log.error("R1 is empty for read-pair %s and %s", read1.name, read2.name)
            sys.exit(1)

        if not read1.sequence:
            log.error("R1 is empty for read-pair %s and %s", read1.name, read2.name)
            sys.exit(1)

        # region Perform sanity checks on R1 and R2 read-pairs. --------------------------------------------------------------------------------

        # Check if R1 and R2 have the same name (correctly paired) and correct length for the first 5.000 reads.
        # Do not check for the entire file as this is too time-consuming.
        if qc["n_pairs"] <= 5000:
            if read1.name != read2.name:
                log.error("Raw R1 and R2 .fastq files are not correctly paired, read-names should be identical, exception:\n(R1) %s\n(R2) %s", read1.name, read2.name)
                sys.exit(1)

            # Check if R1 is empty.
            if not read1.sequence:
                log.error("R1 is empty for read-pair %s and %s", read1.name, read2.name)
                sys.exit(1)

            # Check if R1 is shorter than expected 34nt.
            if len(read1.sequence) != 34:
                log.error("(Read #%d) - R1 is not expected read length (34) read-pair %s (read %d)", qc["n_pairs"], read1.name, qc["n_pairs"] + 1)
                sys.exit(1)

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Retrieve the p5 and p7 indexes -------------------------------------------------------------------------------------------------

        sequence_p7_raw, sequence_p5_raw = read1.comment.split(":")[-1].split("+")

        try:
            name_p5 = dict_p5[sequence_p5_raw]
        except KeyError:
            name_p5 = find_closest_match(sequence_p5_raw, dict_p5)
            if name_p5 != None:
                qc["n_corrected_p5"] += 1
            else:
                qc["n_uncorrectable_p5"] += 1
        try:
            name_p7 = dict_p7[sequence_p7_raw]
        except KeyError:
            name_p7 = find_closest_match(sequence_p7_raw, dict_p7)
            if name_p7 != None:
                qc["n_corrected_p7"] += 1
            else:
                qc["n_uncorrectable_p7"] += 1

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Retrieve the ligation barcode ---------------------------------------------------------------------------------------------------

        # Retrieve the ligation barcode from R1.
        # This ligation barcode varies in length from 9 to 10nt.. So check for both.
        sequence_ligation = ""
        sequence_ligation_9nt = read1.sequence[0:9]
        sequence_ligation_10nt = read1.sequence[0:10]

        try:
            name_ligation = dict_ligation[sequence_ligation_9nt]
            sequence_ligation = sequence_ligation_9nt
            length_ligation = 9
        except KeyError:
            length_ligation = 10
            try:
                name_ligation = dict_ligation[sequence_ligation_10nt]
                sequence_ligation = sequence_ligation_10nt
            except KeyError:
                sequence_ligation = sequence_ligation_10nt
                name_ligation = find_closest_match(sequence_ligation, dict_ligation)
                if name_ligation != None:
                    qc["n_corrected_ligation"] += 1
                else:
                    sequence_ligation = sequence_ligation_10nt
                    qc["n_uncorrectable_ligation"] += 1

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Retrieve the RT and UMI barcodes -----------------------------------------------------------------------------------------------

        # Check length of the ligation barcode to determine the location of the other barcodes.
        if length_ligation == 10:
            # Retrieve the RT barcode from R1 (last 10 bp).
            sequence_rt_raw = read1.sequence[-10:]

            # Retrieve the UMI from R1 (next 8 bp after the ligation barcode).
            sequence_umi = read1.sequence[10:18]
        else:
            # Retrieve the RT barcode from R1 (last 10 bp, minus one).
            sequence_rt_raw = read1.sequence[-11:-1]

            # Retrieve the UMI from R1 (next 8 bp after the ligation barcode).
            sequence_umi = read1.sequence[9:17]

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Lookup the RT barcode ----------------------------------------------------------------------------------------------------------

        try:
            name_rt = dict_rt[sequence_rt_raw]
        except KeyError:
            name_rt = find_closest_match(sequence_rt_raw, dict_rt)
            if name_rt != None:
                qc["n_corrected_rt"] += 1
            else:
                qc["n_uncorrectable_rt"] += 1

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Writing output files -----------------------------------------------------------------------------------------------------------

        # If p5, p7, ligation or RT barcode could not be found, discard the read-pair.
        if name_p5 == None or name_p7 == None or name_ligation == None or name_rt == None:
            log_discarded.info(
                "\t".join(
                    (
                        str(read1.name or "?"),
                        str(name_p5 or sequence_p5_raw),
                        str(name_p7 or sequence_p7_raw),
                        str(name_ligation or sequence_ligation_10nt),
                        str(name_rt or sequence_rt_raw),
                        sequence_umi,
                    )
                )
            )
            fh_discarded_r1.write(str(read1) + "\n")
            fh_discarded_r2.write(str(read2) + "\n")

            qc["n_pairs_failure"] += 1

        else:
            # Retrieve the matching sample_name from dict_rt_barcodes.
            sample = dict_rt_barcodes[name_rt]

            # Set the read-name of R2 to include the various barcodes.
            # P5<i>-P7<i>|R2|LIG|RT_UMI
            # VH00211:236:AACK2KFM5:1:1101:28873:1000|P5A01-P7D10|LIGN|P01-D04_GCGAGCGT
            read2.set_name("{}|P5{}-P7{}|{}|{}_{}".format(read2.name, name_p5, name_p7, name_ligation, name_rt, sequence_umi))

            # Write the read-pair to the correct sample file.
            dict_fh_out[sample]["R1"].write(str(read1) + "\n")
            dict_fh_out[sample]["R2"].write(str(read2) + "\n")

            # Count as a successful read-pair.
            qc["n_pairs_success"] += 1

            # Keep track of correct read-pairs per sample.
            samples_dict[sample]["n_pairs_success"] += 1

            # Add the match_rt as key to the dictionary if it doesn't exist yet.
            if name_rt not in samples_dict[sample]["rt"]:
                samples_dict[sample]["rt"][name_rt] = 0

            samples_dict[sample]["rt"][name_rt] += 1

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Logging -------------------------------------------------------------------------------------------------------------------------

        # Print running statistics.
        if qc["n_pairs"] % 1000000 == 0:
            log.info("Processed %d read-pairs (%d discarded)", qc["n_pairs"], qc["n_pairs_failure"])

        # endregion --------------------------------------------------------------------------------------------------------------------------------

    # Close the input file handlers.
    fh_r1.close()
    fh_r2.close()

    # Close the output file handlers.
    fh_discarded_r1.close()
    fh_discarded_r2.close()

    for sample in dict_fh_out:
        for fh in dict_fh_out[sample].values():
            fh.close()

    # Save QC to pickle file.
    with open(os.path.join(path_out, "qc.pickle"), "wb") as fh:
        pickle.dump(qc, fh)
        pickle.dump(samples_dict, fh)


def printLogging_reads(log, qc, samples_dict):
    log.info("Processed %d read-pairs", qc["n_pairs"])
    log.info("  - %d read-pairs with correct p5, p7, ligation and RT barcode", qc["n_pairs_success"])
    log.info("     - %d read-pairs with corrected p5 barcode", qc["n_corrected_p5"])
    log.info("     - %d read-pairs with corrected p7 barcode", qc["n_corrected_p7"])
    log.info("     - %d read-pairs with corrected ligation barcode", qc["n_corrected_ligation"])
    log.info("     - %d read-pairs with corrected RT barcode", qc["n_corrected_rt"])
    log.info("Discarded read-pairs:")
    log.info("  - %d read-pairs without correct p5, p7, ligation and RT barcode", qc["n_pairs_failure"])
    log.info("     - %d read-pairs with uncorrectable p5 barcode", qc["n_uncorrectable_p5"])
    log.info("     - %d read-pairs with uncorrectable p7 barcode", qc["n_uncorrectable_p7"])
    log.info("     - %d read-pairs with uncorrectable ligation barcode", qc["n_uncorrectable_ligation"])
    log.info("     - %d read-pairs with uncorrectable RT barcode", qc["n_uncorrectable_rt"])

    log.info("Sample statistics:")
    for sample in samples_dict:
        log.info("  - %s: %d read-pairs", sample, samples_dict[sample]["n_pairs_success"])
        log.info("     - %d RT barcodes", len(samples_dict[sample]["rt"]))
        for rt in samples_dict[sample]["rt"]:
            log.info("        - %s: %d", rt, samples_dict[sample]["rt"][rt])


def init_logger():
    # Logging parameters.
    log = logging.getLogger(__name__)

    console = Console(force_terminal=True)
    ch = RichHandler(show_path=False, console=Console(width=255), show_time=True)
    formatter = logging.Formatter("snakemake-sciseq: %(message)s")
    ch.setFormatter(formatter)
    log.addHandler(ch)
    log.propagate = False

    # Set the verbosity level.
    log.setLevel(logging.INFO)

    return log


def main(arguments):
    description = """
    Performs 2 levels of demultiplexing on R1/R2.fastq(.gz) files generated by sci-RNA-seq v3 protocol, based on:
        1) p5/p7 PCR index, 2) RT-barcode (sample).
    
    It requires that the p5 and p7 barcodes are present in the read headers of the .fastq files (bcl2fastq):
    @<read name> 1:N:0:ACGGNNGGCC+NTCATGGNGC
                      |----p7---|+|----p5----|: p5 is reverse-complemented.

    The R1 sequence should adhere to the following scheme:
    First 9 or 10nt:  Ligation barcode
    Next 8nt:    UMI
    Next 6nt:    Primer
    Last 10nt:   RT Barcode (sample-specific)

    Anatomy of R1:
    |ACTTGATTGT| |GAGAGCTC| |CGTGAA| |AGGTTAGCAT|
    |-LIGATION-| |---UMI--| |Primer| |----RT----|
    """

    # Setup argument parser.
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--r1", required=True, type=str, help="(.fq) Input fastq (R1).")
    parser.add_argument("--r2", required=True, type=str, help="(.fq) Input fastq (R2).")
    parser.add_argument("--sequencing_name", required=True, type=str, help="(str) Sequencing sample name.")
    parser.add_argument("--samples", required=True, type=str, help="(str) Path to sample-sheet.")
    parser.add_argument("--barcodes", required=True, type=str, help="(str) Path to barcodes file.")
    parser.add_argument("--out", required=True, type=str, help="(str) Path to output directory.")

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Display help and exit.")

    # Parse arguments.
    args = parser.parse_args()

    # Initialize logging.
    log = init_logger()

    # Open sample-sheet.
    samples = pd.read_csv(args.samples, sep="\t", dtype=str)
    samples = samples.query("sequencing_name == @args.sequencing_name")

    # Open barcode-sheet.
    barcodes = pd.read_csv(args.barcodes, sep="\t", dtype=str)

    # Generate output directory if not exists.
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    # Run the program.
    sciseq_sample_demultiplexing(log=log, sequencing_name=args.sequencing_name, samples=samples, barcodes=barcodes, path_r1=args.r1, path_r2=args.r2, path_out=args.out)

    log.info("Sample demultiplexing for %s is finished.".format(args.sequencing_name.name))

    # Close loggers.
    logging.shutdown()


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()


# import cProfile
# args = parser.parse_args(
#     [
#         "--sequencing_name",
#         "230609",
#         "--samples",
#         "~/jvanriet/git/snakemake-sciseq/workflow/examples/example_samplesheet.tsv",
#         "--barcodes",
#         "~/jvanriet/git/snakemake-sciseq/workflow/examples/barcodes.tsv",
#         "--r1",
#         "/omics/groups/OE0538/internal/projects/sexomics/runJob/fastq/230609/raw/R1_1-of-20.fastq.gz",
#         "--r2",
#         "/omics/groups/OE0538/internal/projects/sexomics/runJob/fastq/230609/raw/R2_1-of-20.fastq.gz",
#         "--out",
#         "/omics/groups/OE0538/internal/projects/sexomics/runJob/fastq/230609/demux_scatter/1-of-20",
#     ]
# )
# cProfile.run('sciseq_sample_demultiplexing(log=log, sequencing_name=args.sequencing_name, samples=samples, barcodes=barcodes, path_r1=args.r1, path_r2=args.r2, path_out=args.out)')
