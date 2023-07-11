__version__ = "0.1"

# Import modules.
import argparse
import gzip
import logging
import os
import pandas as pd
import pickle
import pysam
import sys

from Levenshtein import distance as hamming
from rich.console import Console
from rich.logging import RichHandler
from sanity_checks import retrieve_p5_barcodes, retrieve_p7_barcodes


def find_closest_match(sequence, comparison):
    """Find the closest match between a sequence and a dictionary of sequences.
    If there is only one match with a hamming distance of 1, return the name and sequence.
    If there are multiple matches with a hamming distance of 1, return None.

    Parameters:
        sequence (str): Sequence to compare.
        comparison (dict): Dictionary of sequences to compare to.

    Returns:
        str: Sequence of the closest match.
        str: Name of the closest match.
    """

    # Calculate the hamming distance between sequence and all keys in dict.
    distances = {k: hamming(sequence, k, score_cutoff=1) for k in comparison.keys()}
    distances = {k: v for k, v in distances.items() if v == 1}

    # If there is only one key with a distance of 1, return the name and sequence.
    if len(distances) == 1:
        sequence = next(iter(distances))
        return sequence, comparison[sequence]
    else:
        return None, None


def add_uncorrectable_sequence(sequence, dict_sequence):
    """Add an uncorrectable sequence to the dictionary of uncorrectable sequences.
    If the sequence is not in the dictionary, add it with a count of 1.

    Parameters:
        sequence (str): Sequence to add.
        dict_sequence (dict): Dictionary of uncorrectable sequences.

    Returns:
        None
    """

    if sequence not in dict_sequence:
        dict_sequence[sequence] = 1
    else:
        dict_sequence[sequence] += 1


def sciseq_sample_demultiplexing(log: logging.Logger, sequencing_name: str, samples: pd.DataFrame, barcodes: pd.DataFrame, path_r1: str, path_r2: str, path_out: str):
    """
    Performs demultiplexing of the raw fastq files based on the PCR indexes (p5, p7) and RT barcode to produce sample-specific R1 and R2 files.
    The ligation barcode can be either 9nt or 10nt long and this can affect the location of the UMI and RT barcodes.

    Read-pairs are discarded if:
        - The R1 read is empty.
        - The R1 read is shorter than 34nt.
        - The R1 read is longer than 34nt.
        - One of the barcodes (p5, p7, ligation and/or RT) is not found within R1.

    Parameters:
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

    # Open a file handler for the discarded reads log (gzip).
    path_log_discarded = os.path.join(path_out, "log_" + sequencing_name + "_discarded_reads.tsv.gz")

    try:
        fh_discarded_log = gzip.open(path_log_discarded, "wt")
    except OSError:
        log.error("Could not generate the discarded reads log file, please check the path:\n%s", path_log_discarded)
        sys.exit(1)

    # Header of discard log.
    fh_discarded_log.write("read_name\tp5\tp7\tligation\trt\tumi\n")

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

    # region Metrics for sci-dash --------------------------------------------------------------------------------------------------------------

    qc = {}
    qc["version"] = __version__
    qc["sequencing_name"] = sequencing_name
    qc["n_pairs"] = 0  # Total number of initial read-pairs.
    qc["n_pairs_success"] = 0  # Total number of read-pairs with correct RT, p5, p7 and ligation barcodes.
    qc["n_pairs_failure"] = 0  # Total number of discarded read-pairs due to various reason.

    qc["n_corrected_p5"] = 0  # Total number of read-pairs with 1bp mismatch in p5.
    qc["n_corrected_p7"] = 0  # Total number of read-pairs with 1bp mismatch in p7.
    qc["n_corrected_ligation"] = 0  # Total number of read-pairs with 1bp mismatch in ligation.
    qc["n_corrected_rt"] = 0  # Total number of read-pairs with 1bp mismatch in RT.

    # For all (succesfull) reads, keep track of the number of times each index/barcode.
    qc["p5_index_counts"] = {k: 0 for k in dict_p5.values()}
    qc["p7_index_counts"] = {k: 0 for k in dict_p7.values()}
    qc["ligation_barcode_counts"] = {k: 0 for k in dict_ligation.values()}
    qc["rt_barcode_counts"] = {}

    # Keep track of uncorrectable barcode(s) schemes within each R1.
    # Structure: ["p5", "p7", "ligation", "rt"]
    # Example: qc["uncorrectables"][True, False, True, False] to indicate that the p7 and RT barcodes are uncorrectable.
    qc["uncorrectables_sankey"] = {}

    # Add uncorrectables barcode combinations to the dict.
    for p5 in [True, False]:
        for p7 in [True, False]:
            for ligation in [True, False]:
                for rt in [True, False]:
                    qc["uncorrectables_sankey"][p5, p7, ligation, rt] = 0

    qc["n_uncorrectable_p5"] = 0  # Total number of read-pairs with >1bp mismatch in p5.
    qc["n_uncorrectable_p7"] = 0  # Total number of read-pairs with >1bp mismatch in p7.
    qc["n_uncorrectable_ligation"] = 0  # Total number of read-pairs with >1bp mismatch in ligation.
    qc["n_uncorrectable_rt"] = 0  # Total number of read-pairs with >1bp mismatch in RT.

    # Keep track of the sequence of p5, p7, ligation and RT barcodes that are not found in the respective lookup tables.
    # Structure: qc["uncorrectable_p5"]["ATCG"] = 10
    qc["uncorrectable_p5"] = {}
    qc["uncorrectable_p7"] = {}
    qc["uncorrectable_ligation"] = {}
    qc["uncorrectable_rt"] = {}

    # Sample-specific dictionary to store (key: sample_name):
    # - No. of (succesfull) read-pairs (has p5+p7+ligation+RT)
    samples_dict = {k: {"n_pairs_success": 0} for k in samples["sample_name"].unique()}

    # endregion --------------------------------------------------------------------------------------------------------------------------------

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
            sequence_p5 = sequence_p5_raw
        except KeyError:
            sequence_p5, name_p5 = find_closest_match(sequence_p5_raw, dict_p5)
            if name_p5 != None:
                qc["n_corrected_p5"] += 1
                qc["p5_index_counts"][name_p5] += 1
            else:
                qc["n_uncorrectable_p5"] += 1
                add_uncorrectable_sequence(sequence_p5_raw, qc["uncorrectable_p5"])

        try:
            name_p7 = dict_p7[sequence_p7_raw]
            sequence_p7 = sequence_p7_raw
        except KeyError:
            sequence_p7, name_p7 = find_closest_match(sequence_p7_raw, dict_p7)
            if name_p7 != None:
                qc["n_corrected_p7"] += 1
            else:
                qc["n_uncorrectable_p7"] += 1
                add_uncorrectable_sequence(sequence_p7_raw, qc["uncorrectable_p7"])

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
        except KeyError:
            try:
                name_ligation = dict_ligation[sequence_ligation_10nt]
                sequence_ligation = sequence_ligation_10nt
            except KeyError:
                sequence_ligation = sequence_ligation_10nt
                sequence_ligation, name_ligation = find_closest_match(sequence_ligation, dict_ligation)
                if name_ligation != None:
                    qc["n_corrected_ligation"] += 1
                else:
                    sequence_ligation = sequence_ligation_10nt
                    qc["n_uncorrectable_ligation"] += 1
                    add_uncorrectable_sequence(sequence_ligation, qc["uncorrectable_ligation"])

        length_ligation = len(sequence_ligation)

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
            sequence_rt = sequence_rt_raw
        except KeyError:
            sequence_rt, name_rt = find_closest_match(sequence_rt_raw, dict_rt)
            if name_rt != None:
                qc["n_corrected_rt"] += 1
            else:
                qc["n_uncorrectable_rt"] += 1
                add_uncorrectable_sequence(sequence_rt_raw, qc["uncorrectable_rt"])

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Writing output files -----------------------------------------------------------------------------------------------------------

        # If p5, p7, ligation or RT barcode could not be found, discard the read-pair.
        if name_p5 == None or name_p7 == None or name_ligation == None or name_rt == None:
            fh_discarded_log.write(
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
                + "\n"
            )
            fh_discarded_r1.write(str(read1) + "\n")
            fh_discarded_r2.write(str(read2) + "\n")

            qc["n_pairs_failure"] += 1

            # Add the uncorrectable scheme to the uncorrectable dict.
            sankey_index = [name_p5, name_p7, name_ligation, name_rt]
            sankey_index = [True if x != None else False for x in sankey_index]
            qc["uncorrectables_sankey"][tuple(sankey_index)] += 1

        else:
            # Retrieve the matching sample_name from dict_rt_barcodes.
            sample = dict_rt_barcodes[name_rt]

            # Set the sequence of R1 to the cellular barcode.
            if length_ligation == 10:
                cell_barcode = sequence_p7 + sequence_p5 + sequence_ligation + sequence_rt
                read1.sequence = cell_barcode + sequence_umi
            else:
                cell_barcode = sequence_p7 + sequence_p5 + sequence_ligation + "G" + sequence_rt
                read1.sequence = cell_barcode + sequence_umi
                

            # Set the quality of R1 to random good quality.
            read1.quality = "F" * len(read1.sequence)

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

            # Count the occurence of the barcodes.
            qc["p5_index_counts"][name_p5] += 1
            qc["p7_index_counts"][name_p7] += 1
            qc["ligation_barcode_counts"][name_ligation] += 1

            # Split up RT based on plate and well.
            name_rt_split = name_rt.split("-")

            if name_rt_split[0] not in qc["rt_barcode_counts"]:
                qc["rt_barcode_counts"][name_rt_split[0]] = {}
                qc["rt_barcode_counts"][name_rt_split[0]][name_rt_split[1]] = 1
            else:
                if name_rt_split[1] not in qc["rt_barcode_counts"][name_rt_split[0]]:
                    qc["rt_barcode_counts"][name_rt_split[0]][name_rt_split[1]] = 1
                else:
                    qc["rt_barcode_counts"][name_rt_split[0]][name_rt_split[1]] += 1

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Logging -------------------------------------------------------------------------------------------------------------------------

        # Print running statistics.
        if qc["n_pairs"] % 1000000 == 0:
            log.info("Processed %d read-pairs (%d discarded)", qc["n_pairs"], qc["n_pairs_failure"])

        # endregion --------------------------------------------------------------------------------------------------------------------------------

    # Write the barcodes to whitelist files.
    with open(os.path.join(path_out, sequencing_name + "_whitelist_p5.txt"), "w") as fh:
        for sequence, barcode in dict_p5.items():
            fh.write(sequence + "\n")

    with open(os.path.join(path_out, sequencing_name + "_whitelist_p7.txt"), "w") as fh:
        for sequence, barcode in dict_p7.items():
            fh.write(sequence + "\n")
    
    with open(os.path.join(path_out, sequencing_name + "_whitelist_ligation.txt"), "w") as fh:
        for sequence, barcode in dict_ligation.items():
            if len(sequence) == 10:
                fh.write(sequence + "\n")
            else:
                fh.write(sequence + "G" + "\n")

    with open(os.path.join(path_out, sequencing_name + "_whitelist_rt.txt"), "w") as fh:
        for sequence, barcode in dict_rt.items():
            fh.write(sequence + "\n")

    # Close the input file handlers.
    fh_r1.close()
    fh_r2.close()

    # Close the output file handlers.
    fh_discarded_r1.close()
    fh_discarded_r2.close()
    fh_discarded_log.close()

    for sample in dict_fh_out:
        for fh in dict_fh_out[sample].values():
            fh.close()

    # Save QC to pickle file.
    with open(os.path.join(path_out, "qc.pickle"), "wb") as fh:
        pickle.dump(qc, fh)
        pickle.dump(samples_dict, fh)

    log.info("Done: %d read-pairs processed (%d discarded)", qc["n_pairs"], qc["n_pairs_failure"])


def init_logger():
    # Logging parameters.
    log = logging.getLogger(__name__)

    ch = RichHandler(show_path=False, console=Console(width=255), show_time=True)
    formatter = logging.Formatter("sci-rocket: %(message)s")
    ch.setFormatter(formatter)
    log.addHandler(ch)
    log.propagate = False

    # Set the verbosity level.
    log.setLevel(logging.INFO)

    return log


def main(arguments):
    description = """
    Performs demultiplexing on R1/R2.fastq(.gz) files generated by sci-RNA-seq v3 protocol, based on:
        - RT-barcode (sample)
    
    The R1 sequence is modified to a fixed length sequence (48nt) which includes all (corrected) barcodes: p5(10nt), p7(10nt), ligation(10nt), RT(10nt) and UMI (8nt) for downstream processing.
    The read names for R2 are modified to include the barcodes and UMI.
    
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
    parser.add_argument("-v", "--version", action="version", version=__version__, help="Display version and exit.")

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


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()
