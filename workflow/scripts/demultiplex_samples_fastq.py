# Testing --------------------------------------------------------------------------------------------------------------------------------
import logging
import pandas as pd

from rich.console import Console
from rich.logging import RichHandler

# Logging parameters.
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

console = Console(force_terminal=True)
ch = RichHandler(show_path=False, console=console, show_time=True)
formatter = logging.Formatter("snakemake-sciseq: %(message)s")
ch.setFormatter(formatter)
log.addHandler(ch)
log.propagate = False

# Set the verbosity level.
log.setLevel(logging.INFO)


sequencing_name = "AS-951357"
path_r1 = "/omics/groups/OE0538/internal/projects/sexomics/data/rawData/fastq/230225/rawReads/AS-951357-LR-67093_R1.fastq.gz"
path_r2 = "/omics/groups/OE0538/internal/projects/sexomics/data/rawData/fastq/230225/rawReads/AS-951357-LR-67093_R2.fastq.gz"

path_barcodes = "/home/j103t/jvanriet/git/snakemake-sciseq/workflow/examples/barcodes.tsv"
path_samples = "/home/j103t/jvanriet/git/snakemake-sciseq/workflow/examples/example_samplesheet.tsv"

path_out = "/home/j103t/test/"

# Import sample sheet.
samples = pd.read_csv(path_samples, sep="\t")

# Read the barcodes.
barcodes = pd.read_csv(path_barcodes, sep="\t", comment="#")


# --------------------------------------------------------------------------------------------------------------------------------
import pysam
import sys
import difflib
import os


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
        samples (pd.DataFrame): Sample sheet.
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

    # Get the possible list of unique samples contained in the sequencing run.
    samples_sequencing = samples.query("sequencing_name == @sequencing_name")["sample_name"].unique()

    # Open file handlers for all sample-specific output files.
    dict_fh_out = {}

    for sample in samples_sequencing:
        # Generate the output paths.
        path_r1_out = os.path.join(path_out, sample + "_R1.fastq.gz")
        path_r2_out = os.path.join(path_out, sample + "_R2.fastq.gz")
        path_r1_discarded = os.path.join(path_out, sample + "_R1_discarded.fastq.gz")
        path_r2_discarded = os.path.join(path_out, sample + "_R2_discarded.fastq.gz")

        try:
            # Touch the output files to create an initial empty file.
            open(path_r1_out, "a").close()
            open(path_r2_out, "a").close()
            open(path_r1_discarded, "a").close()
            open(path_r2_discarded, "a").close()

            # Open file handlers for output R1 and R2 files.
            dict_fh_out[sample] = {}
            dict_fh_out[sample]["R1"] = pysam.FastxFile(path_r1_out)
            dict_fh_out[sample]["R2"] = pysam.FastxFile(path_r2_out)
            dict_fh_out[sample]["R1_discarded"] = pysam.FastxFile(path_r1_discarded)
            dict_fh_out[sample]["R2_discarded"] = pysam.FastxFile(path_r2_discarded)

        except OSError:
            log.error(
                "Could not generate the sample-specific demultiplexed output files, please check the paths:\n(R1) %s\n(R2) %s\n(R1 discarded) %s\n(R2 discarded) %s",
                path_r1_out,
                path_r2_out,
                path_r1_discarded,
                path_r2_discarded,
            )
            sys.exit(1)

    # endregion --------------------------------------------------------------------------------------------------------------------------------

    # region Generate RT lookup tables --------------------------------------------------------------------------------------------------------------------------------

    # Get the possible list of unique RT barcodes contained in the sequencing run.
    rt_barcodes_sequencing = samples.query("sequencing_name == @sequencing_name")["barcode_rt"].unique()

    # Select the RT barcode present in the sequencing run.
    barcodes_rt = barcodes.query("type == 'rt' & barcode in @rt_barcodes_sequencing")
    barcodes_ligation = barcodes.query("type == 'ligation'")

    # Generate dicts for fast lookup.
    dict_rt = dict(zip(barcodes_rt["sequence"], barcodes_rt["barcode"]))
    dict_ligation = dict(zip(barcodes_ligation["sequence"], barcodes_ligation["barcode"]))

    # endregion --------------------------------------------------------------------------------------------------------------------------------

    # region Initialize counters --------------------------------------------------------------------------------------------------------------------------------

    # Keep track of the number of read-pairs and QC metrics.
    n_pairs = 0  # Total number of initial read-pairs.
    n_pairs_success = 0  # Total number of read-pairs with matching RT and ligation barcodes.
    n_pairs_failure = 0  # Total number of discarded read-pairs due to various reason.

    n_empty_r1 = 0  # Total number of empty R1 reads.
    n_short_r1 = 0  # Total number of R1 reads shorter than expected 34nt.
    n_longer_r1 = 0  # Total number of R1 reads longer than expected 34nt.
    n_corrected_ligation = 0  # Total number of read-pairs with 1bp mismatch in ligation barcode.
    n_incorrected_ligation = 0  # Total number of read-pairs with mismatching ligation barcode (after correction).

    n_unmatched_pairs = 0  # Total number of read-pairs with mismatching RT barcode to any sample.
    n_corrected_pairs = 0  # Total number of read-pairs with match to sample when correcting for 1bp mismatch.

    n_contains_N_ligation = 0  # Total number of read-pairs with N in ligation barcode.
    n_contains_N_rt = 0  # Total number of read-pairs with N in RT barcode.
    n_contains_N_UMI = 0  # Total number of read-pairs with N in UMI.

    # Dictionary to store the number of (succesfull) read-pairs, no. and type of identified RT barcodes, total UMI and total unique UMI per sample.
    samples_dict = {}

    # endregion --------------------------------------------------------------------------------------------------------------------------------

    # Iterate over the read-pairs.
    for read1, read2 in zip(fh_r1, fh_r2):
        # Increase the number of read-pairs.
        n_pairs += 1

        # region Perform sanity checks on R1 and R2 read-pairs. --------------------------------------------------------------------------------

        # Check if R1 and R2 have the same name (correctly paired) for the first 1.000 reads.
        # Do not check for the entire file as this is too time-consuming.
        if n_pairs <= 1000:
            if read1.name != read2.name:
                log.error("Raw R1 and R2 .fastq files are not correctly paired, read-names should be identical, exception:\n(R1) %s\n(R2) %s", read1.name, read2.name)
                sys.exit(1)

        # Check if R1 is empty.
        if not read1.sequence:
            log.warning("R1 is empty for read-pair %s and %s", read1.name, read2.name)
            n_empty_r1 += 1
            n_pairs_failure += 1
            continue

        # Check if R1 is shorter than expected 34nt.
        if len(read1.sequence) != 34:
            log.warning("(Read #%d) - R1 is not expected read length (34) read-pair %s (read %d)", n_pairs, read1.name, n_pairs + 1)
            if len(read1.sequence) < 34:
                n_short_r1 += 1
            else:
                n_longer_r1 += 1

            n_pairs_failure += 1
            continue

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Retrieve the ligation barcode ---------------------------------------------------------------------------------------------------

        # Retrieve the ligation barcode from R1.
        # This ligation barcode varies in length from 9 to 10nt.. So check for both.
        sequence_ligation = ""
        sequence_ligation_9nt = read1.sequence[0:9]
        sequence_ligation_10nt = read1.sequence[0:10]

        try:
            match_ligation = dict_ligation[sequence_ligation_9nt]
            sequence_ligation = sequence_ligation_9nt
            log.debug("(Read #%d - R1: %s) Match ligation: %s = %s", n_pairs, read1.name, sequence_ligation_9nt, match_ligation)
        except KeyError:
            try:
                match_ligation = dict_ligation[sequence_ligation_10nt]
                sequence_ligation = sequence_ligation_10nt
                log.debug("(Read #%d - R1: %s) Match ligation: %s = %s", n_pairs, read1.name, sequence_ligation_10nt, match_ligation)
            except KeyError:
                # Try to correct for 1bp mismatch (using the 9nt).
                # If there is only a single match with a score of 0.9 or higher, use this as the corrected RT barcode.
                closest_matching = difflib.get_close_matches(sequence_ligation_10nt, dict_ligation.keys(), n=2, cutoff=0.9)

                if len(closest_matching) == 1:
                    match_ligation = dict_ligation[closest_matching[0]]
                    log.debug("(Read #%d - R1: %s) Corrected ligation: %s -> %s", n_pairs, read1.name, sequence_ligation, closest_matching[0])
                    sequence_ligation = closest_matching[0]
                    n_corrected_ligation += 1
                else:
                    sequence_ligation = sequence_ligation_10nt
                    n_incorrected_ligation += 1

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Retrieve the RT and UMI barcodes -----------------------------------------------------------------------------------------------

        # Check length of the ligation barcode to determine the location of the other barcodes.
        if len(sequence_ligation) == 10:
            # Retrieve the RT barcode from R1 (last 10 bp).
            sequence_rt = read1.sequence[-10:]

            # Retrieve the UMI from R1 (next 8 bp after the ligation barcode).
            sequence_umi = read1.sequence[10:18]
        else:
            # Retrieve the RT barcode from R1 (last 10 bp, minus one).
            sequence_rt = read1.sequence[-11:-1]

            # Retrieve the UMI from R1 (next 8 bp after the ligation barcode).
            sequence_umi = read1.sequence[9:17]

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Sanity checks on the barcodes. --------------------------------------------------------------------------------------------------

        if "N" in sequence_ligation:
            log.debug("(Read #%d - R1: %s) Contains Ns - Ligation: %s", n_pairs, read1.name, sequence_ligation)
            n_contains_N_ligation += 1

        if "N" in sequence_rt:
            log.debug("(Read #%d - R1: %s) Contains Ns - RT: %s", n_pairs, read1.name, sequence_rt)
            n_contains_N_rt += 1

        if "N" in sequence_umi:
            log.debug("(Read #%d - R1: %s) Contains Ns - UMI: %s", n_pairs, read1.name, sequence_umi)
            n_contains_N_UMI += 1

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Lookup the ligation barcode in the dictionary. --------------------------------------------------------------------------------

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Lookup the RT barcode in the dictionary. --------------------------------------------------------------------------

        try:
            match_rt = dict_rt[sequence_rt]
            log.debug("(Read #%s - R1: %s) Match RT: %s = %s", n_pairs, read1.name, sequence_rt, match_rt)
        except KeyError:
            # Try to correct for 1bp mismatch.
            # If there is only a single match with a score of 0.9 or higher, use this as the corrected RT barcode.
            closest_matching = difflib.get_close_matches(sequence_rt, dict_rt.keys(), n=2, cutoff=0.9)

            if len(closest_matching) == 1:
                match_rt = dict_rt[closest_matching[0]]
                log.debug("(Read #%d - R1: %s) Corrected RT (1nt): %s -> %s (%s)", n_pairs, read1.name, sequence_rt, closest_matching[0], match_rt)
                n_corrected_pairs += 1
            else:
                log.debug("(Read #%d - R1: %s) Unknown RT sequence - No single close match: %s", n_pairs, read1.name, sequence_rt)
                n_unmatched_pairs += 1
                continue

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region Write the read-pair to the correct sample file. --------------------------------------------------------------------------------

        # Retrieve the matching sample_name from the sample sheet based on RT barcode.
        sample_name = samples.query("barcode_rt == @match_rt")["sample_name"].values[0]
        log.debug("(Read #%d - R1: %s) Matched RT to: %s", n_pairs, read1.name, sample_name)
        n_pairs_success += 1

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        if n_pairs % 1000 == 0:
            log.info("Processed %d read-pairs", n_pairs)
            log.info("  - %d read-pairs with correct RT barcode", n_pairs_success)
            log.info("     - %d read-pairs with corrected RT barcode", n_corrected_pairs)
            log.info("     - %d read-pairs with corrected ligation barcode", n_corrected_ligation)
            log.info("     - %d read-pairs with non-matching ligation barcode", n_incorrected_ligation)
            log.info("  - %d read-pairs with no matching RT barcode", n_unmatched_pairs)
            log.info("  - %d read-pairs with empty R1", n_empty_r1)
            log.info("  - %d read-pairs with short R1", n_short_r1)
            log.info("  - %d read-pairs with long R1:", n_longer_r1)
            log.info("  - %d read-pairs with N in ligation barcode", n_contains_N_ligation)
            log.info("  - %d read-pairs with N in RT barcode", n_contains_N_rt)
            log.info("  - %d read-pairs with N in UMI", n_contains_N_UMI)
 
    # Close the input file handlers.
    fh_r1.close()
    fh_r2.close()

    # Close the output file handlers.
    for sample in dict_fh_out:
        for fh in dict_fh_out[sample].values():
            fh.close()

sciseq_sample_demultiplexing(log, sequencing_name, samples, barcodes, path_r1, path_r2, path_out)
