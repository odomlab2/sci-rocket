import logging
import sys
import pandas as pd
import difflib

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
log.setLevel(logging.DEBUG)

# Testing --------------------------------------------------------------------------------------------------------------------------------

path_r1 = "/omics/groups/OE0538/internal/projects/sexomics/data/rawData/fastq/230225/rawReads/AS-951357-LR-67093_R1.fastq.gz"
path_r2 = "/omics/groups/OE0538/internal/projects/sexomics/data/rawData/fastq/230225/rawReads/AS-951357-LR-67093_R2.fastq.gz"

path_barcodes = "/home/j103t/jvanriet/git/snakemake-sciseq/workflow/examples/barcodes.tsv"
path_samples = "/home/j103t/jvanriet/git/snakemake-sciseq/workflow/examples/example_samplesheet.tsv"

# Import sample sheet.
samples = pd.read_csv(path_samples, sep="\t")

# Read the barcodes.
barcodes = pd.read_csv(path_barcodes, sep="\t", comment="#")

# Select the ligation and RT barcodes.
barcodes_rt = barcodes.query("type == 'rt'")
barcodes_ligation = barcodes.query("type == 'ligation'")

# Generate dicts for fast lookup.
dict_rt = dict(zip(barcodes_rt["sequence"], barcodes_rt["barcode"]))
dict_ligation = dict(zip(barcodes_ligation["sequence"], barcodes_ligation["barcode"]))

# Output files.
path_r1_out = "asd_R1.fq.gz"
path_r2_out = "asd_R2.fq.gz"
path_r1_discarded = "asd_R1_discarded.fq.gz"
path_r2_discarded = "asd_R2_discarded.fq.gz"

# --------------------------------------------------------------------------------------------------------------------------------
import pysam


def sciseq_sample_demultiplexing(log, sample, samples, barcodes, path_r1, path_r2, path_r1_out, path_r2_out, path_r1_discarded, path_r2_discarded):
    """
    Performs demultiplexing of the raw fastq files based on the RT and ligation barcodes to produce sample-specific R1 and R2 files.

    Args:
        log (logging.Logger): Logger.
        sample (str): Sample name.
        samples (pd.DataFrame): Sample sheet.
        barcodes (pd.DataFrame): Barcode sheet.
        path_r1 (str): Path to R1 fastq file.
        path_r2 (str): Path to R2 fastq file.
        path_r1_out (str): Path to output R1 fastq file.
        path_r2_out (str): Path to output R2 fastq file.
        path_r1_discarded (str): Path to output discarded R1 fastq file.
        path_r2_discarded (str): Path to output discarded R2 fastq file.

    Returns:
        None
    """
    pass


# Open file handlers for input R1 and R2 files.
try:
    fh_r1 = pysam.FastxFile(path_r1)
    fh_r2 = pysam.FastxFile(path_r2)
except OSError:
    log.error("Could not find the R1 and R2 .fastq files, please check the paths:\n(R1) %s\n(R2) %s", path_r1, path_r2)
    sys.exit(1)

# Open file handlers for output files.
try:
    fh_r1_out = pysam.FastxFile(path_r1_out)
    fh_r2_out = pysam.FastxFile(path_r2_out)
    fh_r1_discarded_out = pysam.FastxFile(path_r1_discarded)
    fh_r2_discarded_out = pysam.FastxFile(path_r1_discarded)
except OSError:
    log.error(
        "Could not generate the sample-specific demultiplexed output files, please check the paths:\n(R1) %s\n(R2) %s\n(R1 discarded) %s\n(R2 discarded) %s",
        path_r1_out,
        path_r2_out,
        path_r1_discarded,
        path_r2_discarded,
    )
    sys.exit(1)

# Keep track of the number of read-pairs and QC metrics.
n_pairs = 0  # Total number of initial read-pairs.
n_pairs_success = 0  # Total number of read-pairs with matching RT and ligation barcodes.
n_pairs_failure = 0  # Total number of discarded read-pairs due to various reason.

n_empty_r1 = 0  # Total number of empty R1 reads.
n_short_r1 = 0  # Total number of R1 reads shorter than expected 34nt.
n_longer_r1 = 0  # Total number of R1 reads longer than expected 34nt.

n_unmatched_pairs = 0  # Total number of read-pairs with mismatching RT barcode to any sample.
n_corrected_pairs = 0  # Total number of read-pairs with match to sample when correcting for 1bp mismatch.

n_contains_N_ligation = 0  # Total number of read-pairs with N in ligation barcode.
n_contains_N_rt = 0  # Total number of read-pairs with N in RT barcode.
n_contains_N_UMI = 0  # Total number of read-pairs with N in UMI.

# Dictionary to store the number of (succesfull) read-pairs, no. and type of identified RT barcodes, total UMI and total unique UMI per sample.
samples_dict = {}

log.info("Starting sample-based demultiplexing of R1/R2:\n(R1) %s\n(R2) %s", path_r1, path_r2)

# Iterate over the read-pairs.
for read1, read2 in zip(fh_r1, fh_r2):
    # region Perform sanity checks on R1 and R2 read-pairs. --------------------------------------------------------------------------------

    # Check if R1 and R2 have the same name (correctly paired) for the first 10.000 reads.
    # Do not check for the entire file as this is too time-consuming.
    if n_pairs <= 10000:
        if read1.name != read2.name:
            log.error("Raw R1 and R2 .fastq files are not correctly paired, read-names should be identical, exception:\n(R1) %s\n(R2) %s", read1.name, read2.name)
            sys.exit(1)

    # Check if R1 is empty.
    if not read1.sequence:
        log.warning("R1 is empty for read-pair %s and %s", read1.name, read2.name)
        n_empty_r1 += 1
        n_pairs_failure += 1
        break

    # Check if R1 is shorter than expected 34nt.
    if len(read1.sequence) != 34:
        log.warning("R1 is not expected read length (34) read-pair %s (read %d)", read1.name, n_pairs + 1)
        if len(read1.sequence) < 34:
            n_short_r1 += 1
        else:
            n_longer_r1 += 1

        n_pairs_failure += 1
        break

    # endregion --------------------------------------------------------------------------------------------------------------------------------

    # region Dissect R1 to determine the ligation barcode, RT barcode and UMI. ---------------------------------------------------------------

    # Retrieve the ligation barcode from R1 (first 10 bp).
    sequence_ligation = read1.sequence[0:10]

    # Retrieve the RT barcode from R1 (last 10 bp).
    sequence_rt = read1.sequence[-10:]

    # Retrieve the UMI from R1 (next 8 bp after the ligation barcode).
    sequence_umi = read1.sequence[10:18]

    # Check Ns in the ligation barcode.
    if "N" in sequence_ligation:
        log.warning("(R1: %s) Contains Ns - Ligation: %s", read1.name, sequence_ligation)
        n_contains_N_ligation += 1

    if "N" in sequence_rt:
        log.warning("(R1: %s) Contains Ns - RT: %s", read1.name, sequence_rt)
        n_contains_N_rt += 1

    if "N" in sequence_umi:
        log.warning("(R1: %s) Contains Ns - UMI: %s", read1.name, sequence_umi)
        n_contains_N_UMI += 1

    # endregion --------------------------------------------------------------------------------------------------------------------------------

    # region Lookup the RT barcode in the dictionary. --------------------------------------------------------------------------

    try:
        match_rt = dict_rt[sequence_rt]
        log.debug("(R1: %s) Match RT: %s = %s", read1.name, sequence_rt, match_rt)
    except KeyError:
        log.warning("(R1: %s) Unknown sequence - RT: %s", read1.name, sequence_rt)

        # Try to correct for 1bp mismatch.
        # If there is only a single match with a score of 0.9 or higher, use this as the corrected RT barcode.
        closest_matching = difflib.get_close_matches(sequence_rt, dict_rt.keys(), n=2, cutoff=0.9)

        if len(closest_matching) == 1:
            log.debug("(R1: %s) Corrected RT: %s -> %s", read1.name, sequence_rt, closest_matching[0])
            match_rt = dict_rt[closest_matching[0]]
            log.debug("(R1: %s) Match RT: %s = %s", read1.name, sequence_rt, match_rt)
            n_corrected_pairs += 1
        else:
            n_unmatched_pairs += 1
            break

    # endregion --------------------------------------------------------------------------------------------------------------------------------

    # region Write the read-pair to the correct sample file. --------------------------------------------------------------------------------

    # Retrieve the matching sample_name from the sample sheet based on RT barcode.
    sample_name = samples.query("barcode_rt == @match_rt")["sample_name"].values[0]

    # endregion --------------------------------------------------------------------------------------------------------------------------------

    # Increase the number of read-pairs.
    n_pairs += 1


# Close the file handlers.
fh_r1.close()
fh_r2.close()
