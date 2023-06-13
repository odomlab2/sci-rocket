import logging
import sys
import pandas as pd

from rich.console import Console
from rich.logging import RichHandler

# Logging parameters.
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

console = Console(force_terminal=True)
ch = RichHandler(show_path=False, console=console, show_time=True)
formatter = logging.Formatter('snakemake-sciseq: %(message)s')
ch.setFormatter(formatter)
log.addHandler(ch)
log.propagate = False

import pysam

path_r1 = "/omics/groups/OE0538/internal/projects/sexomics/data/rawData/fastq/230225/rawReads/AS-951357-LR-67093_R1.fastq.gz"
path_r2 = "/omics/groups/OE0538/internal/projects/sexomics/data/rawData/fastq/230225/rawReads/AS-951357-LR-67093_R2.fastq.gz"

path_barcodes = "/home/j103t/jvanriet/git/snakemake-sciseq/workflow/examples/barcodes.tsv"
path_samples = "/home/j103t/jvanriet/git/snakemake-sciseq/workflow/examples/example_samplesheet.tsv"

# Import sample sheet.
log.debug("Importing sample sheet to retrieve RT barcodes for each sample.")
samples = pd.read_csv(path_samples, sep="\t")

# Read the barcodes.
log.debug("Importing sample sheet to retrieve RT barcodes for each sample.")
barcodes = pd.read_csv(path_barcodes, sep="\t", comment="#")

# Select the ligation and RT barcodes.
barcodes_rt = barcodes.query("type == 'rt'")
barcodes_ligation = barcodes.query("type == 'ligation'")

# Generate dicts for fast lookup.
dict_rt = dict(zip(barcodes_rt["barcode"], barcodes_rt["sequence"]))
dict_ligation = dict(zip(barcodes_ligation["barcode"], barcodes_ligation["sequence"]))

# Generate file handlers for the R1 and R2 .fastq files.

# Try to open the file handlers.
try:
    fh_r1 = pysam.FastxFile(path_r1)
    fh_r2 = pysam.FastxFile(path_r2)
except FileNotFoundError:
    log.error("Could not find the R1 and R2 .fastq files, please check the paths:\n(R1) %s\n(R2) %s", path_r1, path_r2)
    sys.exit(1)

# Keep track of the number of read-pairs and QC metrics.
n_pairs = 0 # Total number of initial read-pairs.
n_pairs_success = 0 # Total number of read-pairs with matching RT and ligation barcodes.
n_pairs_failure = 0 # Total number of discarded read-pairs due to various reason.

n_empty_r1 = 0 # Total number of empty R1 reads.
n_short_r1 = 0 # Total number of R1 reads shorter than expected 34nt.
n_longer_r1 = 0 # Total number of R1 reads longer than expected 34nt.

n_unmatched_pairs = 0 # Total number of read-pairs with mismatching RT barcode to any sample.
n_corrected_pairs = 0 # Total number of read-pairs with match to sample when correcting for 1bp mismatch.

# Dictionary to store the number of (succesfull) read-pairs, no. and type of identified RT barcodes, total UMI and total unique UMI per sample.
samples_dict = {}


# Iterate over the read-pairs.
for read1, read2 in zip(fh_r1, fh_r2):

    # Check if R1 and R2 have the same name (correctly paired) for the first 10.000 reads.
    # Do not check for the entire file as this is too time-consuming.
    if n_pairs <= 10000:
        if  read1.name != read2.name:
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

    # Increase the number of read-pairs.
    n_pairs += 1



# Close the file handlers.
fh_r1.close()
fh_r2.close()