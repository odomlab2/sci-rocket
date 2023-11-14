import os
import sys
import argparse
import logging
from rich.console import Console
from rich.logging import RichHandler
import pickle
import pandas as pd
import numpy as np
import seaborn as sns


def write_cell_hashing_table(qc, out):
    """
    Write the following cell-based metrics:
        - Total no. of hash reads.
        - Total no. of distinct UMI  per hashing barcode.

    Parameters:
        qc (dict): Dictionary containing the QC, incl. hashing metrics.
        out (str): Path to output directory.
    
    Returns:
        None
    """

    # Generate output directory if not exists.
    if not os.path.exists(out):
        os.makedirs(out)
    
    # Open file-handlers to no. of total hashing reads and no. distinct UMI for the hashes, per cell.
    path_cell_hash_table = os.path.join(out, qc["sequencing_name"] + "_hashing_metrics.tsv")
    fh_hashing = open(path_cell_hash_table, "w")

    # Write header.
    fh_hashing.write("sequencing_name\thash_barcode\tcell_barcode\tn_hash\tn_hash_umi\n")
    sequencing_name = qc["sequencing_name"]

    for hash_barcode in qc["hashing"]:
        for cell_barcode in qc["hashing"][hash_barcode]["counts"]:
            qc["hashing"][hash_barcode]["counts"][cell_barcode]["n_umi"] = len(qc["hashing"][hash_barcode]["counts"][cell_barcode]["umi"])

            # Write metrics to file.
            fh_hashing.write(f"{sequencing_name}\t{hash_barcode}\t{cell_barcode}\t{qc['hashing'][hash_barcode]['counts'][cell_barcode]['count']}\t{qc['hashing'][hash_barcode]['counts'][cell_barcode]['n_umi']}\n")

    # Close file-handlers.
    fh_hashing.close()

def determine_background_distribution(qc, out):
    """
    Per cell, we determine the cell-based threshold for the background distribution of hashing barcodes vs. true hashing barcodes.

    Parameters:
        qc (dict): Dictionary containing the QC, incl. hashing metrics.
        out (str): Path to output directory.

    Returns:
        None
    """

    # Get all cellular barcodes.
    cell_set = set()
    
    for hash_barcode in qc["hashing"]:
        for cell_barcode in qc["hashing"][hash_barcode]["counts"]:
            cell_set.add(cell_barcode)

    # Initialize dict storing hash counts per cell.
    cell_hash_distribution = dict.fromkeys(cell_set, 0)

    # For each cell, count the number of hashing barcodes.
    for hash_barcode in qc["hashing"]:
        for cell_barcode in qc["hashing"][hash_barcode]["counts"]:
            cell_hash_distribution[cell_barcode] += qc["hashing"][hash_barcode]["counts"][cell_barcode]["count"]      
    
    # Plot distribution.
    plt = sns.displot(cell_hash_distribution.values(), kde=True, log_scale=True, bins=100)
    plt.set(xlabel="Density (cells)", ylabel=r"$\mathregular{log_{10}(No. hashing reads)}$")

    # Add log ticks to x-axis.
    plt.set(xscale="log")
    plt.set(xticks=[1, 10, 100, 1000, 10000, 100000, 1000000])

    # Add ticks to y-axis (log10).
    plt.set(yscale="log")
    plt.set(yticks=[1, 10, 100, 1000, 10000, 100000, 1000000])

    # Add light-grey grid.
    plt.set(style="whitegrid")
    
    # Change font size of axis titles + bold.
    plt.set_axis_labels(fontsize=12, fontweight="bold")

    # Save plot.
    plt.savefig(os.path.join(out, qc["sequencing_name"] + "_hashing_distribution.png"), dpi=300)
   

def init_logger():
    """
    Initializes the logger.

    Parameters:
        None

    Returns:
        log (logging.Logger): Logger object.
    """

    log = logging.getLogger(__name__)

    ch = RichHandler(show_path=False, console=Console(width=255), show_time=True)
    formatter = logging.Formatter("sci-hashing: %(message)s")
    ch.setFormatter(formatter)
    log.addHandler(ch)
    log.propagate = False

    # Set the verbosity level.
    log.setLevel(logging.INFO)

    return log


def main(arguments):
    description = """
    Calculates hashing metrics and background signal distributions.
    """

    # Setup argument parser.
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--pickle", required=True, type=str, help="(.pickle) Pickled dictionary containing the qc metrics of sci-rocket.")
    parser.add_argument("--out", required=True, type=str, help="(str) Path to output directory.")

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Display help and exit.")

    # Parse arguments.
    args = parser.parse_args()

    # Check if hashing was performed (skip script if not).
    if "hashing" in qc:

        # Initialize logging.
        log = init_logger()

        # Load the pickled dictionary
        log.info("Loading pickled dictionary.")
        with open(args.pickle, "rb") as handle:
            qc = pickle.load(handle)
            
        log.info("Calculating hashing metrics.")
        write_cell_hashing_table(qc, args.out)


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()

# # Read pickle
# with open("/omics/groups/OE0538/internal/users/l375s/hash_testing/runJob/e3_zhash/demux_reads/e3_zhash_qc.pickle", "rb") as handle:
#     qc = pickle.load(handle)

# out = "/omics/groups/OE0538/internal/users/l375s/hash_testing/runJob/e3_zhash/demux_reads/"

# qc["hashing"]["20uM_hash_P7_B1"]["counts"]["A01_A01_LIG74_P01-A01"]
# # e3.zhash        A01_A01_P01-A01_LIG74   20uM_hash_P7_B1 573

# # UMI found in bbi-sci for this sample.
# # zgrep -E "A01_A01_P01-A01_LIG74" 00.hash.gz | grep -E "20uM_hash_P7_B1" | grep -E "GGCGGCTA"

# "CCGTCTAA" in qc["hashing"]["20uM_hash_P7_B1"]["counts"]["A01_A01_LIG74_P01-A01"]["umi"]
