import argparse
import glob
import json
import os
import pandas as pd
import pickle
import sys

def write_cell_hashing_table(qc, out):
    """
    Write the following cell-based metrics:
        - Total no. of hash reads.
        - Total no. of distinct UMI  per hashing barcode.

    Parameters:
        qc (dict): Dictionary containing the QC, incl. hashing metrics.
        out (str): Path to output hashing metrics.

    Returns:
        (dict): Dictionary of the hashing metrics.
    """

    # Create a list of dictionaries.
    array_hashing = []
    
    # qc["hashing"][hashing_sample]["counts"][hashing_name][cellular_barcode] = {"umi": set(), "count": 0}
    for hashing_sample in qc["hashing"]:
        for hashing_name in qc["hashing"][hashing_sample]["counts"]:
            for cellular_barcode in qc["hashing"][hashing_sample]["counts"][hashing_name]:
                array_hashing.append(
                    {
                        "experiment_name": qc["experiment_name"],
                        "hashing_sample": hashing_sample,
                        "hashing_barcode": hashing_name,
                        "cell_barcode": cellular_barcode,
                        "count": qc["hashing"][hashing_sample]["counts"][hashing_name][cellular_barcode]["count"],
                        "n_umi": len(qc["hashing"][hashing_sample]["counts"][hashing_name][cellular_barcode]["umi"])
                    }
                )

    # Convert to pandas dataframe.
    df_hashing = pd.DataFrame(columns=["experiment_name", "hashing_sample", "hash_barcode", "cell_barcode", "count", "n_umi"], data=array_hashing)

    # Order on total hash count.
    df_hashing = df_hashing.sort_values(by=["count"], ascending=False)

    # Generate folder.
    if not os.path.exists(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out))

    # Write to file.
    df_hashing.to_csv(out, sep="\t", index=False, header=True, encoding="utf-8", mode="w")

    # Transform df_hashing into a dictionary with the hashing_sample column as key.
    dict_hashing = df_hashing.to_dict(orient="index")

    # Add the hashing_sample metrics
    for hashing_sample in qc["hashing"]:
        dict_hashing[hashing_sample]["n_correct"] = qc["hashing"][hashing_sample]["n_correct"]
        dict_hashing[hashing_sample]["n_corrected"] = qc["hashing"][hashing_sample]["n_corrected"]
        dict_hashing[hashing_sample]["n_correct_upstream"] = qc["hashing"][hashing_sample]["n_correct_upstream"]

    return dict_hashing

def combine_logs(path_pickle, path_star, path_hashing):
    """
    Combine the demuxxing logs with the STAR logs for the sci-dash.

    Parameters:
        path_pickle (str): Path to the pickled dictionaries.
        path_star (str): Path to the STAR output folder.
        path_hashing (str): Path to store the hashing metrics.

    Returns:
        qc (dict): Dictionary containing the demuxxing statistics (in JSON format).
    """

    # region Import demuxxing statistics. ---------------------------------------------------------------------------------
    qc_json = {}

    # Load the pickled dictionary
    with open(path_pickle, "rb") as handle:
        qc = pickle.load(handle)

    # endregion

    # region Calculate summary statistics. --------------------------------------------------------------------------------
    if qc != None:

        # Take over all the keys from the pickle which are not nested.
        for key in qc:
            if not isinstance(qc[key], dict):
                qc_json[key] = qc[key]

        # Determine the top recurrent uncorrectable barcodes.
        qc_json["top_uncorrectables"] = {}
        top_n = 15
        qc_json["top_uncorrectables"]["p5"] = sorted(qc["uncorrectable_p5"].items(), key=lambda x: x[1], reverse=True)[:top_n]
        qc_json["top_uncorrectables"]["p7"] = sorted(qc["uncorrectable_p7"].items(), key=lambda x: x[1], reverse=True)[:top_n]
        qc_json["top_uncorrectables"]["ligation"] = sorted(qc["uncorrectable_ligation"].items(), key=lambda x: x[1], reverse=True)[:top_n]
        qc_json["top_uncorrectables"]["rt"] = sorted(qc["uncorrectable_rt"].items(), key=lambda x: x[1], reverse=True)[:top_n]

        # Name the key/value pairs in the top_uncorrectables dictionary.
        for key in qc_json["top_uncorrectables"]:
            qc_json["top_uncorrectables"][key] = [{"barcode": barcode, "frequency": frequency} for barcode, frequency in qc_json["top_uncorrectables"][key]]

        # Name the key/value pairs in the plate counts. Split up the key in row and col.
        def transform_plate_counts(dict):
            return [{"row": key[0], "col": key[1:].lstrip("0"), "frequency": dict[key]} for key in dict]

        qc_json["p5_index_counts"] = transform_plate_counts(qc["p5_index_counts"])
        qc_json["p7_index_counts"] = transform_plate_counts(qc["p7_index_counts"])

        # Transform the rt_barcode_counts dictionary to a list of dictionaries.
        qc_json["rt_barcode_counts"] = {}

        for key in qc["rt_barcode_counts"]:
            qc_json["rt_barcode_counts"][key] = transform_plate_counts(qc["rt_barcode_counts"][key])

        # Transform the ligation_barcode_counts dictionary to a list of dictionaries.
        qc_json["ligation_barcode_counts"] = [{"barcode": barcode, "frequency": qc["ligation_barcode_counts"][barcode]} for barcode in qc["ligation_barcode_counts"]]

        # Remove all unnecessary data. ------------------------------------------------------------------------------------

        # Remove the ligation_barcode_counts with zero counts.
        qc_json["ligation_barcode_counts"] = [barcode for barcode in qc_json["ligation_barcode_counts"] if barcode["frequency"] > 0]

        # Convert the uncorrectables_sankey.
        qc_json["uncorrectables_sankey"] = [{"source": str(key), "value": qc["uncorrectables_sankey"][key]} for key in qc["uncorrectables_sankey"]]

    # endregion ----------------------------------------------------------------------------------------------------------

    # region Import STAR statistics. -------------------------------------------------------------------------------------

    # Per sample, load the STAR Log.final.out and STARSolo GeneFull summaries.
    qc_json["sample_succes"] = qc["sample_succes"]
    for sample in qc_json["sample_succes"]:

        # Find the STARsolo GeneFull summary that matches the sample name.
        path_solo = glob.glob(path_star + sample + "_*/*/Summary.csv", recursive=True)

        # Load the STARsolo GeneFull summary file and extract several statistics.
        with open(path_solo[0], "r") as handle:
            for line in handle:
                line = line.split(",")
                if line[0] == "Sequencing Saturation":
                    qc_json["sample_succes"][sample]["sequencing_saturation"] = float(line[1].strip())
                elif line[0] == "Estimated Number of Cells":
                    qc_json["sample_succes"][sample]["estimated_cells"] = int(line[1].strip())

        # Load the CellReads.stats file and extract several sample-wise statistics.
        path_cellreads = glob.glob(path_star + sample + "_*/*/CellReads.stats", recursive=True)
        df_cellreads = pd.read_csv(path_cellreads[0], sep="\t", header=0, index_col=0)

        # Summarize all cells.
        df_cellreads_summed = df_cellreads.sum(axis=0)
        df_cellreads_mean = df_cellreads.mean(axis=0)

        # Add the sample-wise statistics to the dictionary.
        qc_json["sample_succes"][sample]["total_mapped_reads"] = int(df_cellreads_summed.genomeU + df_cellreads_summed.genomeM)
        qc_json["sample_succes"][sample]["total_unique_reads"] = int(df_cellreads_summed.genomeU)
        qc_json["sample_succes"][sample]["total_multimapped_reads"] = int(df_cellreads_summed.genomeM)
        qc_json["sample_succes"][sample]["total_correct_reads_genes"] = int(df_cellreads_summed.countedU + df_cellreads_summed.countedM)
        qc_json["sample_succes"][sample]["total_exonic_reads"] = int(df_cellreads_summed.exonic)
        qc_json["sample_succes"][sample]["total_intronic_reads"] = int(df_cellreads_summed.intronic)
        qc_json["sample_succes"][sample]["total_intergenic_reads"] = int(df_cellreads_summed.genomeU + df_cellreads_summed.genomeM - df_cellreads_summed.exonic - df_cellreads_summed.intronic)
        qc_json["sample_succes"][sample]["total_mitochondrial_reads"] = int(df_cellreads_summed.mito)
        qc_json["sample_succes"][sample]["total_exonicAS_reads"] = int(df_cellreads_summed.exonicAS)
        qc_json["sample_succes"][sample]["total_intronicAS_reads"] = int(df_cellreads_summed.intronicAS)
        qc_json["sample_succes"][sample]["mean_reads_per_cell"] = int(df_cellreads_mean.genomeU + df_cellreads_mean.genomeM)
        qc_json["sample_succes"][sample]["mean_genes_per_cell"] = int(df_cellreads_mean.nGenesUnique)
        qc_json["sample_succes"][sample]["mean_umis_per_cell"] = int(df_cellreads_mean.nUMIunique)

    # endregion ----------------------------------------------------------------------------------------------------------

    # region Write / import hashing statistics. --------------------------------------------------------------------------
    if "hashing" in qc:
        df_hashing = write_cell_hashing_table(qc, path_hashing)
        qc_json["hashing"] = df_hashing
    else:
        qc_json["hashing"] = {}
    # endregion ----------------------------------------------------------------------------------------------------------

    # Return JSON-structured dict.
    return qc_json


def main(arguments):
    # Setup argument parser.
    parser = argparse.ArgumentParser(description="Combine the scattered demultiplexing files into a JSON structure for the sci-dash.", add_help=False)
    parser.add_argument("--path_pickle", required=True, type=str, help="(str) Path to demultiplexing qc pickle.")
    parser.add_argument("--path_star", required=True, type=str, help="(str) Path to the star alignment folder.")
    parser.add_argument("--path_out", required=True, type=str, help="(str) Path to store JSON structure.")
    parser.add_argument("--path_hashing", required=True, type=str, help="(str) Path to store hashing metrics (if applicable).")

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Display help and exit.")

    # Parse arguments.
    args = parser.parse_args()

    # Combine the demuxxing logs with the STAR logs for the sci-dash.
    qc_json = combine_logs(args.path_pickle, args.path_star, args.path_hashing)

    # Write the JSON structure to file.
    if not os.path.exists(os.path.dirname(args.path_out)):
        os.makedirs(os.path.dirname(args.path_out))

    with open(args.path_out, "w") as handle:
        handle.write("var data = ")
        json.dump(qc_json, handle, indent=4)

    # Close the file
    handle.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()

# path_pickle='/omics/groups/OE0538/internal/projects/liver_fcg/data/sci-rocket/sx42b/demux_reads/sx42b_qc.pickle'
# path_hashing='/home/j103t/test/sx42b_hashing.tsv'
# path_star='/omics/groups/OE0538/internal/projects/liver_fcg/data/sci-rocket/sx42b/alignment/'