import argparse
import sys
import os
import pickle
import glob
import json


def combine_logs(path_pickle, path_star):
    """
    Combine the demuxxing logs with the STAR logs for the sci-dash.

    Parameters:
        path_pickle (str): Path to the pickled dictionaries.
        path_star (str): Path to the STAR output folder.

    Returns:
        None
    """

    # region Import demuxxing statistics. ---------------------------------------------------------------------------------
    data_demux = None
    data_samples = None

    # Load the pickled dictionary

    with open(path_pickle, "rb") as handle:
        data_demux = pickle.load(handle)
        data_samples = pickle.load(handle)

    # endregion

    # region Calculate summary statistics. --------------------------------------------------------------------------------
    if data_demux != None and data_samples != None:
        # Determine the top recurrent uncorrectable barcodes.
        data_demux["top_uncorrectables"] = {}
        top_n = 15
        data_demux["top_uncorrectables"]["p5"] = sorted(data_demux["uncorrectable_p5"].items(), key=lambda x: x[1], reverse=True)[:top_n]
        data_demux["top_uncorrectables"]["p7"] = sorted(data_demux["uncorrectable_p7"].items(), key=lambda x: x[1], reverse=True)[:top_n]
        data_demux["top_uncorrectables"]["ligation"] = sorted(data_demux["uncorrectable_ligation"].items(), key=lambda x: x[1], reverse=True)[:top_n]
        data_demux["top_uncorrectables"]["rt"] = sorted(data_demux["uncorrectable_rt"].items(), key=lambda x: x[1], reverse=True)[:top_n]

        # Name the key/value pairs in the top_uncorrectables dictionary.
        for key in data_demux["top_uncorrectables"]:
            data_demux["top_uncorrectables"][key] = [{"barcode": barcode, "frequency": frequency} for barcode, frequency in data_demux["top_uncorrectables"][key]]

        # Name the key/value pairs in the plate counts. Split up the key in row and col.
        def transform_plate_counts(dict):
            return [{"row": key[0], "col": key[1:].lstrip("0"), "frequency": dict[key]} for key in dict]

        data_demux["p5_index_counts"] = transform_plate_counts(data_demux["p5_index_counts"])
        data_demux["p7_index_counts"] = transform_plate_counts(data_demux["p7_index_counts"])

        for key in data_demux["rt_barcode_counts"]:
            data_demux["rt_barcode_counts"][key] = transform_plate_counts(data_demux["rt_barcode_counts"][key])

        # Transform the ligation_barcode_counts dictionary to a list of dictionaries.
        data_demux["ligation_barcode_counts"] = [{"barcode": barcode, "frequency": data_demux["ligation_barcode_counts"][barcode]} for barcode in data_demux["ligation_barcode_counts"]]

        # Remove all unnecessary data. ------------------------------------------------------------------------------------

        del data_demux["uncorrectable_p5"]
        del data_demux["uncorrectable_p7"]
        del data_demux["uncorrectable_ligation"]
        del data_demux["uncorrectable_rt"]

        # Remove the ligation_barcode_counts with zero counts.
        data_demux["ligation_barcode_counts"] = [barcode for barcode in data_demux["ligation_barcode_counts"] if barcode["frequency"] > 0]

        # Convert the uncorrectables_sankey.
        data_demux["uncorrectables_sankey"] = [{"source": str(key), "value": data_demux["uncorrectables_sankey"][key]} for key in data_demux["uncorrectables_sankey"]]

    # endregion ----------------------------------------------------------------------------------------------------------

    # region Import STAR statistics. -------------------------------------------------------------------------------------

    # Per sample, load the STAR Log.final.out and STARSolo GeneFull summaries.
    for sample in data_samples:
        # Find the STAR log that matches the sample name.
        path_log = glob.glob(path_star + sample + "_*Log.final.out")[0]

        # Find the STARsolo GeneFull summary that matches the sample name.
        path_solo = glob.glob(path_star + sample + "_*/*/Summary.csv", recursive=True)

        # Load the STAR Log.final.out file and extract several statistics.
        with open(path_log, "r") as handle:
            for line in handle:
                line = line.lstrip().split("|")

                if line[0].startswith("Number of input reads"):
                    data_samples[sample]["total_reads"] = int(line[1].strip())
                elif line[0].startswith("Uniquely mapped reads number"):
                    data_samples[sample]["uniq_reads"] = int(line[1].strip())
                elif line[0].startswith("Number of reads mapped to multiple loci"):
                    data_samples[sample]["multimapped_reads"] = int(line[1].strip())
                elif line[0].startswith("Number of reads mapped to too many loci"):
                    data_samples[sample]["filtered_reads"] = int(line[1].strip())
                elif line[0].startswith("Number of reads unmapped"):
                    data_samples[sample]["filtered_reads"] += int(line[1].strip())
                elif line[0].startswith("Number of chimeric reads"):
                    data_samples[sample]["chimeric_reads"] = int(line[1].strip())

        # Load the STARsolo GeneFull summary file and extract several statistics.
        with open(path_solo[0], "r") as handle:
            for line in handle:
                line = line.split(",")
                if line[0] == "Sequencing Saturation":
                    data_samples[sample]["solo_sequencing_saturation"] = float(line[1].strip())
                elif line[0] == "Reads Mapped to GeneFull: Unique GeneFull":
                    data_samples[sample]["solo_uniq_gene_reads"] = data_samples[sample]["uniq_reads"] * float(line[1].strip())
                elif line[0] == "Estimated Number of Cells":
                    data_samples[sample]["solo_estimated_cells"] = int(line[1].strip())
                elif line[0] == "Unique Reads in Cells Mapped to GeneFull":
                    data_samples[sample]["solo_uniq_gene_reads_in_cells"] = int(line[1].strip())
                elif line[0] == "Mean Reads per Cell":
                    data_samples[sample]["solo_mean_reads_per_cell"] = int(line[1].strip())
                elif line[0] == "Mean UMI per Cell":
                    data_samples[sample]["solo_mean_umi_per_cell"] = int(line[1].strip())
                elif line[0] == "Mean GeneFull per Cell":
                    data_samples[sample]["solo_mean_gene_per_cell"] = int(line[1].strip())

    # endregion ----------------------------------------------------------------------------------------------------------

    # Convert to JSON structure and write to file.
    if data_demux != None and data_samples != None:
        data_demux["samples_qc"] = data_samples

    return data_demux


def main(arguments):
    # Setup argument parser.
    parser = argparse.ArgumentParser(description="Combine the scattered demultiplexing files into a JSON structure for the sci-dash.", add_help=False)
    parser.add_argument("--path_pickle", required=True, type=str, help="(str) Path to demultiplixing qc pickle.")
    parser.add_argument("--path_star", required=True, type=str, help="(str) Path to the star alignment folder.")
    parser.add_argument("--path_out", required=True, type=str, help="(str) Path to store JSON structure.")

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Display help and exit.")

    # Parse arguments.
    args = parser.parse_args()

    # Combine the demuxxing logs with the STAR logs for the sci-dash.
    qc = combine_logs(args.path_pickle, args.path_star)

    # Write the JSON structure to file.
    if not os.path.exists(os.path.dirname(args.path_out)):
        os.makedirs(os.path.dirname(args.path_out))

    with open(args.path_out, "w") as handle:
        handle.write("var data = ")
        json.dump(qc, handle, indent=4)

    # Close the file
    handle.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()
