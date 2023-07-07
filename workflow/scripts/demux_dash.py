import argparse
import sys
import os
import pickle
import json


def combine_logs(path_pickle, path_star, path_out):
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

        def transform_plate_counts_rt(dict):
            return [{"row": key[0], "col": key[1:].lstrip("0"), "frequency": dict[key]} for key in dict]

        data_demux["p5_index_counts"] = transform_plate_counts(data_demux["p5_index_counts"])
        data_demux["p7_index_counts"] = transform_plate_counts(data_demux["p7_index_counts"])

        for key in data_demux["rt_barcode_counts"]:
            data_demux["rt_barcode_counts"][key] = transform_plate_counts(data_demux["rt_barcode_counts"][key])

        # Transform the ligation_barcode_counts dictionary to a list of dictionaries.
        data_demux["ligation_barcode_counts"] = [{"barcode": barcode, "frequency": data_demux["ligation_barcode_counts"][barcode]} for barcode in data_demux["ligation_barcode_counts"]]

        # Remove all unnecessary data. ------------------------------------------------------------------------------------
        
        del(data_demux["uncorrectable_p5"])
        del(data_demux["uncorrectable_p7"])
        del(data_demux["uncorrectable_ligation"])
        del(data_demux["uncorrectable_rt"])

        # Remove the ligation_barcode_counts with zero counts.
        data_demux["ligation_barcode_counts"] = [barcode for barcode in data_demux["ligation_barcode_counts"] if barcode["frequency"] > 0]

        # Remove the CB from the data_samples.
        for sample in data_samples:
            del(data_samples[sample]["CB"])

        # Convert the uncorrectables_sankey.        
        data_demux["uncorrectables_sankey"] = [{"source": str(key), "value": data_demux["uncorrectables_sankey"][key]} for key in data_demux["uncorrectables_sankey"]]

    # endregion ----------------------------------------------------------------------------------------------------------


    # region Import STAR statistics. -------------------------------------------------------------------------------------

    # Per sample, load the STAR Log.final.out and STARSolo GeneFull summaries.


    # endregion ----------------------------------------------------------------------------------------------------------
    # # Convert to JSON structure and write to file.
    # if qc != None and data_samples != None:
    #     qc["samples_qc"] = data_samples

    #     # Create path if it does not exist
    #     if not os.path.exists(os.path.dirname(path_out)):
    #         os.makedirs(os.path.dirname(path_out))

    #     with open(path_out, "w") as handle:
    #         handle.write("var data = ")
    #         json.dump(qc, handle, indent=4)

    #     # Close the file
    #     handle.close()

    return None


def main(arguments):
    # Setup argument parser.
    parser = argparse.ArgumentParser(description="Combine the scattered demultiplexing files into a JSON structure for the sci-dash.", add_help=False)
    parser.add_argument("--path_scatter", required=True, type=str, help="(str) Path to the scatter_demux folder.")
    parser.add_argument("--path_out", required=True, type=str, help="(str) Path to store JSON structure.")

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Display help and exit.")

    # Parse arguments.
    args = parser.parse_args()

    # Combine the demuxxing logs with the STAR logs for the sci-dash.
    combine_logs(args.path_scatter, args.path_out)

    # Write the JSON structure to file.


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()
