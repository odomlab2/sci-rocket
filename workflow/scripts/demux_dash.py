import argparse
import sys
import os
import pickle
import json


def combine_pickle(pickle_dict, combined_dict):
    """
    Combines two pickled dictionaries.

    Parameters:
        pickle_dict (dict): Dictionary to be combined.
        combined_dict (dict): Dictionary to combine into.

    Returns:
        combined_dict (dict): Combined dictionary.
    """
    for key in pickle_dict:
        if key in ["sequencing_name", "version"]:
            continue
        else:
            # Check if value is a dictionary
            if isinstance(pickle_dict[key], dict):
                for index in pickle_dict[key]:
                    if index not in combined_dict[key]:
                        combined_dict[key][index] = pickle_dict[key][index]
                    else:
                        if isinstance(combined_dict[key][index], int):
                            combined_dict[key][index] += pickle_dict[key][index]
                        elif isinstance(pickle_dict[key], set):
                            combined_dict[key].update(pickle_dict[key])

            elif isinstance(pickle_dict[key], int):
                combined_dict[key] += pickle_dict[key]

            elif isinstance(pickle_dict[key], set):
                combined_dict[key].update(pickle_dict[key])

    return combined_dict


def combine_scattered(path_scatter, path_out):
    """
    Combine the pickled qc files into a single JSON structure for the sci-dashboard.

    Parameters:
        path_scatter (str): Path to the scatter_demux folder.
        path_out (str): Path to store JSON output.

    Returns:
        None
    """

    # Find all qc.pickle files
    paths_pickle = []
    for root, dirs, files in os.walk(path_scatter):
        for file in files:
            if file.endswith("qc.pickle"):
                paths_pickle.append(os.path.join(root, file))

    # region Combine the qc.pickle files from the scattered run. (qc, sample_dict) --------------------------------------
    qc = None
    sample_dict = None

    # Load the pickled dictionary
    for path_pickle in paths_pickle:
        with open(path_pickle, "rb") as handle:
            print(f"Loading {path_pickle}")
            
            # Combine the qc dictionaries.
            qc_pickle = pickle.load(handle)
            

            # If the combined dictionary is empty, add the pickled dictionary
            if qc is None:
                qc = qc_pickle
            else:
                qc = combine_pickle(qc_pickle, qc)

            # Combine the sample_dict dictionaries.
            sample_dict_pickle = pickle.load(handle)
            
            if sample_dict is None:
                sample_dict = sample_dict_pickle
            else:
                sample_dict = combine_pickle(sample_dict_pickle, sample_dict)

    # endregion

    # region Calculate summary statistics. --------------------------------------------------------------------------------
    if qc != None and sample_dict != None:
        # Calculate the total number of unique cells using all unique cell_barcodes.
        qc["n_cells"] = len(set().union(*[sample_dict[sample]["cells"] for sample in sample_dict]))

        for sample in sample_dict:
            # Calculate the number of unique cells.
            sample_dict[sample]["n_cells"] = len(sample_dict[sample]["cells"])

            # Calculate the number of unique UMIs.
            sample_dict[sample]["n_UMIs"] = len(set().union(*[sample_dict[sample]["cells"][cell]["UMIs"] for cell in sample_dict[sample]["cells"]]))

            # For each sample, calculate the number of cells with >=100 and >=1000 UMI.
            sample_dict[sample]["n_cells_umi_100"] = len([cell for cell in sample_dict[sample]["cells"] if sample_dict[sample]["cells"][cell]["total_UMIs"] >= 100])
            sample_dict[sample]["n_cells_umi_1000"] = len([cell for cell in sample_dict[sample]["cells"] if sample_dict[sample]["cells"][cell]["total_UMIs"] >= 1000])

            # Calculate the duplication rate for each sample by dividing the number of unique UMIs by the total number of UMIs.
            dup_rate_cell = [ 1 - len(sample_dict[sample]["cells"][cell]["UMIs"]) / sample_dict[sample]["cells"][cell]["total_UMIs"] for cell in sample_dict[sample]["cells"]]            
            sample_dict[sample]["dup_rate"] = sum(dup_rate_cell) / len(dup_rate_cell)

        # Calculate the n_cells_umi_100 and n_cells_umi_1000 over all samples.
        qc["n_cells_umi_100"] = sum([sample_dict[sample]["n_cells_umi_100"] for sample in sample_dict])
        qc["n_cells_umi_1000"] = sum([sample_dict[sample]["n_cells_umi_1000"] for sample in sample_dict])

        # Determine the top recurrent uncorrectable barcodes.
        qc["top_uncorrectables"] = {}
        top_n = 15
        qc["top_uncorrectables"]["p5"] = sorted(qc["uncorrectable_p5"].items(), key=lambda x: x[1], reverse=True)[:top_n]
        qc["top_uncorrectables"]["p7"] = sorted(qc["uncorrectable_p7"].items(), key=lambda x: x[1], reverse=True)[:top_n]
        qc["top_uncorrectables"]["ligation"] = sorted(qc["uncorrectable_ligation"].items(), key=lambda x: x[1], reverse=True)[:top_n]
        qc["top_uncorrectables"]["rt"] = sorted(qc["uncorrectable_rt"].items(), key=lambda x: x[1], reverse=True)[:top_n]

        # Name the key/value pairs in the top_uncorrectables dictionary.
        for key in qc["top_uncorrectables"]:
            qc["top_uncorrectables"][key] = [{"barcode": barcode, "frequency": frequency} for barcode, frequency in qc["top_uncorrectables"][key]]

        # Name the key/value pairs in the plate counts. Split up the key in row and col.
        def transform_plate_counts(dict):
            return [{"row": key[0], "col": key[1:].lstrip("0"), "frequency": dict[key]} for key in dict]

        def transform_plate_counts_rt(dict):
            return [{"row": key[0], "col": key[1:].lstrip("0"), "frequency": dict[key]} for key in dict]

        qc["p5_index_counts"] = transform_plate_counts(qc["p5_index_counts"])
        qc["p7_index_counts"] = transform_plate_counts(qc["p7_index_counts"])

        for key in qc["rt_barcode_counts"]:
            qc["rt_barcode_counts"][key] = transform_plate_counts(qc["rt_barcode_counts"][key])

        # Transform the ligation_barcode_counts dictionary to a list of dictionaries.
        qc["ligation_barcode_counts"] = [{"barcode": barcode, "frequency": qc["ligation_barcode_counts"][barcode]} for barcode in qc["ligation_barcode_counts"]]

        # Remove all unnecessary data. ------------------------------------------------------------------------------------
        
        del(qc["uncorrectable_p5"])
        del(qc["uncorrectable_p7"])
        del(qc["uncorrectable_ligation"])
        del(qc["uncorrectable_rt"])

        # Remove the ligation_barcode_counts with zero counts.
        qc["ligation_barcode_counts"] = [barcode for barcode in qc["ligation_barcode_counts"] if barcode["frequency"] > 0]

        for sample in sample_dict:
            del(sample_dict[sample]["cells"])

        # Convert the uncorrectables_sankey.        
        qc["uncorrectables_sankey"] = [{"source": str(key), "value": qc["uncorrectables_sankey"][key]} for key in qc["uncorrectables_sankey"]]

    # endregion ----------------------------------------------------------------------------------------------------------

    # Convert to JSON structure and write to file.
    if qc != None and sample_dict != None:
        qc["samples_qc"] = sample_dict

        # Create path if it does not exist
        if not os.path.exists(os.path.dirname(path_out)):
            os.makedirs(os.path.dirname(path_out))

        with open(path_out, "w") as handle:
            handle.write("var data = ")
            json.dump(qc, handle, indent=4)

        # Close the file
        handle.close()

    return None


def main(arguments):
    # Setup argument parser.
    parser = argparse.ArgumentParser(description="Combine the scattered demultiplexing files into a JSON structure for the sci-dash.", add_help=False)
    parser.add_argument("--path_scatter", required=True, type=str, help="(str) Path to the scatter_demux folder.")
    parser.add_argument("--path_out", required=True, type=str, help="(str) Path to store JSON structure.")

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Display help and exit.")

    # Parse arguments.
    args = parser.parse_args()

    # Combine the pickles
    combine_scattered(args.path_scatter, args.path_out)


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()
