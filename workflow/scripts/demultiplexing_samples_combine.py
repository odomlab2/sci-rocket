import argparse
import sys
import os

import pickle
import json


def combine_pickles(path_scatter, path_log):
    # Find all qc.pickle files
    paths_pickle = []
    for root, dirs, files in os.walk(path_scatter):
        for file in files:
            if file.endswith("qc.pickle"):
                paths_pickle.append(os.path.join(root, file))

    # Combine the qc.pickle files
    qc = None
    sample_dict = None

    # Load the pickled dictionary
    for path_pickle in paths_pickle:
        with open(path_pickle, "rb") as handle:
            # Combine the qc dictionaries.
            qc_pickle = pickle.load(handle)

            # If the combined dictionary is empty, add the pickled dictionary
            if qc is None:
                qc = qc_pickle
            else:
                # For each key in the dictionary, add the values to the combined dictionary
                for key in qc_pickle:
                    if key == "sequencing_name":
                        continue
                    else:
                        qc[key] += qc_pickle[key]

            # Combine the sample_dict dictionaries.
            sample_dict_pickle = pickle.load(handle)

            if sample_dict is None:
                sample_dict = sample_dict_pickle
            else:
                for sample in sample_dict_pickle:
                    sample_dict[sample]["n_pairs_success"] += sample_dict_pickle[sample]["n_pairs_success"]
                    for key in sample_dict_pickle[sample]["rt"]:
                        if key in sample_dict[sample]["rt"]:
                            sample_dict[sample]["rt"][key] += sample_dict_pickle[sample]["rt"][key]
                        else:
                            sample_dict[sample]["rt"][key] = sample_dict_pickle[sample]["rt"][key]

    # Print the QC to a JSON.
    if qc != None and sample_dict != None:

        # Sort the uncorrectable barcodes by frequency.
        qc["uncorrectable_p5"] = {k: v for k, v in sorted(qc["uncorrectable_p5"].items(), key=lambda item: item[1], reverse=True)}
        qc["uncorrectable_p7"] = {k: v for k, v in sorted(qc["uncorrectable_p7"].items(), key=lambda item: item[1], reverse=True)}
        qc["uncorrectable_ligation"] = {k: v for k, v in sorted(qc["uncorrectable_ligation"].items(), key=lambda item: item[1], reverse=True)}
        qc["uncorrectable_rt"] = {k: v for k, v in sorted(qc["uncorrectable_rt"].items(), key=lambda item: item[1], reverse=True)}

        # Keep the top 50 uncorrectable barcodes.
        qc["uncorrectable_p5"] = dict(list(qc["uncorrectable_p5"].items())[:50])
        qc["uncorrectable_p7"] = dict(list(qc["uncorrectable_p7"].items())[:50])
        qc["uncorrectable_ligation"] = dict(list(qc["uncorrectable_ligation"].items())[:50])
        qc["uncorrectable_rt"] = dict(list(qc["uncorrectable_rt"].items())[:50])

        qc["sample_dict"] = sample_dict

        # Create path if it does not exist
        if not os.path.exists(os.path.dirname(path_log)):
            os.makedirs(os.path.dirname(path_log))
            
        with open(path_log, "w") as handle:
            json.dump(qc, handle, indent=4)

            # Add a description to the JSON file
            description = {}
            description["Description"] = {}
            description["Description"]["sequencing_name"] = "Name of the sequencing run."
            description["Description"]["n_pairs"] = "Total number of initial read-pairs."
            description["Description"]["n_pairs_success"] = "Total number of read-pairs with correct RT, p5, p7 and ligation barcodes."
            description["Description"]["n_pairs_failure"] = "Total number of discarded read-pairs due to various reason."
            description["Description"]["n_corrected_p5"] = "Total number of read-pairs with 1bp mismatch in p5."
            description["Description"]["n_corrected_p7"] = "Total number of read-pairs with 1bp mismatch in p7."
            description["Description"]["n_corrected_ligation"] = "Total number of read-pairs with 1bp mismatch in ligation."
            description["Description"]["n_corrected_rt"] = "Total number of read-pairs with 1bp mismatch in RT."
            description["Description"]["n_uncorrectable_p5"] = "Total number of read-pairs with >1bp mismatch in p5."
            description["Description"]["n_uncorrectable_p7"] = "Total number of read-pairs with >1bp mismatch in p7."
            description["Description"]["n_uncorrectable_ligation"] = "Total number of read-pairs with >1bp mismatch in ligation."
            description["Description"]["n_uncorrectable_rt"] = "Total number of read-pairs with >1bp mismatch in RT."
            description["Description"]["sample_dict"] = "Dictionary with sample names as keys and a dictionary showing the number of read-pairs for each RT as values."

            json.dump(description, handle, indent=4)

        # Close the file
        handle.close()

    return None


def main(arguments):
    # Setup argument parser.
    parser = argparse.ArgumentParser(description="Combine the scattered demultiplexing qc files into a JSON log", add_help=False)
    parser.add_argument("--path_scatter", required=True, type=str, help="(str) Path to the scatter_demux folder.")
    parser.add_argument("--path_log", required=True, type=str, help="(str) Path to store JSON log.")

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Display help and exit.")

    # Parse arguments.
    args = parser.parse_args()

    # Combine the pickles
    combine_pickles(args.path_scatter, args.path_log)


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()


