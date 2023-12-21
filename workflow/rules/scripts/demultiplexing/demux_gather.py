import argparse
import sys
import os
import pickle

from collections import defaultdict

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
        # Skip the experiment_name and version as these are already defined and identical.
        if key in ["experiment_name", "version"]:
            continue

        # Merge the hashing metrics (if applicable)
        if key == "hashing":
            for hashing_sample in pickle_dict["hashing"]:
                if hashing_sample not in combined_dict["hashing"]:
                    combined_dict["hashing"][hashing_sample] = pickle_dict["hashing"][hashing_sample]
                else:
                    # Merge the overall hashing metrics.
                    combined_dict["hashing"][hashing_sample]["n_correct"] += pickle_dict["hashing"][hashing_sample]["n_correct"]
                    combined_dict["hashing"][hashing_sample]["n_corrected"] += pickle_dict["hashing"][hashing_sample]["n_corrected"]
                    combined_dict["hashing"][hashing_sample]["n_correct_upstream"] += pickle_dict["hashing"][hashing_sample]["n_correct_upstream"]

                    # Merge the hashing metrics per hashing_name.
                    for hashing_name in pickle_dict["hashing"][hashing_sample]["counts"]:
                        if hashing_name not in combined_dict["hashing"][hashing_sample]["counts"]:
                            combined_dict["hashing"][hashing_sample]["counts"][hashing_name] = pickle_dict["hashing"][hashing_sample]["counts"][hashing_name]
                        else:
                            for cellular_barcode in pickle_dict["hashing"][hashing_sample]["counts"][hashing_name]:
                                if cellular_barcode not in combined_dict["hashing"][hashing_sample]["counts"][hashing_name]:
                                    combined_dict["hashing"][hashing_sample]["counts"][hashing_name][cellular_barcode] = pickle_dict["hashing"][hashing_sample]["counts"][hashing_name][cellular_barcode]
                                else:
                                    combined_dict["hashing"][hashing_sample]["counts"][hashing_name][cellular_barcode]["umi"].update(pickle_dict["hashing"][hashing_sample]["counts"][hashing_name][cellular_barcode]["umi"])
                                    combined_dict["hashing"][hashing_sample]["counts"][hashing_name][cellular_barcode]["count"] += pickle_dict["hashing"][hashing_sample]["counts"][hashing_name][cellular_barcode]["count"]
                                    
        elif key == "sample_succes":
            for sample in pickle_dict["sample_succes"]:
                if sample not in combined_dict["sample_succes"]:
                    combined_dict["sample_succes"][sample] = pickle_dict["sample_succes"][sample]
                else:
                    combined_dict["sample_succes"][sample]["n_pairs_success"] += pickle_dict["sample_succes"][sample]["n_pairs_success"]
        
        # Merge everything else.
        else:
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


def combine_scattered(path_demux_scatter, path_out):
    """
    Combine the pickled qc files into a single pickle.

    Parameters:
        path_demux_scatter (str): Path to the scatter_demux folder.
        path_out (str): Path to store pickle object.

    Returns:
        None
    """

    # Find all qc.pickle files
    paths_pickle = []
    for root, dirs, files in os.walk(path_demux_scatter):
        for file in files:
            if file.endswith("qc.pickle"):
                paths_pickle.append(os.path.join(root, file))

    # region Combine the qc.pickle files from the scattered run. (qc, sample_dict) --------------------------------------
    qc = None

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

    # endregion

    # Combined pickles.
    pickle.dump(qc, open(path_out, "wb"))


def main(arguments):
    # Setup argument parser.
    parser = argparse.ArgumentParser(description="Combine the scattered pickles into a single pickle.", add_help=False)
    parser.add_argument("--path_demux_scatter", required=True, type=str, help="(str) Path to the scatter_demux folder.")
    parser.add_argument("--path_out", required=True, type=str, help="(str) Path to store combined pickle.")

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Display help and exit.")

    # Parse arguments.
    args = parser.parse_args()

    # Combine the pickles
    combine_scattered(args.path_demux_scatter, args.path_out)


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()
