import argparse
import sys
import os
import pickle


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
        
        if key == "hashing":
            # Merge the hashing metrics.
            for hash_barcode in pickle_dict[key]:
                if hash_barcode not in combined_dict[key]:
                    combined_dict[key][hash_barcode] = pickle_dict[key][hash_barcode]
                else:
                    for cell_barcode in pickle_dict[key][hash_barcode]["counts"]:
                        if cell_barcode not in combined_dict[key][hash_barcode]["counts"]:
                            combined_dict[key][hash_barcode]["counts"][cell_barcode] = pickle_dict[key][hash_barcode]["counts"][cell_barcode]
                        else:
                            combined_dict[key][hash_barcode]["counts"][cell_barcode]["count"] += pickle_dict[key][hash_barcode]["counts"][cell_barcode]["count"]
                            combined_dict[key][hash_barcode]["counts"][cell_barcode]["umi"].update(pickle_dict[key][hash_barcode]["counts"][cell_barcode]["umi"])

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

    # region Calculate additional hashing metrics (if used). --------------------------------------------------------------

    if "hashing" in qc:
        # We calculate the following metrics:
        # - Total no. of hash reads per cell.
        # - Total no. of unique hash/UMI combinations per cell.
        # - Total no. of hash/UMI combinations per cell per hash barcode.

        # Open file-handlers to store hashing metrics.
        path_hashing = os.path.join(path_out, qc["sequencing_name"] + "_hashing_metrics.txt")
        fh_hashing = open(path_hashing, "w")

        # Write header.
        fh_hashing.write("sequencing_name\thash_barcode\tcell_barcode\tn_hash\tn_hash_umi\n")
        sequencing_name = qc["sequencing_name"]

        for hash_barcode in qc["hashing"]:
            for cell_barcode in qc["hashing"][hash_barcode]["counts"]:
                qc["hashing"][hash_barcode]["counts"][cell_barcode]["n_umi"] = len(qc["hashing"][hash_barcode]["counts"][cell_barcode]["umi"])
                del qc["hashing"][hash_barcode]["counts"][cell_barcode]["umi"]

                # Write metrics to file.
                fh_hashing.write(f"{sequencing_name}\t{hash_barcode}\t{cell_barcode}\t{qc['hashing'][hash_barcode]['counts'][cell_barcode]['count']}\t{qc['hashing'][hash_barcode]['counts'][cell_barcode]['n_umi']}\n")

    # Close file-handlers.
    fh_hashing.close()

    # endregion

    # Combined pickles.
    with open(path_out, "wb") as fh:
        pickle.dump(qc, fh)
        pickle.dump(sample_dict, fh)


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

# Read pickle
with open("e3_zhash_qc.pickle", "rb") as handle:
    qc = pickle.load(handle)
    sample_dict = pickle.load(handle)
