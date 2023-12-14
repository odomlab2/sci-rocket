import argparse
import sys
import pandas as pd

def convert_barcodes(path_starsolo_barcodes, path_barcodes, path_out):
    """
    Converts the sequence-based barcodes to their respective barcode name.

    Parameters:
        path_starsolo_barcodes (str): Path to the STARSolo barcodes.tsv.
        path_barcodes (str): Path to the barcoding scheme.
        path_out (str): Path to store converted barcodes.tsv.

    Returns:
        None
    """

    # Open barcode-sheet.
    barcodes = pd.read_csv(path_barcodes, sep="\t", dtype=str)

    # Add extra G to the ligation sequence if sequence is 9nt long.
    barcodes.loc[barcodes["type"] == "ligation", "sequence"] = barcodes.loc[barcodes["type"] == "ligation", "sequence"].apply(lambda x: x + "G" if len(x) == 9 else x)

    # Generate dictionary with barcode names as keys and sequences as values.
    barcodes_dict = dict.fromkeys(barcodes.type.unique())

    for barcode_type in barcodes_dict.keys():
        barcodes_dict[barcode_type] = dict(zip(barcodes.loc[barcodes["type"] == barcode_type, "sequence"], barcodes.loc[barcodes["type"] == barcode_type, "barcode"]))

    # The p5 index needs to be reverse complemented.
    barcodes_dict["p5"] = {k[::-1].translate(str.maketrans("ATCG", "TAGC")): v for k, v in barcodes_dict["p5"].items()}

    # Open filehandler to output.
    fh_out = open(path_out, "w")

    # Read STARSolo barcodes, line by line.
    with open(path_starsolo_barcodes, "r") as f:
        for line in f:
            # Get barcode sequences (separator is _).
            barcodes = line.strip().split("_")

            # For each barcode, get the barcode name.
            barcodes_converted = "_".join(
                [
                barcodes_dict["p7"][barcodes[0]],
                barcodes_dict["p5"][barcodes[1]],
                barcodes_dict["ligation"][barcodes[2]],
                barcodes_dict["rt"][barcodes[3]]
                ]
            )

            # Write to output.
            fh_out.write(barcodes_converted + "\n")

def main(arguments):
    # Setup argument parser.
    parser = argparse.ArgumentParser(description="Convert the sequences from STARSolo barcodes.tsv into their respective barcode name.", add_help=False)
    parser.add_argument("--starsolo_barcodes", required=True, type=str, help="(str) Path to the STARSolo barcodes.tsv.")
    parser.add_argument("--barcodes", required=True, type=str, help="(str) Path to the barcoding scheme.")
    parser.add_argument("--out", required=True, type=str, help="(str) Path to store converted barcodes.tsv.")
    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Display help and exit.")

    # Parse arguments.
    args = parser.parse_args()

    # Combine the pickles
    convert_barcodes(path_starsolo_barcodes = args.starsolo_barcodes, path_barcodes = args.barcodes, path_out = args.out)


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()

# path_starsolo_barcodes="/home/j103t/odomLab/liver_fcg/data/sci-rocket/sx42b/alignment/sx42_15_mouse_Solo.out/GeneFull/raw/barcodes.tsv"
# path_barcodes="/home/j103t/jvanriet/git/sci-rocket/workflow/examples/example_barcodes.tsv"
# path_out="~/test/barcodes_converted.tsv"