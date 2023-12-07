import pandas as pd
import argparse
import sys


def main(arguments):
    
    parser = argparse.ArgumentParser(description='Join the H1/H2/UA count tables of nmi_tools.', add_help=False)
    parser.add_argument('--h1', required=True, type=argparse.FileType('r'), help='Count table of H1.')
    parser.add_argument('--h2', required=True, type=argparse.FileType('r'), help='Count table of H2.')
    parser.add_argument('--ua', required=True, type=argparse.FileType('r'), help='Count table of UA.')
    parser.add_argument("--barcodes", required=True, type=str, help="(str) Path to the barcoding scheme.")
    parser.add_argument("-o", "--output", action='store', type=argparse.FileType('w'), help='Path to write the joined table.')
    args = parser.parse_args()

    # For each potential germline variant, determine the genotype within the F1 cross-hybrid.
    # If only one haplotype is specified, the other is assumed to be identical to the reference genome (GRCm39/B6).
    join_tables(args)

def join_tables(args):

    # Join H1/H2/UA. -----------------------------------------------------------

    file1 = pd.read_table(args.h1.name, sep='\t', header=0, names=['gene', 'cell', 'H1'])
    file2 = pd.read_table(args.h2.name, sep='\t', header=0, names=['gene', 'cell', 'H2'])

    # Join H1/H2.
    joined = pd.merge(file1, file2, on=['gene', 'cell'], how='outer')
    del(file1, file2)

    # Join UA.
    file3 = pd.read_table(args.ua.name, sep='\t', header=0, names=['gene', 'cell', 'UA'])
    joined = pd.merge(joined, file3, on=['gene', 'cell'], how='outer')

    # Remove rows without a gene name.
    joined = joined[joined["gene"] != "-"]

    # Fill NA with 0.
    joined = joined.fillna(0)
    
    # Convert barcodes. --------------------------------------------------------

    barcodes = pd.read_csv(args.barcodes, sep="\t", dtype=str)

    # Add extra G to the ligation sequence if sequence is 9nt long.
    barcodes.loc[barcodes["type"] == "ligation", "sequence"] = barcodes.loc[barcodes["type"] == "ligation", "sequence"].apply(lambda x: x + "G" if len(x) == 9 else x)

    # Generate dictionary with barcode names as keys and sequences as values.
    barcodes_dict = dict.fromkeys(barcodes.type.unique())

    for barcode_type in barcodes_dict.keys():
        barcodes_dict[barcode_type] = dict(zip(barcodes.loc[barcodes["type"] == barcode_type, "sequence"], barcodes.loc[barcodes["type"] == barcode_type, "barcode"]))

    # The p5 index needs to be reverse complemented.
    barcodes_dict["p5"] = {k[::-1].translate(str.maketrans("ATCG", "TAGC")): v for k, v in barcodes_dict["p5"].items()}

    # Split cell barcodes into p7, p5, ligation and RT using _.
    joined["p7"] = joined["cell"].apply(lambda x: x.split("_")[0])
    joined["p5"] = joined["cell"].apply(lambda x: x.split("_")[1])
    joined["ligation"] = joined["cell"].apply(lambda x: x.split("_")[2])
    joined["rt"] = joined["cell"].apply(lambda x: x.split("_")[3])

    # For each barcode, get the barcode name.
    joined["p7"] = joined["p7"].apply(lambda x: barcodes_dict["p7"][x])
    joined["p5"] = joined["p5"].apply(lambda x: barcodes_dict["p5"][x])
    joined["ligation"] = joined["ligation"].apply(lambda x: barcodes_dict["ligation"][x])
    joined["rt"] = joined["rt"].apply(lambda x: barcodes_dict["rt"][x])

    # Combine the barcodes into one column (cell).
    joined["cell"] = joined[["p7", "p5", "ligation", "rt"]].apply(lambda x: "_".join(x), axis=1)

    # Remove the separate barcode columns.
    joined = joined.drop(["p7", "p5", "ligation", "rt"], axis=1)
    
    # Write output. ------------------------------------------------------------

    joined.to_csv(args.output.name, sep='\t', index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
    sys.exit()
