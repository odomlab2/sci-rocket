import pandas as pd
import argparse
import sys


def main(arguments):
    
    parser = argparse.ArgumentParser(description='Join the H1/H2/UA count tables of nmi_tools.', add_help=False)
    parser.add_argument('--h1', required=True, type=argparse.FileType('r'), help='Count table of H1.')
    parser.add_argument('--h2', required=True, type=argparse.FileType('r'), help='Count table of H2.')
    parser.add_argument('--ua', required=True, type=argparse.FileType('r'), help='Count table of UA.')
    parser.add_argument("-o", "--output", action='store', type=argparse.FileType('w'), help='Path to write the joined table.')
    args = parser.parse_args()

    # For each potential germline variant, determine the genotype within the F1 cross-hybrid.
    # If only one haplotype is specified, the other is assumed to be identical to the reference genome (GRCm39/B6).
    join_tables(args)

def join_tables(args):
    file1 = pd.read_table(args.h1.name, sep='\t', header=0, names=['gene', 'cell', 'H1'])
    file2 = pd.read_table(args.h2.name, sep='\t', header=0, names=['gene', 'cell', 'H2'])

    # Join H1/H2.
    joined = pd.merge(file1, file2, on=['gene', 'cell'], how='outer')
    del(file1, file2)

    # Join UA.
    file3 = pd.read_table(args.ua.name, sep='\t', header=0, names=['gene', 'cell', 'UA'])
    joined = pd.merge(joined, file3, on=['gene', 'cell'], how='outer')

    # Write output.
    joined.to_csv(args.output.name, sep='\t', index=False)

if __name__ == '__main__':
    main(sys.argv[1:])
    sys.exit()
