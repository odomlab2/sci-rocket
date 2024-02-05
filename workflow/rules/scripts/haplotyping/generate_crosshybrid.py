#!/usr/bin/env python3
__version__ = "0.1"

import sys
import argparse

import pysam

# Initialize logging.
from rich.console import Console
console = Console()
error_console = Console(stderr=True, style="bold red")

def main(arguments):
    
    parser = argparse.ArgumentParser(description='Determine parentally-distinct cross-hybrid haplotypes from a VCF containing one or multiple genotypes.', add_help=False)
    parser.add_argument('--haplotype', required=True, type=argparse.FileType('r'), help='VCF containing haplotype-specific variants.')
    parser.add_argument('--h1', required=True, type=str, help='First genotype (Must be present in --mgp).')
    parser.add_argument('--h2', required=False, type=str, help='(Optional) Secondary genotype (Must be present in --mgp).')
    parser.add_argument('--cross_name', required=False, default='F1', type=str, help='Name of the predicted cross-hybrid.')
    parser.add_argument('--highconfidence', required=False, action='store_true', help='Only select high-confidence genotypes (FI==1 field).')
    parser.add_argument('--no_indels', required=False, action='store_true', help='Only use SNPs, disregard InDels.')
    parser.add_argument('--out', required=True, type=argparse.FileType('w'), help='Path to output VCF with predicted F1 genotypes.')
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display help and exit.')
    parser.add_argument('-v', '--version', action="version", version="%(prog)s: v" + __version__)

    args = parser.parse_args()

    # For each potential germline variant, determine the genotype within the F1 cross-hybrid.
    # If only one haplotype is specified, the other is assumed to be identical to the reference genome (GRCm39/B6).
    generate_hybrid_genotype(args)

    
def generate_hybrid_genotype(args):
    """
    For each SNP or InDel, generate the corresponding cross-hybrid genotype.
    Only keep informative genotypes (i.e. all possible inherited alleles are distinct to each parental haplotype).
    I.e. identical haplotypes between both parents are not informative and are therefore not included in the output VCF.

    In the case of two distinct heterozygous alternative alleles (0/1 vs. 0/2), the genotype is set to all possible genotypes (0/1/2).
    
    Examples:
    H1: | 0/0 | 1/1 | 1/1 | 0/1
    H2: | 1/1 | 0/0 | 2/2 | 0/2

    F1: | 0/1 | 1/0 | 1/2 | 0/1/2

    Args:
        args (argparse.Namespace): Arguments passed to the script.
    """
    
    fVCF_haplotype = pysam.VariantFile(args.haplotype.name)

    # Check if the specified genotypes and required fields are available.
    if(args.h1 not in fVCF_haplotype.header.samples):
        error_console.print("Genotype 1 (" + args.h1 + ") not present in haplotype VCF.")
        error_console.print("Available: " + ", ".join(fVCF_haplotype.header.samples))
        sys.exit(1)

    if(args.h2 and args.h2 not in fVCF_haplotype.header.samples):
        error_console.print("Genotype 2 (" + args.h2 + ") not present in haplotype VCF.")
        error_console.print("Available: " + ", ".join(fVCF_haplotype.header.samples))
        sys.exit(1)

    # Check if the GT field is present.
    if(not "GT" in fVCF_haplotype.header.formats):
        error_console.print("GT (GENO) field not present in haplotype VCF.")
        sys.exit(1)
    
    # Subset the genotyping VCF to only contain the specified genotypes.
    if(args.h2):
        fVCF_haplotype.subset_samples([args.h1, args.h2])
    else:
        fVCF_haplotype.subset_samples([args.h1])

    # Open output VCF and add the cross-hybrid genotype.
    fVCF_out = pysam.VariantFile(args.out.name, "w", header=fVCF_haplotype.header)
    fVCF_out.header.add_meta("FORMAT", items=[("ID", "HP"), ("Number", "."), ("Type", "String"), ("Description", "Haplotype identifier.")])

    fVCF_out.header.add_sample(args.cross_name)

    # Keep track of the number of genotypes.
    qc = dict.fromkeys(['n_input', 'n_candidates', 'n_informative_H1', 'n_informative_H2', 'n_informative_hetero', 'n_informative_hetero_0/1', 'n_informative_multi', 'n_noninformative_F1', 'n_indels'], 0)

    with console.status("[bold green]Reading variants...") as status:
        for v in fVCF_haplotype:

            H1 = False
            H2 = False

            # Skip InDels if specified.
            is_indel = any([len(a) > 1 for a in v.alleles])
            if(is_indel and args.no_indels):
                continue
            
            qc['n_input'] += 1

            # Check if haplotype-specific variant (not homozygous reference for both strains).
            if(v.samples[args.h1]['GT'] != (None, None) and v.samples[args.h1]['GT'] != (0, 0)):
                qc['n_informative_H1'] += 1
                H1 = True

            if(args.h2 and (v.samples[args.h2]['GT'] != (None, None) and v.samples[args.h2]['GT'] != (0, 0))):
                qc['n_informative_H2'] += 1
                H2 = True
                
            if(not H1 and not H2):
                continue

            qc['n_candidates'] += 1

            # Count InDels passing filters.
            if(is_indel and not args.no_indels):
                qc['n_indels'] += 1

            # Retrieve the haplotypes.
            geno_h1 = v.samples[args.h1]['GT']
            if(args.h2):
                geno_h2 = v.samples[args.h2]['GT']
            else:
                geno_h2 = tuple([0] * len(v.alleles))

            # Clean haplotypes (None to 0).
            geno_h1 = tuple([0 if(g == None) else g for g in geno_h1])
            geno_h2 = tuple([0 if(g == None) else g for g in geno_h2])

            # Check if high-confidence, otherwise set to zero.
            if(args.highconfidence):
                if(v.samples[args.h1]['FI'] == None or v.samples[args.h1]['FI'] == 0):
                    geno_h1 = (0, 0)
                if(args.h2 and (v.samples[args.h2]['FI'] == None or v.samples[args.h2]['FI'] == 0)):
                    geno_h2 = (0, 0)

            # Check if each haplotype contains at least one alternative allele not found in the other haplotype.
            if(len(set(geno_h1) - set(geno_h2)) == 0 or len(set(geno_h2) - set(geno_h1)) == 0):
                qc['n_noninformative_F1'] += 1
                continue
            
            # Determine the genotype of the F1 hybrid. Keep the order of the alleles.            
            geno_h1_uniq = list(set(geno_h1))
            geno_h2_uniq = list(set(geno_h2))
            geno_f1 = tuple(geno_h1_uniq + geno_h2_uniq)
            
            # Convert to integers.
            geno_f1 = tuple([int(g) for g in geno_f1])
                
            # Heterozygous (diploid) F1 genotype.
            if(len(geno_f1) == 2 and geno_f1[0] != geno_f1[1]):
                qc['n_informative_hetero'] += 1
                
                if(geno_f1 == (0, 1) or geno_f1 == (1, 0)):
                    qc['n_informative_hetero_0/1'] += 1

            # Heterozygous (multi-allelic) F1 genotype.
            else:
                qc['n_informative_multi'] += 1

            # Generate new VCF record for output.
            r = fVCF_out.new_record(contig=v.chrom, start=v.start, stop = v.stop, alleles = v.alleles)

            # Add GT of H1/H2 and F1.
            r.samples[args.h1]['GT'] = geno_h1
            if(args.h2): 
                r.samples[args.h2]['GT'] = geno_h2
            r.samples[args.cross_name]['GT'] = geno_f1

            # Add an artifial HP of 1 variants per block.
            block_number = str(qc['n_informative_hetero'])
            r.samples[args.h1]['HP'] = "."
            if(args.h2): 
                r.samples[args.h2]['HP'] = "."
            r.samples[args.cross_name]['HP'] = tuple([block_number + '-1', block_number + '-2'])

            # Add the F1 genotype to the output VCF.
            fVCF_out.write(r)
            
        console.print("Processed {} variants; {} with at least one alternative allele in supplied haplotypes ({}%).".format(qc['n_input'], qc['n_candidates'], round(qc['n_candidates'] / qc['n_input'] * 100, 1)))
        console.print("- {} of {} were informative in {} ({}%).".format(qc['n_informative_H1'], qc['n_candidates'], args.h1, round(qc['n_informative_H1'] / qc['n_candidates'] * 100, 1)))
        if(args.h2): console.print("- {} of {} were informative in {} ({}%).".format(qc['n_informative_H2'], qc['n_candidates'], args.h2, round(qc['n_informative_H2'] / qc['n_candidates'] * 100, 1)))    

        console.print("{} of {} parental haplotypes are informative in F1 ({}%).".format(qc['n_informative_hetero'] + qc['n_informative_multi'], qc['n_candidates'], round((qc['n_informative_hetero'] + qc['n_informative_multi']) / qc['n_candidates'] * 100, 1)))
        console.print("- {} informative diploid heterozygous variants ({}%).".format(qc['n_informative_hetero'], round(qc['n_informative_hetero'] / qc['n_candidates'] * 100, 1)))
        console.print("-- {} informative diploid heterozygous variants with genotype 0/1 ({}%).".format(qc['n_informative_hetero_0/1'], round(qc['n_informative_hetero_0/1'] / qc['n_informative_hetero'] * 100, 1)))
        console.print("-- {} informative diploid heterozygous variants with other genotypes ({}%).".format(qc['n_informative_hetero'] - qc['n_informative_hetero_0/1'], round( (qc['n_informative_hetero'] - qc['n_informative_hetero_0/1']) / qc['n_informative_hetero'] * 100, 1)))
        console.print("-- {} informative variants were InDels.".format(qc['n_indels']))

        console.print("- {} informative multi-allelic variants ({}%).".format(qc['n_informative_multi'], round(qc['n_informative_multi'] / qc['n_candidates'] * 100, 1)))

        # Close and write VCF.
        fVCF_out.close()

if __name__ == '__main__':
    main(sys.argv[1:])
    sys.exit()
