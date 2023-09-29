import itertools
from pysam import VariantFile
import argparse 
from argparse import ArgumentParser
import numpy as np
import pandas as pd
from tqdm import tqdm



parser = ArgumentParser()
parser.add_argument("-VCF", "--VCF_file", help="Add in the VCF file under this flag", required=True)
parser.add_argument("-C", "--Chrom_info", help="Add in the VCF file under this flag", required=False)
parser.add_argument("-o", "--output", action='store', help="Assembly information, SNP, INDELs, matches, number of lines", required=False)
parser.add_argument("--Use_specific_Contigs", action='store_true', help="Use Specific contigs in -i2 textfile format with each contig on a new line", required=False)
args = parser.parse_args()


#Inputs 
vcf = VariantFile(args.VCF_file)
chrom_file=args.Chrom_info

#chrom_file = ("/Users/frankieswift/OneDrive/RA_Work/Indel_Project/data_set_chroms/GCA_916050605.1_chromosomes.txt")
#vcf = VariantFile("/Volumes/Seagate/Frankie_DTOL_lep_project/outputs/VCF/idPlaAlba1.1.vcf.gz")


#Outputs 
output = args.output


#Reading in the chrom_file and turning it into a table 
if args.input2 and args.Use_specific_Contigs:
chrom=[]
with open(chrom_file) as w:
    for line in w:
        string_chrom, numerical_chrom = line.strip().split('\t')
        if str(numerical_chrom).isdigit() == True:
            chrom.append(string_chrom)

#Set up parameters and dictionary for loop
Assembly_info = dict.fromkeys(['number_of_lines', 'matches', 'mismatches'], 0)
previous = -1


if args.input2 and args.Use_specific_Contigs:
    for chromosome in tqdm(chrom):
        for variant in vcf.fetch(chromosome):
            #Get the number of lines
            if variant.pos != previous:
                Assembly_info['number_of_lines'] += 1
            #Get the number of matches
            if variant.alts == None and variant.info['DP'] == 1:
                Assembly_info['matches'] += 1
            if variant.alts != None and variant.info ['DP'] == 1:
                Assembly_info['mismatches'] += 1
            previous = variant.pos
else:
    for variant in vcf.fetch():
        if variant.pos != previous:
                Assembly_info['number_of_lines'] += 1
        #Get the number of matches
        if variant.alts == None and variant.info['DP'] == 1:
            Assembly_info['matches'] += 1
        if variant.alts != None and variant.info ['DP'] == 1:
            Assembly_info['mismatches'] += 1
        previous = variant.pos


print(Assembly_info)

Assembly_info_df = pd.DataFrame(Assembly_info, index = [0])

Assembly_info_df.to_csv(args.output, index=False)

