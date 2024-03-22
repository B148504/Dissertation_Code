import argparse 
from argparse import ArgumentParser
import pandas as pd
import numpy as np

def parseFasta(string, makeUppercase=False):
    splitString = string.split(">")[1:]
    names = [s.split()[0] for s in splitString]
    seqs = [s[s.index("\n"):].replace("\n","").replace(" ","") for s in splitString]
    if makeUppercase: seqs = [s.upper() for s in seqs]
    return (names,seqs)

parser=ArgumentParser()
parser.add_argument("-fasta", "--fasta_input", help="Add in the fasta file under this flag", required=True)
parser.add_argument("-soft_masked_regions_bed_output", "--output", action='store', help="Bed output of soft masked regions", required=False)
parser.add_argument("-chroms", "--chromosome_file", help="contig/scaffold/chromosomes that you want to include in this analysis", required=False)
args=parser.parse_args()

fasta_file=args.fasta_input
chrom_file=args.chromosome_file
output=args.output

with open(fasta_file) as f:
    allText = f.read()

names, seqs = parseFasta(allText)

names_starts_and_ends = []
if chrom_file:
    file1=open(chrom_file, 'r')
    Lines=file1.readlines()
    chrom=[]
    for line in Lines:
        chrom.append("{}".format(line.strip()))
    for i in range(len(names)):
        if names[i] in chrom:
            print(names[i])
            inside = False
            for j in range(len(seqs[i])):
                if inside == False:
                    if seqs[i][j].islower():
                        inside= True
                        current = [names[i], j]
                else:
                    if not seqs[i][j].islower():
                        current.append(j)
                        names_starts_and_ends.append(current)
                        inside = False
else:
    for i in range(len(names)):
        print(names[i])
        inside = False
        for j in range(len(seqs[i])):
            if inside == False:
                if seqs[i][j].islower():
                    inside= True
                    current = [names[i], j]
            else:
                if not seqs[i][j].islower():
                    current.append(j)
                    names_starts_and_ends.append(current)
                    inside = False


names_starts_and_ends_df = pd.DataFrame(names_starts_and_ends, columns=['Chrom', "Start", 'End'])

names_starts_and_ends_df.to_csv(output, sep ='\t' , index = False, header = False)