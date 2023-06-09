# Dissertation
Scripts I have written in the process my dissertation: Characterising indel diversity in Lepidoptera to infer the evolutionary effects of genetic drift and selection

## Cloning and Requirements
To download the repository + install reqirements use:

```
github clone https://github.com/B148504/Dissertation_Code.git

cd Dissertation_Code

python3 -m pip install -r Dissertation_requirements.txt
```
# Dissertation_full_shell_pipeline.sh
This is the full pipeline of the project in a blank form. The whole script was automated so the path names are arbritary. If any R steps were required before a step could be done this is outlined in the script.

# Assembly_information.py
A script that pulls out number of lines, SNPs, homologous site and Indels for all or a specified list of contigs. Indel positions are recorded in a [Chromosome, Start, End]. Indels of lengths 1-100 are binned for every length. All Indels are also binned by length in ranges 1-100, 101-200 etc. All specific Indel lengths are also recorded. 

### Usage:
```
python3 Assembly_info.py -inputs -outputs
```
The only required input is the VCF file. If you wish to specify specific contigs or chromosomes to be run use -i2 and use a textfile with each contig name on a new row.

### Inputs
| Flags | Input | Required | 
|-|-|-|
|-i| VCF file| Yes |
|-i2|Textfile specifying which contigs to use| No |
|--Use_specific_Contigs|If using -i2 also include this flag | No |

### Outputs
| Flag | Output | Printed to terminal | Output type | Required | 
|-|-|-|-|-|
|-o | Number of lines, number of homozygous sites, number of SNPs, number of Indels| Yes | Textfile | No |
|-o2 | Chromosome, start, end position of all Indels | No | Bed file | No |
|-o3 | Indels of lengths 1-100 binned by length | Yes | Textfile | No |
|-o4 | All Indels binned by length in ranges | Yes | Textfile | No |  

If no output file is specified with a flag the outputs with Yes in printed to terminal will print outputs to terminal but no outputs will be saved

# Indel_sorter.py 
This script takes Indels start and end positions to discern Indel length and then bins Indels between 1-100 and then all Indel lengths are binned in ranges e.g. 1-100, 101-200 etc. 

### Usage:
```
python3 exon_R_indel_sorter.py -inputs -outputs 
```

### Inputs
| Flag | Input | Required | 
|-|-|-|
|-i| Indel positions textfile with Chromosome, Start, End format| Yes |

### Outputs
| Flag | Output | Printed to terminal | Output type | Required | 
|-|-|-|-|-|
|-o | Indels of lengths 1-100 binned by length  | Yes | Textfile | No |
|-o2 | All Indels binned by length in ranges | Yes | Textfile | No |  

If no output file is specified with a flag the outputs with Yes in printed to terminal will print outputs to terminal but no outputs will be saved
# Repeat_finder.py
A script that finds masked repetitive regions in fasta files

### Usage:
```
python3 Repeat_finder.py -inputs -outputs 
```
### Inputs
| Flag | Input | Required | 
|-|-|-|
|-i| reference fasta file| Yes |
|-i2|Textfile specifying which contigs to use| No |

### Outputs
| Flag | Output | Printed to terminal | Output type | Required | 
|-|-|-|-|-|
|-o | Chromosome, start and end position of masked repeat | No | Textfile | No |

If no output file is specified with a flag the outputs with Yes in printed to terminal will print outputs to terminal but no outputs will be saved

# SNP_finder.py
A script that finds polymorphic bases between two haplotypes on specified contigs

### Usage:
```
python3 SNP_finder.py -inputs -outputs 
```
### Inputs
| Flag | Input | Required | 
|-|-|-|
|-i| VCF file| Yes |
|-i2|Textfile specifying which contigs to use| No |
|--Use_specific_Contigs|If using -i2 also include this flag | No |

### Outputs
| Flag | Output | Printed to terminal | Output type | Required | 
|-|-|-|-|-|
|-o | Number of lines, number of homozygous sites, number of SNPs| Yes | Textfile | No |

If no output file is specified with a flag the outputs with Yes in printed to terminal will print outputs to terminal but no outputs will be saved

