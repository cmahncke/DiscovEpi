# Automated Protein Retrieval and Epitope Prediction

***

This algorithm was developed to improve the workload of research about immunogenic t-cell epitopes.
Specified only by organism, functional subcellular location and status aof review whole proteomes
and subproteomes are accessible. For each sequence of the proteomes epitopes are predicted by the
algorithms SYFPEITHI and NetMHCpan. Accessed via webservers the parameters allele, length and
threshold for binding probability can be specified. All steps are automated and require the minimum
of workload. Output will be spreadsheets for UniProt, SYFPEITHI and NetMHCpan each.

***

## Table of contents
1. [General Info](#general-info)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Output](#output)
5. [Contact](#contact)

### General Info
***
The algorithm expects commandline parameter. 

_-h_, _--help_: Calls the help function.

_-e_, _--example_: Executes the algorithm for Organism: Staphylococcus Aureus; Locations: cell membrane,
cell wall, cytoplasm and secreted; Reviewed: yes; Format: tab; Delete redundant: yes; Allele: 
HLA-A*02:01; Length: 9; SYFPEITHI theshold: 20; NetMHCpan thershold: 3.0; Output names: 
example_sa_uniprot.xlsx, example_sa_syfpeithi.xlsx, example_sa_netmhcpan.xlsx.

_-f_, _--full_: Executes the algorithm for individual data. You will be asked to input the parameters for organism, 
location(s) (split by comma), status of review and output format of the UniProt data. Organism 
parameter must exactly occure in the last entry of taxonomy in UniProt. E.x. choose '(strain USA300)'
instead of 'strain USA300' if the strain USA300 / TH1516 is not wanted. Locations are optional to set.
Only reviewed entries are certainly reliable. Tab format will collect proteom data and epitope 
predictions whereas fasta format only produces the .fsa files of the proteomes. At the end of protein
retrieval deleting redundant data for the same protein but different strains is asked. The allele 
parameter must exactly match the international nomenclature until allele position like HLA-A\*02:01.
Epitope length is restrict by the prediction algorithms, look it up beforehand SYFPEITHI's score 
ranges from 0 up to a maximum value that is dependent on the allel. For HLA-A\*02:01 e.x. it is 36.
NetMHCpan's score is the percentile rank of an epitope compared to a set of random natural peptides
and thus is affected by bias and ranges from 0 to 100. At the end of each step a spreadsheet file
named by you or standard will be produced.

### Requirements:
***
Python version:    Python 3.*

Packages: argparse, Counter (from "collections"), numpy, operator, re, requests, sys, time, xlsxwriter

Stable internet connenction.

### Installation:
***
Execute the _installation.py_ file in this directory to install all packages.

The _main.py_ file is now ready for use.

```
$ git clone https://github.com/cmahncke/cmahncke
$ cd .../path_to_the_directory/cmahncke/Code
$ python installation.py
$ python main.py -h
```

### Output:
***
Three spreadsheet files will be produced. 

Each file contains the parameters and frequencies of the sequences on its first page.

The UniProt file also contains information about ID, protein and gene name, length, mass, lineage, 
signal peptide length and amino acid sequence.

The prediction files also contain information about ID, protein name, epitope density over the whole 
sequence and signal peptide and all epitopes scoring better than the threshold with position and score.

### Contact:
Cedric Mahncke

s-cemahn@uni-greifswald.de