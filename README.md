# Automated Protein Data Retrieval and Epitope Prediction

***

This algorithm was developed to improve the workflow of research about immunogenic t-cell epitopes.
Specified only by organism and functional subcellular location whole proteomes and subproteomes are
accessible. For each sequence of the proteome epitopes are predicted by the epitope predictor 
NetMHCpan. Accessed via REST api the parameters mhc allele, epitope length and binding score 
threshold can be specified. All steps are automated and require the minimum of workload.
Output comprises spreadsheets for protein and epitope data each and an epitope map image file.

***

## Table of contents
1. [Requirements](#requirements)
2. [Installation Linux](#Installation-on-Linux-OS)
3. [Installation Windows](#Installation-on-Windows-OS)
4. [Run](#Run-DiscovEpi)
5. [Output](#output)
6. [Contact](#contact)


### Requirements:
***
Python version:    Python 3.*

Packages will be deliviered with the executable for Windows OS and retrieved automatically upon
installation on Linux OS.

A stable internet connenction is required for installation on Linux OS and running the algorithm.

### Installation on Linux OS:
***
Download the directory _DiscovEpi_linux_  from the github repository or clone the respository as
stated below:

```
$ git clone https://github.com/cmahncke/cmahncke
$ cd .../path_to_the_directory/cmahncke/DiscovEpi_linux
$ python installation_gui.py
$ python main_gui.py
```

Execute the _installation.py_ file in this directory to install all packages:
$ python installation.py

The _main.py_ file is now ready for use:
$ ./main_gui.py

### Installation on Windows OS:
***
Download the file _DiscovEpi_win_ from the github repository or clone the repository as stated
above. If you move the files always store the executable in a directory togethter with the 
_internal_ directory. Otherwise the program will not work properly.
Instead you can create a shortcut of the executable and move this to any location e. x. Desktop.

Double click on the executable application file _DiscovEpi_win_ will start the program.

### Run DiscovEpi
***
All input fields are prefilled with variables to retrieve reviewed membrane-associated proteins of
SARS-CoV-2, predict corresponding epitopes of length 9 amino acids for the MHC allele HLA-A*02:01
and create an epitope map where all proteins are displayed up to a length of 1000 amino acids.

The first step is to specify the proteome by choosing an organism and optionally a subcellular 
location the proteins should be associated with. Possible options are given in the drop-down list.
If the organism strain is not clearly selected multiple entries of teh same protein can be
retrieved. This can be prevented by answering the pop-up-window about redundancy deletion with yes.
Different input will result in an error. Additionally, protein entry can be restricted to the
status of review so that only manually curated and thus more reliable data will be retrieved.
The output directory is set to create a folder _epitope_prediction_ in the systems standard
_documents_ directory. This can be changed in the _directory_ text field at the bottom. Data will 
be retrieved upon clicking the _submit_ button.

When protein data is retrieved the _NetMHCpan_ button is enabled. The epitope prediction requires
the user  to input a MHC allele and number of amino acids of epitope length which can be chosen 
from the drop-down list and a numerical threshold to restrict retrieved epitopes to a minimum
binding score. Epitopes will be predicted for each protein automatically upon clicking the _submit_
button.

When epitope data is retrieved the _Epitope_Map_ button is enabled. The program will include every
protein on the epitope map if the number of proteins is not limited in the respective input field. 
Limiting the number of proteins will increase the vertical resolution of the map, especially when 
proteomes comprise hundreds of proteins. Likewise, the length of the displayed proteins as number 
of amino acids can be limited to increase horizontal resolution.  

### Output:
***

The protein table comprises on the first sheet log data, the input parameters and the frequencies
of each protein sequence. The second sheet is named after the subcellular and comprises data for
each protein about UniProt ID, protein name, gene name, protein length, subcellular location, 
taxonomy, signal peptide sequence position and the amino acid sequence. The filename consists of a
prefix _unp_, the organism and subcellular location. In the example _unp_SARS_CoV_2_Membrane.xlsx_.
 
The epitope table also shows log, input and epitope frequencies on the first sheet. On the second,
for each protein the corresponding epitopes are listed with the respective position in the amino
acid sequence and binding score. Additionally, there is the name and ID of the protein and three
metrics to value the proteins given. The protein score is the mean epitope binding score. The 
density is the number of predicted epitopes divided by the number of overall possible peptides of
the epitope length in the respective protein. The signal density describes the same restricted to
the signal peptide range of the protein. The filename differs to the protein file in the prefix
_nmp_ and the additional allele info. In the example _nmp_SARS_CoV_2_Membrane_HLA_A_02_01.xlsx_.

The epitope map output is a .png-file where each line of the heatmap corresponds to one proteins 
amino acid sequence. A light grey background indicates the range of the protein and darker bands
the presence of one or multiple epitopes. If the are overlapping epitopes, the epitope score is
added. The higher the epitope score at a position the darker the grey band. The filename has the
suffix _hmp_, organism and location. In the example _hmp_SARS_CoV_2_Membrane.png_.

### Contact:
Cedric Mahncke
cedric.mahncke@leibniz-liv.de
