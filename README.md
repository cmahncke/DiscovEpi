# DiscovEpi

A Tool to automatically retrieve protein data, predict corresponding epitopes and produce potential epitope binding maps for whole proteome.

# Usage
  - Download the DiscovEpi directory for your operating system.
  - Windows: Find the DiscovEpi.exe file as ready-to-use windows executable and run it.
  - Linux: run the installation script "linux_installation.py" once and run DiscovEpi with the file main_gui.py.
  - The program requires a stable internet connection.

# Description
  This application automates epitope prediction with NetMHCpan for a set of proteins. It takes an organism's name (UniProt name not identifier, e.x. Staphylococcus aureus (strain USA300)) and eventually but not necessarily a term for the cellular location of proteins as input and returns the corresponding set of proteins. Further, one can specify the output directory. 
  
  Protein information are stored in a spreadsheet and also passed to NetMHCpan to predict putative epitopes of each protein.
  
  NetMHCpan takes an allele of the HLA locus, length of epitopes to be predicted and a threshold for binding probability of each protein as input. 
  
  Epitope binding information are stored in a spreadsheet and also passed to heatmap to produce a visual of an epitope map.
  
  The heatmap takes a length cutoff as input. It may occur that short proteins are not displayed clearly if the map tries to include a whole polyproteins of e. g. 7000 amino acids.
  
  In the end there will be two spreadsheets and a .png file of the heatmap in the given directory.

# Contact
  For any correspondance please contact cedric.mahncke@leibniz-liv.de.
