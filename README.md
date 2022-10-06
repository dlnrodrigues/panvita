# PanViTa - Pan Virulence and resisTance Analysis
Hello Friend!

This is a script designed to make comparisons between multiple genomes and specific databases (CARD, BacMet and VFDB).
If you have another database in mind, feel free to send a request to dlnrodrigues@ufmg.br and we will update the code.
The figure created by the script is known as a clustermap and by default uses a Euclidean distance measure to group the data.

IMPORTANT INFORMATION: The PanViTa tool is free to use under registration 20210006 of the Federal University of Minas Gerais, and is therefore the intellectual property of the institution. For this reason unofficial changes to the original script will not be allowed without prior consent.

To use the script it is possible to use the command:

```
python3 panvita.py
```

or

```
python3 panvita.py -h
```

Basically it is possible to use this tool through the command:

```
python3 panvita.py -card -vfdb -bacmet *.gbk
```

IMPORTANT INFORMATION: As input you must use .gbk, .gbf or .gbff files.

## How do I get this tool?
User, you will be able to obtain this script using the simple command:

```
git clone https://github.com/dlnrodrigues/panvita.git
cd panvita
python3 panvita.py -h
```
## What about dependencies?
Yes, we use dependencies to make this tool work. Therefore it will try to obtain all dependencies and databases automatically, however, if it is not possible you will have to obtain them by traditional means.

```
pip install wget
pip install seaborn
pip install pandas
pip install matplotlib
```

If you have conda support, the following commands may work:
```
conda install -c anaconda wget
conda install seaborn
conda install pandas
conda install -c conda-forge matplotlib
```
## How do I use this tool?
PanViTa uses some parameters to work properly. They are listed bellow.
### Databases
```
-bacmet\tAntibacterial Biocide and Metal Resistance Genes Database
-card\tComprehensive Antibiotic Resistance Database
-vfdb\tVirulence Factor Database
```
### Parameters
```
-update\tUpdate databases and dependences
-u\tSame as -update
-help\tPrint this help
-h\tSame as -help
-keep\tMaintains the protein sequences used, as well as the CDS position files
-k\tSame as -keep
-i\tMinimum identity to infer presence (default = 70)
-c\tMinimum coverage to infer presence (default = 70)
-pdf\tFigures will be saved as PDF (default)
-png\tFigures will be saved as PNG (WARNING! High memory consumption)
-g\tDownload the genomes fasta files (require CSV table from NCBI)
-a\tDownload and annote the genomes using PROKKA pipeline (require CSV table from NCBI)
-b\tDownload the genome GenBank files (require CSV table from NCBI)
-s\tKeep the locus_tag as same as the strain (require -b)
-m\tGet the metadata from BioSample IDs (require CSV table from NCBI)
```
