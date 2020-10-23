# Compare
Hello Friend!

This is a script designed to make comparisons between multiple genomes and specific databases (CARD, BacMet and VFDB).
If you have another database in mind, feel free to send a request to dlnrodrigues@ufmg.br and we will update the code.
The figure created by the script is known as a clustermap and by default uses a Euclidean distance measure to group the data.

To use the script it is possible to use the command:

```
python3 compare.py
```

or

```
python3 compare.py -h
```

Basically it is possible to use this tool through the command:

```
python3 compare.py -card -vfdb -bacmet * .gbk
```

IMPORTANT INFORMATION: As input you must use .gbk or .gbf files.

## How to get this tool?
User, you will be able to obtain this script using the simple command:

```
git clone https://github.com/dlnrodrigues/compare.git
cd compare
python3 compare.py -h
```
## What about dependencies?
Yes, we use dependencies to make this tool work. Therefore it will try to obtain all dependencies and databases automatically, however, if it is not possible you will have to obtain them by traditional means.

```
pip install wget
pip install seaborn
pip install pandas
```

If you have conda support, the following commands may work:
```
conda install -c anaconda wget
conda install seaborn
conda install pandas
```
