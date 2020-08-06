#!/usr/bin/env python3

import sys
import os
import shutil
import math
try:
	import wget
except:
	try:
		print("You may not have \'wget\' from python3.\nWe will try to install using pip...")
		os.system("pip install wget")
		import wget
	except:
		print("""It didn't work. Try using ONE of the commands below and then run this script again.

			pip install wget --user

			If you use conda:
			conda install -c anaconda wget
			conda install -c conda-forge python-wget
			conda install -c conda-forge/label/gcc7 python-wget""")
		exit()

if (("-card" not in sys.argv) and ("-bacmet" not in sys.argv) and ("-vfdb" not in sys.argv)) or ("-h" in sys.argv) or ("-help" in sys.argv):
	print('''
Hello user!

This script has the function of comparing multiple genomes against previously selected databases.
The result consists of a clustermap and a presence matrix.

As input use GBF or GBK files derived from Prokka.

USAGE:
python3 '''+sys.argv[0]+''' -card -vfdb -bacmet files.gbk\n
Databases:
-bacmet\tAntibacterial Biocide and Metal Resistance Genes Database
-card\tComprehensive Antibiotic Resistance Database
-vfdb\tVirulence Factor Database

Parameters:
-update\tUpdate databases
-u\tSame as -update
-help\tPrint this help
-h\tSame as -help
-keep\tMaintains the protein sequences used, as well as the CDS position files
-k\tSame as -keep

Contact: dlnrodrigues@ufmg.br
		''')
	exit()
else:
	print('''
Hello user!

This script has the function of comparing multiple genomes against previously selected databases.
The result consists of a clustermap and a presence matrix.

Contact: dlnrodrigues@ufmg.br

Let's continue with your analysis.''')

if ("-update" in sys.argv) or ("-u" in sys.argv):
	if "DB" in os.listdir():
		shutil.rmtree("DB")

print("\nChecking your dependences...")
if "Dependences" not in os.listdir():
	os.mkdir("Dependences")
if "diamond" not in os.listdir("Dependences"):
	print("You may not have DIAMOND.\nWe will try to get it.\n")
	atual = os.listdir()
	diamond = wget.download("http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz")
	os.system("tar -xzf diamond-linux64.tar.gz")
	os.system("chmod 755 diamond")
	os.rename("diamond", "Dependences/diamond")
	for i in os.listdir():
		if i not in atual:
			os.remove(i)
	print("\nDownload complete!")

print("\nChecking your databases...")
if "DB" not in os.listdir():
	os.mkdir("DB")
if "bacmet_2.fasta" not in os.listdir("DB"):
	print("\nDownloading BacMet database...")
	bacmet = wget.download("http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta")
	os.rename(bacmet, "DB/bacmet_2.fasta")
	os.system("Dependences/diamond makedb --in DB/bacmet_2.fasta -d DB/bacmet_2 --quiet")
if "vfdb_core.fasta" not in os.listdir("DB"):
	print("\nDownloading VFDB database...")
	vfdb = wget.download("http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz")
	os.rename(vfdb, "vfdb_core.fasta.gz")
	os.system("gunzip -d vfdb_core.fasta.gz")
	os.rename("vfdb_core.fasta", "DB/vfdb_core.fasta")
	os.system("Dependences/diamond makedb --in DB/vfdb_core.fasta -d DB/vfdb_core --quiet")
if "card_protein_homolog_model.fasta" not in os.listdir("DB"):
	atual = os.listdir()
	print("\nDownloading CARD database...")
	card = wget.download("https://card.mcmaster.ca/latest/data")
	os.system("tar -jxvf "+card)
	os.rename("protein_fasta_protein_homolog_model.fasta", "DB/card_protein_homolog_model.fasta")
	for i in os.listdir():
		if i not in atual:
			os.remove(i)
	os.system("Dependences/diamond makedb --in DB/card_protein_homolog_model.fasta -d DB/card_protein_homolog_model --quiet")
'''if "megares_v2.fasta" not in os.listdir("DB"):
	print("\nDownloading MEGARes database...")
	megares = wget.download("https://megares.meglab.org/download/megares_v2.00/megares_full_database_v2.00.fasta")
	os.rename(megares, "DB/megares_v2.fasta")
	os.system("makeblastdb -in DB/megares_v2.fasta -title megares_v2 -dbtype prot -out DB/megares_v2")'''

try:
	print('\nTrying to import \'Pandas\'')
	import pandas as pd
except:
	print("You may not have \'Pandas\'.\nWe will try to install using pip...")
	os.system("pip install pandas --user")
	import pandas as pd
try:
	print('Trying to import \'Seaborn\'')
	import seaborn as sns
except:
	print("You may not have \'Seaborn\'.\nWe will try to install using pip3...")
	os.system("pip install seaborn --user")
	import seaborn as sns

files = []
print("Separating files")
for i in sys.argv:
	if (i.endswith(".gbk")) or (i.endswith(".gbf")):
		files.append(i)

parameters = []
print("Separating parameters")
for i in sys.argv:
	if ("-card" == i) or ("-bacmet" == i) or ("-vfdb" == i):
		parameters.append(i)

############################################Extrair posições#################################################
def extract_positions(a):
	gbk = open(a, 'rt')
	cds = gbk.readlines()
	gbk.close()
	positions = {}
	for i in range(0, len(cds)):
		if "   CDS   " in cds[i]:
			locus_tag = ""
			position = ""
			for j in range(i, len(cds)):
				if "   CDS   " in cds[j]:
					position = cds[j].replace('\n', '')
					position = position.replace("CDS", "")
					position = position.strip()
					if "complement(" in position:
						position = position.replace("complement(", "")
						position = position.replace(")","")
					position = position.strip()
					position = position.replace("..","\t")
				if "/locus_tag=" in cds[j]:
					locus_tag = cds[j].replace("/locus_tag=", "")
					locus_tag = locus_tag.strip()
					locus_tag = locus_tag.replace("\"", "")
					locus_tag = locus_tag.replace("\n", "")
					positions[locus_tag] = position
					break
	return(positions)

print("\nExtracting CDS positions from GenBank files\n")
strains = []
pos = {}
for i in files:
	print("Extracting from file "+i)
	k = extract_positions(i)
	strain = str(i)[::-1].split('/')[0][::-1].replace(".gbk", "").replace(".gbf", "")
	strains.append(strain)
	pos[i] = k

if "Positions_1" not in os.listdir():
	os.mkdir("Positions_1")
else:
	shutil.rmtree("Positions_1")
	os.mkdir("Positions_1")

for i in strains:
	for j in pos.keys():
		if i in j:
			positions = open("Positions_1/"+i+".tab", 'w')
			for k in pos[j].keys():
				positions.write(k+'\t')
				positions.write(pos[j][k])
				positions.write('\n')
############################################Extrair posições#############################################

#################################################Extrair o faa###########################################
def extract_faa(a): #Extract aminoacid sequence from a gbk file
	gbk = open(a, 'rt')
	cds = gbk.readlines()
	gbk.close()
	final = []
	for i in range(0, len(cds)):
		if "   CDS   " in cds[i]:
			locus_tag = ""
			product = ""
			sequence = []
			for j in range(i, len(cds)):
				if "/locus_tag=" in cds[j]:
					locus_tag = cds[j].replace("/locus_tag=", "")
					locus_tag = locus_tag.strip()
					locus_tag = locus_tag.replace("\"", "")
					locus_tag = locus_tag.replace("\n", "")
				elif "/product=" in cds[j]:
					product = cds[j].replace("/product=", "")
					product = product.strip()
					product = product.replace("\"", "")
					product = product.replace("\n", '')
				elif "/translation=" in cds[j]:
					seq = cds[j].replace("/translation=", "")
					seq = seq.replace("\"", "")
					seq = seq.strip()+"\n"
					sequence.append(seq)
					k = j + 1
					while ("   CDS   " not in cds[k]):
						if ("   gene" not in cds[k]) and ("locus_tag" not in cds[k]) and ("gene=" not in cds[k]) and ("assembly_gap" not in cds[k]) and ("_" not in cds[k]):
							seq = cds[k].replace(" ", "")
							seq = seq.replace("\"", "")
							sequence.append(seq)
						if ("\"" in cds[k]) or ("ORIGIN" in cds[k]) or ("assembly_gap" in cds[k]) or ("gene" in cds[k]) or ("_" in cds[k]):
							break
						k = k + 1
					header = ">"+locus_tag+" "+product+"\n"
					final.append(header)
					for l in sequence:
						final.append(l)
					break
	return(final)

if "faa" not in os.listdir():
	os.mkdir("faa")
else:
	shutil.rmtree("faa")
	os.mkdir("faa")

print("\nExtracting protein sequences form GenBank files\n")

for i in files:
	for j in strains:
		if j in i:
			print("Extracting "+j)
			faa = open("faa/"+j+".faa", 'w')
			k = extract_faa(i)
			for l in k:
				faa.write(l)
#################################################Extrair o faa###########################################


for p in parameters:
	dirinicial = os.listdir()
#################################################Alinhar#################################################
	def align(a, b, c):
		#a = input
		#b = title
		#c = result
		cmd = 'Dependences/diamond blastp -q '+a+' -d '+b+' -o '+c+' --quiet -k 1 -e 5e-6 -f 6 qseqid sseqid pident qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore'
		print(cmd)
		os.system(cmd)

	if "Tabular_1" not in os.listdir():
		os.mkdir("Tabular_1")
	else:
		shutil.rmtree("Tabular_1")
		os.mkdir("Tabular_1")

	if "-card" == p:
		b = "DB/card_protein_homolog_model.dmnd"
	elif "-vfdb" == p:
		b = "DB/vfdb_core.dmnd"
	elif "-bacmet" == p:
		b = "DB/bacmet_2.dmnd"

	print("\nStarting alignments\n")
	for i in strains:
		a = "faa/"+i+".faa"
		c = "Tabular_1/"+i+".tab"
		align(a, b, c)
#################################################Alinhar#################################################

#################################################Mineirar#################################################
	def blastmining(a):
		_in = "Tabular_1/"+a
		_out = _in.replace("Tabular_1", "Tabular_2")
		evalue = 5e-06
		identidade = 70
		cobertura = 70
		print(_out)
		try:
			fileO = open(_in, 'rt')
			file = fileO.readlines()
			fileO.close()
			saida = open(_out, 'w')
			for j in file:
				linha = j.split('\t')
				if float(linha[10]) <= evalue:
					if float(linha[2]) >= identidade:
						if float(linha[3]) >= cobertura:
							saida.write(j)
		except:
			y = 0

	if "Tabular_2" not in os.listdir():
		os.mkdir("Tabular_2")
	else:
		shutil.rmtree("Tabular_2")
		os.mkdir("Tabular_2")

	print("\nMining alingments...\n")
	for i in os.listdir("Tabular_1"):
		blastmining(i)
#################################################Mineirar#################################################

#################################################Extrair keys#################################################
	comp = {}
#################BACMET########################
	if "-bacmet" == p:
		bacmetFile = open('DB/bacmet_2.fasta', 'rt')
		bacmet = bacmetFile.readlines()
		bacmetFile.close()
		for i in bacmet:
			if '>' in i:
				ident = i[i.find('>')+1:i.find('|')]
				if ' ' in ident:
					ident = ident.replace(' ', '')
				gene = i.split('|')
				gene = gene[1].replace("\n", "")
				comp[str(ident)] = str(gene)

#################CARD########################
	if "-card" == p:
		cardFile = open('DB/card_protein_homolog_model.fasta', 'rt')
		card = cardFile.readlines()
		cardFile.close()
		for i in card:
			if '>' in i:
				line = i.split("|")
				ident = line[1]
				gene = line[-1][:line[-1].find(" ")].replace("\n", "")
				comp[str(ident)] = str(gene)

#################VFDB########################
	if "-vfdb" == p:
		vfdbFile = open('DB/vfdb_core.fasta', 'rt')
		vfdb = vfdbFile.readlines()
		vfdbFile.close()
		for i in vfdb:
			if '>' in i:
				ident = i[i.find('>')+1:i.find('(')]
				if ' ' in ident:
					ident = ident.replace(' ', '')
				gene = i[i.find(' ')+2:i.find(')',i.find(' '))].replace("\n", "")
				comp[str(ident)] = str(gene)

#################MEGARes########################
	if "-megares" == p:
		megaresFile = open('DB/megares_v2.fasta', 'rt')
		megares = megaresFile.readlines()
		megaresFile.close()
		for i in megares:
			if '>' in i:
				ident = i[i.find('>')+1:i.find('|')]
				if ' ' in ident:
					ident = ident.replace(' ', '')
				gene = i.split('|')
				gene = gene[4].replace("\n", "")
				comp[str(ident)] = str(gene)
#################################################Extrair keys#################################################

#################################################Gerar a matriz#################################################
	if "-card" == p:
		titulo = "matriz_card.csv"
	elif "-vfdb" == p:
		titulo = "matriz_vfdb.csv"
	elif "-bacmet" in sys.argv:
		titulo = "matriz_bacmet.csv"
	elif "-megares" == p:
		titulo = "matriz_megares.csv"

	print("\nGenerating the presence and identity matrix...")

	dicl = {}
	totalgenes = []
	for i in os.listdir("Tabular_2"):
		linhagem = i[:-4]
		i = "Tabular_2/"+i
		file = open(i, 'rt')
		linhas = file.readlines()
		file.close()
		genes = {}
		for j in linhas:
			linha = j[:-1]
			linha = linha.split('\t')
			for k in comp.keys():
				if k in linha[1]:
					gene = comp[k]
					break
			identidade = linha[2]
			genes[gene] = identidade
			if gene not in totalgenes:
				totalgenes.append(gene)
		dicl[str(linhagem)] = genes
	saida = open(titulo, 'w')
	saida.write('Strains;')
	for i in totalgenes:
		saida.write(i+';')
	saida.write('\n')
	for j in dicl.keys(): #j = linhagens
		saida.write(j+';')
		for k in totalgenes: #k = todos os genes possíveis
			if k in dicl[j].keys(): #dicl[j] todos os genes possíveis em todas as linhagens
				saida.write(str(dicl[j][k])+';')
			else:
				saida.write('0;')
		saida.write('\n')
	saida.close()
#################################################Gerar a matriz#################################################

#################################################Gerar positions#################################################
	if "Positions" not in os.listdir():
		os.mkdir("Positions")
	else:
		shutil.rmtree("Positions")
		os.mkdir("Positions")
	print("\nExtracting positions from specific factors...")
	for i in strains:
		pos = {}
		positions = "Positions_1/"+i+".tab"
		positions2 = "Positions/"+i+".tab"
		tab = open(positions, 'rt')
		arq = tab.readlines()
		tab.close()
		for j in arq:
			line = j.split("\t")
			try:
				pos[line[0]] = [line[1], line[2]]
			except:
				u = 0
		file = "Tabular_2/"+i+".tab"
		tab = open(file, 'rt')
		file = tab.readlines()
		tab.close()
		final = {}
		for j in file:
			linha = j.split('\t')
			for k in comp.keys():
				if k in linha[1]:
					gene = comp[k]
					break
			final[linha[0]] = gene, pos[linha[0]]
		out = open(positions2, 'w')
		for j in final.keys():
			if "-card" == p:
				color = "blue\n"
			elif "-vfdb" == p:
				color = "red\n"
			elif "-bacmet" == p:
				color = "green\n"
			elif "-megares" == p:
				color = "yellow\n"
			out.write(final[j][1][0]+'\t'+final[j][1][1][:-1]+'\t'+final[j][0]+'\t'+color)
	out.close()
#################################################Gerar positions#################################################

#################################################Gerar figura###################################################
	if p == "-card":
		out = "clustermap_card.pdf"
		color = "Blues"
	elif p == "-vfdb":
		out = "clustermap_vfdb.pdf"
		color = "Reds"
	elif p == "-bacmet":
		out = "clustermap_bacmet.pdf"
		color = "Greens"

	df = pd.read_csv(titulo, sep=';')
	df = df.set_index('Strains')
	headers = list(df.columns.values)
	lines = list(df.index.values)
	for i in headers:
		if "Unnamed:" in i:
			df = df.drop(columns=[i])
	x = math.ceil(len(headers)*0.65)
	y = math.ceil(len(lines)*0.65)
	if x >= 100:
		x = 100
	if y >= 100:
		y = 100
	print("\nPlotting final clustermap...")
	p2 = sns.clustermap(df, figsize=(x,y), cmap=color)
	p2.savefig(out, format='pdf', dpi=200, bbox_inches="tight")
#################################################Gerar figura###################################################

#################################################Organizando#################################################
	print("\nGrouping the results...\n")
	diratual = os.listdir()
	if "Results_1" not in diratual:
		os.mkdir("Results_1")
		atual = "Results_1"
	else:
		ind = 2
		atual = "Results_"+str(ind)
		while "Results_"+str(ind) in diratual:
			ind = ind + 1
			atual = "Results_"+str(ind)
		os.mkdir(atual)

	for i in diratual:
		if i not in dirinicial:
			shutil.move(i, atual)
#################################################Organizando#################################################

#################################################Removendo arquivos intermediarios#################################################
if ("-keep" not in sys.argv) and ("-k" not in sys.argv):
	shutil.rmtree("Positions_1")
	shutil.rmtree("faa")
#################################################Removendo arquivos intermediarios#################################################

print('''That's all, folks!
Thank you for using this script.
Do not forget to quote the databases used.\n''')
if "-bacmet" in sys.argv:
	print("BacMet\thttps://doi.org/10.1093/nar/gkt1252\t2014")
if "-card" in sys.argv:
	print("CARD\thttps://doi.org/10.1093/nar/gkz935\t2020")
if "-vfdb" in sys.argv:
	print("VFDB\thttps://doi.org/10.1093/nar/gky1080\t2019")
print('')
