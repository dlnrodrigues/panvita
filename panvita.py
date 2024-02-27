#!/usr/bin/env python3

import sys
import re
import os
import shutil
import math
import time
import gzip
import string
import random
from datetime import datetime
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
		
version = ("1.1.8")

if ("-v" in sys.argv) or ("-version" in sys.argv):
	print("-----------------------------------------------")
	print("PanViTa - Pan Virulence and resisTance Analysis")
	print("https://doi.org/10.3389/fbinf.2023.1070406")
	print("version", version)
	print("-----------------------------------------------")
	exit()

if (("-card" not in sys.argv) and ("-bacmet" not in sys.argv) and ("-vfdb" not in sys.argv) and ("-u" not in sys.argv) and ("-update" not in sys.argv) and ("-g" not in sys.argv) and ("-a" not in sys.argv) and ("-m" not in sys.argv) and ("-b" not in sys.argv)) or ("-h" in sys.argv) or ("-help" in sys.argv):
	print('''
Hello user!

PanViTa has the function of comparing multiple genomes against previously selected databases.
The result consists of a clustermap and a presence matrix.

As input use GBF or GBK files derived from Prokka or available on NCBI.
WARNING! Files from NCBI MUST have .gbf or .gbff extension.

USAGE:
python3 '''+sys.argv[0]+''' -card -vfdb -bacmet files.gbk\n
Databases:
-bacmet\tAntibacterial Biocide and Metal Resistance Genes Database
-card\tComprehensive Antibiotic Resistance Database
-vfdb\tVirulence Factor Database

Parameters:
-update\tUpdate databases and dependences
-u\tSame as -update
-help\tPrint this help
-h\tSame as -help
-v\tPrint version and exit
-keep\tMaintains the protein sequences used, as well as the CDS position files
-k\tSame as -keep
-i\tMinimum identity to infer presence (default = 70)
-c\tMinimum coverage to infer presence (default = 70)
-d\tForce to use DIAMOND from system
-pdf\tFigures will be saved as PDF (default)
-png\tFigures will be saved as PNG (WARNING! High memory consumption)
-g\tDownload the genomes fasta files (require CSV table from NCBI)
-a\tDownload and annote the genomes using PROKKA pipeline (require CSV table from NCBI)
-b\tDownload the genome GenBank files (require CSV table from NCBI)
-s\tKeep the locus_tag as same as the strain (require -b)
-m\tGet the metadata from BioSample IDs (require CSV table from NCBI)

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

##################################Functions start################################
def checkDP():
	print("\nChecking your dependences...")
	home = os.path.expanduser('~')
	config_file = home+"/.panvita.dp.paths"
	if os.path.exists(config_file) == True:
		with open(config_file, "rt") as file:
			cfile = file.readlines()
			dppath = cfile[0].strip()
			file.close()
	else:
		dppath = os.getcwd()+"/Dependences/"
		with open(config_file, "w") as file:
			file.write(dppath)
			file.close()
	if "Dependences" not in os.listdir(str(dppath.replace("Dependences/",""))):
		os.mkdir(dppath)
	if "diamond" not in os.listdir(dppath):
		print("You may not have DIAMOND.\nWe will try to get it.\n")
		atual = os.listdir()
		diamond = wget.download("http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz")
		os.system("tar -xzf diamond-linux64.tar.gz")
		os.system("chmod 755 diamond")
		os.rename("diamond", dppath+"diamond")
		for i in os.listdir():
			if i not in atual:
				os.remove(i)
		print("\nDIAMOND - Download complete!")
	return(dppath)

def checkDB():
	print("\nChecking your databases...")
	home = os.path.expanduser('~')
	config_file = home+"/.panvita.db.paths"	
	if os.path.exists(config_file) == True:
		with open(config_file, "rt") as file:
			cfile = file.readlines()
			dbpath = cfile[0].strip()
			file.close()
	else:
		dbpath = os.getcwd()+"/DB/"
		with open(config_file, "w") as file:
			file.write(dbpath)
			file.close()
	if "DB" not in os.listdir(str(dbpath.replace("DB/",""))):
		os.mkdir(dbpath)

	if ("bacmet_2.fasta" not in os.listdir(dbpath)) or ("bacmet_2.txt" not in os.listdir(dbpath)):
		print("\nDownloading BacMet database...")
		bacmet = wget.download("http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta")
		os.rename(bacmet, dbpath+"bacmet_2.fasta")
		print("")
		os.system(diamond_exe+" makedb --in "+dbpath+"bacmet_2.fasta -d "+dbpath+"bacmet_2 --quiet")
		print("\nDownloading BacMet annotation file...")
		bacmet_an = wget.download("http://bacmet.biomedicine.gu.se/download/BacMet2_EXP.753.mapping.txt")
		os.rename(bacmet_an, dbpath+"bacmet_2.txt")
	if "vfdb_core.fasta" not in os.listdir(dbpath):
		print("\nDownloading VFDB database...")
		vfdb = wget.download("http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz")
		os.rename(vfdb, "vfdb_core.fasta.gz")
		os.system("gunzip -d vfdb_core.fasta.gz")
		os.rename("vfdb_core.fasta", dbpath+"vfdb_core.fasta")
		print("")
		os.system(diamond_exe+" makedb --in "+dbpath+"vfdb_core.fasta -d "+dbpath+"vfdb_core --quiet")
	if "card_protein_homolog_model.fasta" not in os.listdir(dbpath):
		atual = os.listdir()
		print("\nDownloading CARD database...")
		card = wget.download("https://card.mcmaster.ca/latest/data")
		os.system("tar -jxvf "+card)
		os.rename("protein_fasta_protein_homolog_model.fasta", dbpath+"card_protein_homolog_model.fasta")
		os.rename("aro_index.tsv", dbpath+"aro_index.tsv")
		for i in os.listdir():
			if i not in atual:
				os.remove(i)
		print("")
		os.system(diamond_exe+" makedb --in "+dbpath+"card_protein_homolog_model.fasta -d "+dbpath+"card_protein_homolog_model --quiet")
	if "latlon.csv" not in os.listdir(dbpath):
		print("\nDownloading coordinates keys file...")
		latlon = wget.download("https://raw.githubusercontent.com/dlnrodrigues/panvita/dlnrodrigues-Supplementary/latlon.csv")
		os.rename(latlon, dbpath+"latlon.csv")
		print("")
	'''if "megares_v2.fasta" not in os.listdir(dbpath):
		print("\nDownloading MEGARes database...")
		megares = wget.download("https://megares.meglab.org/download/megares_v2.00/megares_full_database_v2.00.fasta")
		os.rename(megares, dbpath+"megares_v2.fasta")
		os.system("makeblastdb -in "+dbpath+"megares_v2.fasta -title megares_v2 -dbtype prot -out "+dbpath+"megares_v2")'''
	return(dbpath)
def getMeta(a):
	mlst = shutil.which("mlst")
	out = open(a, "w")
	if type(mlst) != str:
		out.write("Species\tStrain\tHost\tHost Disease\tGeografic Localization\tIsolation Source\tGenomic Size\tGC%\tRelease Date\tBioSample\n")
	else:
		out.write("Species\tStrain\tHost\tHost Disease\tGeografic Localization\tIsolation Source\tGenomic Size\tGC%\tRelease Date\tST\tBioSample\n")
	out.close()
	for i in dic2.keys():
		out = open(a, "a")
		spex = dic2[i][0].split(" ")
		spex = spex[0]+" "+spex[1]
		out.write(str(spex)+"\t")
		temp = re.sub("((?![\.A-z0-9_-]).)", "_", str(dic2[i][3]))
		temp3 = re.sub("((?![\.A-z0-9_-]).)", "_", str(dic2[i][3]))
		out.write(temp+"\t")
		ID = "https://www.ncbi.nlm.nih.gov/biosample/"+str(i)
		print(ID)
		file = wget.download(ID)
		print("")
		html = open(file, 'rt')
		txt = html.readlines()
		html.close()
		x = ""
		for j in txt:
			temp = j.strip().replace("\n", "")
			x = x + temp
		if x.find("<th>host</th><td>") != -1:
			host = (x[(x.find("<th>host</th><td>")+17):(x.find("<", (x.find("<th>host</th><td>")+17)))]).strip().capitalize()
		else:
			host = "NA"
		out.write(host+"\t")
		if x.find("host disease</th><td>") != -1:
			disease = (x[(x.find("host disease</th><td>")+21):(x.find("<", (x.find("host disease</th><td>")+21)))]).strip().capitalize()
		else:
			disease = "NA"
		out.write(disease+"\t")
		if x.find("geo_loc_name=") != -1:
			geo = (x[(x.find("geo_loc_name=")+13):(x.find("&", (x.find("geo_loc_name=")+13)))]).strip()
		else:
			geo = "NA"
		out.write(geo+"\t")
		if x.find(">isolation source</th><td>") != -1:
			source = (x[(x.find(">isolation source</th><td>")+26):(x.find("<", (x.find(">isolation source</th><td>")+26)))]).strip().capitalize()
		else:
			source = "NA"
		out.write(source+"\t")
		os.remove(file)
		out.write(str(dic2[i][4])+"\t")
		out.write(str(dic2[i][5])+"\t")
		out.write(str(dic2[i][6][:dic2[i][6].find("T")])+"\t")
		if type(mlst) == str:
			if "-b" in sys.argv:
				try:
					print(temp3)
					os.system("mlst --quiet "+temp3+".gbf > .temp.temp")
					sttemp = open(".temp.temp", "rt")
					st = sttemp.readlines()[0].split("\t")[2]
					sttemp.close()
					os.remove(".temp.temp")
					out.write(str(st)+"\t")
				except:
					out.write("NA\t")
		out.write(str(i)+"\n")
		out.close()

	meta = pd.read_csv("meta_data.tsv", sep="\t")
	countries = meta["Geografic Localization"].tolist()
	latlon = dbpath+"latlon.csv"
	countries_keys = pd.read_csv(latlon, sep=";", index_col="homecontinent")

	unique = []
	k = []
	for i in countries:
		if type(i) == str:
			if ":" not in i:
				string = i
				k.append(string)
			else:
				string = i.split(":")[0]
				k.append(string)
			if string not in unique:
				unique.append(string)

	cont = []
	for i in unique:
		cont.append(k.count(i))

	data = {}
	homelat = []
	homelon = []
	cont = []
	temp_unique = []
	for i in unique:
		if i in countries_keys.index.values.tolist():
			homelat.append(countries_keys["homelat"][i])
			homelon.append(countries_keys["homelon"][i])
			cont.append(k.count(i))
			temp_unique.append(i)
		else:
			print(f"Unidentified country: {i}")
	data = {"homecontinent": temp_unique,
			"homelat": homelat,
			"homelon": homelon,
			"n": cont}
	data = pd.DataFrame.from_dict(data)
	data.to_csv("meta_data_countries_count.csv")

	try:
		plt.figure(figsize=(20, 15))
		plt.rcParams["figure.figsize"]=20,15;

		m=Basemap(llcrnrlon=-180, llcrnrlat=-65, urcrnrlon=180, urcrnrlat=80)
		m.drawmapboundary(fill_color='#A6CAE0', linewidth=0)
		m.fillcontinents(color='green', alpha=0.3)
		m.drawcoastlines(linewidth=0.1, color="white")

		data['labels_enc'] = pd.factorize(data['homecontinent'])[0]

		m.scatter(
			x=data['homelon'], 
			y=data['homelat'], 
			s=data['n']*10, 
			alpha=0.4, 
			c=data['labels_enc'], 
			cmap="plasma")
		 
		#plt.text(-175, -62,'Isolates geographical localization', ha='left', va='bottom', size=, color='#555555' )

		plt.savefig('meta_data_maps.png', dpi=300, bbox_inches="tight")
	except:
		erro_string = "ERROR: It was'nt possible to plot the final map.\nPlease check the log file.\n"
		erro.append(erro_string)
		print(erro_string)

def getNCBI_GBF():
	gbff = []
	attempts = []
	removal = []
	all_strains = []
	for i in dic.keys():
		if type(dic[i][1]) == str:
			ftp = dic[i][1]
			file = "/"+ftp[:: -1].split("/")[0][:: -1]+"_genomic.gbff.gz"
			ftp = ftp + file
		elif type(dic[i][2]) == str:
			ftp = dic[i][2]
			file = "/"+ftp[:: -1].split("/")[0][:: -1]+"_genomic.gbff.gz"
			ftp = ftp + file
		else:
			erro_string = f"ERROR: It was'nt possible to download the file {str(i)}.\nIt dosen't have a FTP accession.\nPlease check the log file.\n"
			erro.append(erro_string)
			print(erro_string)
			continue
		genus = dic[i][0].split(" ")[0]
		species = dic[i][0].split(" ")[1]
		strain = re.sub("((?![\.A-z0-9_-]).)", "_", str(dic[i][3]))
		if "-s" in sys.argv:
			ltag = strain
		else:
			ltag = genus[0]+species+"_"+strain
		attempts.append((ftp, species, genus, strain, ltag))
		indic = 0
	while len(attempts) != 0:
		try:
			print(f"Strain {attempts[0][3]}: attempt {indic + 1}\n{attempts[0][0]}")
			file = wget.download(attempts[0][0])
			print("\n")
			newfile = re.sub("((?![\.A-z0-9_-]).)", "_", str(file))
			os.rename(file, newfile)
			removal.append(newfile)
			file = newfile
			os.system(f"gunzip {file}")
			file = file.replace(".gz", "")
			dic[i] = (dic[i][0], file)
			try:
				if attempts[0][4] not in all_strains:
					temp_string = attempts[0][4]
					os.rename(file, attempts[0][4]+".gbf")
					all_strains.append(attempts[0][4])
				else:
					temp_string = f"{attempts[0][4]}_dup_{''.join(random.choices(string.ascii_uppercase + string.digits, k=5))}"
					os.rename(file, f"{temp_string}.gbf")
			except:
				time.sleep(3)
				if attempts[0][4] not in all_strains:
					temp_string = attempts[0][4]
					os.rename(file, attempts[0][4]+".gbf")
					all_strains.append(attempts[0][4])
				else:
					temp_string = f"{attempts[0][4]}_dup_{''.join(random.choices(string.ascii_uppercase + string.digits, k=5))}"
					os.rename(file, f"{temp_string}.gbf")
			gbff.append(f"./{temp_string}.gbf")
			attempts.pop(0)
			indic = 0
		except:
			indic = indic + 1
			if indic > 4:
				attempts.pop(0)
				erro_string = f"ERROR: It was'nt possible to download the GenBank file for the strain {attempts[0][3]} even after 5 attempts.\nPlease check the internet conection, the log and input files.\n"
				erro.append(erro_string)
				print(erro_string)
				indic = 0
				continue
			else:
				continue
	for i in os.listdir("."):
		if i.endswith(".tmp"):
			os.remove(i)
		elif i in removal:
			os.remove(i)
	return(gbff)

def getNCBI_FNA():
	prokka = shutil.which("prokka")
	pkgbf = []
	pk = []
	all_strains = []
	for i in dic3.keys():
		attempts = 1
		while True:
			genus = dic3[i][0].split(" ")[0]
			species = dic3[i][0].split(" ")[1]
			strain = re.sub("((?![\.A-z0-9_-]).)", "_", str(dic3[i][3]))
			if strain not in all_strains:
				all_strains.append(strain)
			else:
				strain = f"{strain}_dup_{''.join(random.choices(string.ascii_uppercase + string.digits, k=5))}"
			if "-s" in sys.argv:
				new_file = strain
			else:
				new_file = genus[0]+species+"_"+strain
			if f"{new_file}.fna" in os.listdir():
				file = new_file+".fna"
				ltag = genus[0]+species[0]+"_"+strain
				temp = "./"+ltag+"/"+strain+".gbf"
				pkgbf.append(temp)
				cmd = ("prokka --addgenes --force --species "+species+" --genus "+genus+" --strain "+strain+" "+file+" --prefix "+ltag+" --outdir "+ltag+" --locustag "+ltag)
				print(f"Skipping file {new_file}.fna. File already exists.\n")
				if ltag not in os.listdir():
					pk.append(cmd)
				print(pk)
				break
			try:
				if type(dic3[i][1]) == str:
					ftp = dic3[i][1]
					file = "/"+ftp[:: -1].split("/")[0][:: -1]+"_genomic.fna.gz"
					ftp = ftp + file
					print("Strain",dic3[i][3],f"Attempt: {attempts}","\n",ftp)
				elif type(dic3[i][2]) == str:
					ftp = dic3[i][2]
					file = "/"+ftp[:: -1].split("/")[0][:: -1]+"_genomic.fna.gz"
					ftp = ftp + file
					print("Strain",dic3[i][3],f"Attempt: {attempts}","\n",ftp)
				else:
					erro_string = "ERROR: It was'nt possible to download the file "+str(i)+".\nPlease check the log file.\n"
					erro.append(erro_string)
					print(erro_string)
					break
				file = wget.download(ftp)
				print("\n")
				os.system("gunzip "+str(file))
				file = file.replace(".gz", "")
				dic3[i] = (dic3[i][0], file)
				try:
					file = os.rename(file, new_file+".fna")
				except:
					time.sleep(1)
					file = os.rename(file, new_file+".fna")
				file = new_file+".fna"
				ltag = genus[0]+species[0]+"_"+strain
				temp = "./"+ltag+"/"+strain+".gbf"
				pkgbf.append(temp)
				cmd = ("prokka --addgenes --force --species "+species+" --genus "+genus+" --strain "+strain+" "+file+" --prefix "+ltag+" --outdir "+ltag+" --locustag "+ltag)
				if ltag not in os.listdir():
					pk.append(cmd)
				break
			except:
				if attempts < 5:
					attempts = attempts + 1
				else:
					erro_string = f"ERROR: It was'nt possible to download the Fasta file for the strain {dic3[i][3]} even after 5 attempts.\nPlease check the internet conection, the log and input files.\n"
					print(erro_string)
					erro.append(erro_string)
					break
	uscript = open("PROKKA.sh","w")
	for i in pk:
		uscript.write(i+"\n")
	uscript.close()
	if "-a" in sys.argv:
		if type(prokka) == str:
			for i in pk:
				print(i)
				os.system(i)
		else:
			erro_string = "Sorry but we didn't find PROKKA in your computer.\nBe sure that the installation was performed well.\nThe annotation will not occur.\nIf you install PROKKA some day, you can use a script we made specially for you!\n"
			print(erro_string)
			erro.append(erro_string)
			pkgbf = [""]
	return(pkgbf)

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
					if cds[j].count("\"") == 2:
						seq = cds[j].replace("/translation=", "")
						seq = seq.replace("\"", "")
						seq = seq.strip()+"\n"
						header = ">"+locus_tag+" "+product+"\n"
						final.append(header)
						final.append(seq)
						break
					else:
						seq = cds[j].replace("/translation=", "")
						seq = seq.replace("\"", "")
						seq = seq.strip()+"\n"
						sequence.append(seq)
						k = j + 1
						if "\"" in cds[k]:
							seq = cds[k].replace(" ", "")
							seq = seq.replace("\"", "")
							sequence.append(seq)
						else:							
							while ("\"" not in cds[k]):
								seq = cds[k].replace(" ", "")
								seq = seq.replace("\"", "")
								sequence.append(seq)
								k = k + 1
							seq = cds[k].replace(" ", "")
							seq = seq.replace("\"", "")
							sequence.append(seq)						
						header = ">"+locus_tag+" "+product+"\n"
						final.append(header)
						for l in sequence:
							final.append(l)
						break
	return(final)

def extract_positions(a):
	gbk = open(a, 'rt')
	cds = gbk.readlines()
	gbk.close()
	positions = {}
	lenght = 0
	totalcds = 0
	for i in range(0, len(cds)):
		if ("   CDS   " in cds[i]) and ("   ::" not in cds[i]):
			locus_tag = ""
			position = ""
			for j in range(i, len(cds)):
				if "   CDS   " in cds[j]:
					totalcds = totalcds + 1
					position = cds[j].replace('\n', '')
					position = position.replace("CDS", "")
					position = position.strip()
					if (">" in position) or ("<" in position):
						position = position.replace(">", "").replace("<","")
						#break
					if "complement(join(" in position:
						position = position.replace("complement(join(", "")
						position = position.replace(")","").strip()
						position = position.split("..")
						if "," in position[0]:
							temp = position[0].split(",")
							position[0] = temp[0]
						if totalcds == 1:
							try:
								position = [int(position[0])+lenght, int(position[1])+lenght]
							except:
								position = [int(position[0])+lenght, int(position[2])+lenght]
					elif "complement(" in position:
						position = position.replace("complement(", "")
						position = position.replace(")","").strip()
						position = position.split("..")
						position = [int(position[0])+lenght, int(position[1])+lenght]
					elif "join(" in position:
						position = position.replace("join(", "")
						position = position.replace(")","").strip()
						position = position.split("..")
						if "," in position[0]:
							temp = position[0].split(",")
							position[0] = temp[0]
							print(position)
						if totalcds == 1:
							try:
								position = [int(position[0])+lenght, int(position[1])+lenght]
							except:
								position = [int(position[0])+lenght, int(position[2])+lenght]
					else:
						position = position.replace("<","").replace(">","").strip()
						position = position.split("..")
						position = [int(position[0])+lenght, int(position[1])+lenght]
				if "/locus_tag=" in cds[j]:
					locus_tag = cds[j].replace("/locus_tag=", "")
					locus_tag = locus_tag.strip()
					locus_tag = locus_tag.replace("\"", "")
					locus_tag = locus_tag.replace("\n", "")
					positions[locus_tag] = str(position[0])+"\t"+str(position[1])
					break
		if "CONTIG " in cds[i]:
			lenght = lenght + int(cds[i].strip().replace(")","").split("..")[-1])
	return(positions)

def align(a, b, c):
	#a = input
	#b = title
	#c = result
	cmd = diamond_exe+' blastp -q '+a+' -d '+b+' -o '+c+' --quiet -k 1 -e 5e-6 -f 6 qseqid sseqid pident qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore'
	print(cmd)
	os.system(cmd)

def blastmining(a):
	_in = "Tabular_1/"+a
	_out = _in.replace("Tabular_1", "Tabular_2")
	evalue = 5e-06
	identidade = 70
	cobertura = 70
	if "-i" in sys.argv:
		try:
			identidade = float(sys.argv[sys.argv.index("-i")+1])
		except:
			cobertura = 70
	if "-c" in sys.argv:
		try:
			cobertura = float(sys.argv[sys.argv.index("-c")+1])
		except:
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

def writeERR(erro_list):
	erro_file = "panvita_error_"+datetime.now().strftime("%d-%m-%Y_%H-%M-%S")+".err"
	with open(erro_file, "w") as erro_out:
		string = "PanViTa "+str(datetime.now().strftime("%d-%m-%Y_%H-%M-%S"))+"\n"
		erro_out.write(string)
		erro_out.write("The following lines contains the errors from the last runing.\nPlease, check them carefully.\n")
		for i in erro_list:
			erro_out.write(i)
		erro_out.write("\nIf the lines appear empty, no error has been computed for the current run.\n")
		erro_out.close()
	return(erro_file)
##################################Functions end################################

erro = []

if ("-update" in sys.argv) or ("-u" in sys.argv):
	home = os.path.expanduser('~')
	config_file = home+"/.panvita.dp.paths"
	if ".panvita.dp.paths" in os.listdir(home):
		with open(config_file, "rt") as file:
			rm_path = file.readlines()[0]
			file.close()
		os.remove(config_file)
		shutil.rmtree(rm_path)
		dppath = checkDP()
	else:
		dppath = checkDP()
	diamond_exe = dppath+"diamond"
	if ("-d" in sys.argv) or ("-diamond" in sys.argv):
		if type(shutil.which("diamond-aligner")) == str:
			diamond_exe = shutil.which("diamond-aligner")
		elif type(shutil.which("diamond")) == str:
			diamond_exe = shutil.which("diamond")
		else:
			erro_string = "\nWe could'nt locate DIAMOND on your system.\nWe'll try to use the default option.\nPlease verify the alignment outputs.\n"
			print(erro_string)
			erro.append(erro_string)

	config_file = home+"/.panvita.db.paths"
	if ".panvita.db.paths" in os.listdir(home):
		with open(config_file, "rt") as file:
			rm_path = file.readlines()[0]
			file.close()
		os.remove(config_file)
		shutil.rmtree(rm_path)
		dbpath = checkDB()
	else:
		dbpath = checkDB()
	if ("-card" not in sys.argv) and ("-bacmet" not in sys.argv) and ("-vfdb" not in sys.argv):
		print('''
Your databases and dependences were updated. ^^

Contact: dlnrodrigues@ufmg.br
''')
		exit()

dppath = checkDP()
diamond_exe = dppath+"diamond"
if ("-d" in sys.argv) or ("-diamond" in sys.argv):
	if type(shutil.which("diamond-aligner")) == str:
		diamond_exe = shutil.which("diamond-aligner")
	elif type(shutil.which("diamond")) == str:
		diamond_exe = shutil.which("diamond")
	else:
		erro_string = "\nWe could'nt locate DIAMOND on your system.\nWe'll try to use the default option.\nPlease verify the alignment outputs.\n"
		print(erro_string)
		erro.append(erro_string)


####################################################
dbpath = checkDB()
print("\nDependences path:\n",dppath)
print("Selected DIAMOND:\n",diamond_exe)
print("Databases path:\n",dbpath)
#######################################################################################################################################################

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
try:
	print('Trying to import \'Matplotlib\'')
	import matplotlib.pyplot as plt
except:
	print("You may not have \'Matplotlib\'.\nWe will try to install using pip3...")
	os.system("python -m pip install -U matplotlib")
	import matplotlib.pyplot as plt
try:
	print('Trying to import \'SciPy\'')
	import scipy
except:
	print("You may not have \'SciPy\'.\nWe will try to install using pip...")
	os.system("pip install scipy --user")
	import scipy
try:
	print('Trying to import \'Basemap\'\n')
	from mpl_toolkits.basemap import Basemap
except:
	print("You may not have \'Basemap\'.\nWe will try to install using pip...")
	os.system("pip install basemap --user")
	from mpl_toolkits.basemap import Basemap
#############################################NEW################################################################
if ("-a" in sys.argv) or ("-b" in sys.argv) or ("-g" in sys.argv) or ("-m" in sys.argv):
	for i in sys.argv:
		if i.endswith(".csv"):
			table = i
	df = pd.read_csv(table, sep = ',')
	species = df["#Organism Name"].tolist()
	strains = df["Strain"].tolist()
	biosample = df["BioSample"].tolist()
	size = df["Size(Mb)"].tolist()
	GC = df["GC%"].tolist()
	refseq = df["RefSeq FTP"].tolist()
	genbank = df["GenBank FTP"].tolist()
	data = df["Release Date"].tolist()
	dic = {}
	dic2 = {}
	dic3 = {}
	ind = 0
	while ind in range(0, len(strains)):
		if "-b" in sys.argv:
			dic[biosample[ind]] = (species[ind], refseq[ind], genbank[ind], strains[ind], size[ind], GC[ind], data[ind])
		if "-m" in sys.argv:
			dic2[biosample[ind]] = (species[ind], refseq[ind], genbank[ind], strains[ind], size[ind], GC[ind], data[ind])
		if ("-a" in sys.argv) or ("-g" in sys.argv):
			dic3[biosample[ind]] = (species[ind], refseq[ind], genbank[ind], strains[ind], size[ind], GC[ind], data[ind])
		ind = ind + 1
##############################################NEW############################################################
if "-b" in sys.argv:
	gbff = getNCBI_GBF()

if ("-a" in sys.argv) or ("-g" in sys.argv):
	pkgbf = getNCBI_FNA()

if "-m" in sys.argv:
	getMeta("meta_data.tsv")
##############################################NEW############################################################

files = []
print("Separating files")
for i in sys.argv:
	if (i.endswith(".gbk")) or (i.endswith(".gbf")) or (i.endswith(".gbff")):
		files.append(i)
for i in sys.argv:
	if (i.endswith(".csv")):
		prokaryotes = i
		break
if ("-a" in sys.argv) and ("-b" not in sys.argv):
	if pkgbf[0] != "":
		for i in pkgbf:
			files.append(i)
elif "-b" in sys.argv:
	for i in gbff:
		files.append(i)

parameters = []
print("Separating parameters")
for i in sys.argv:
	if ("-card" == i) or ("-bacmet" == i) or ("-vfdb" == i):
		parameters.append(i)

if len(parameters) == 0:
	if "-b" in sys.argv:
		pasta = "GenBank_"+datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
		os.mkdir(pasta)
		for i in gbff:
			try:
				shutil.move(i, pasta)
			except:
				erro_string = "It was not possible to move the file "+str(i)+" to the final directory.\nPlease, chack the output path.\n"
				print(erro_string)
				erro.append(erro_string)
	if len(erro) > 0:
		final_erro = writeERR(erro)
		print("\nNumber of errors reported: "+str(len(erro)))
		print("Please check the "+str(final_erro)+" file.\n")
	else:
		print("\nNo error reported!")
	print("\nThat's all folks...\nThank you so much for using this program!\n")
	exit()
############################################Extrair posições#################################################

print("\nExtracting CDS positions from GenBank files\n")
strains = []
pos = {}
tempfiles = []
for i in files:
	try:
		k = extract_positions(i)
	except:
		try:
			f = i.replace(".gbff", ".gbk")
			k = extract_positions(f)
		except:
			try:
				f = i.replace(".gbf", ".gbk")
				k = extract_positions(f)
			except:
				if os.path.exists(i) == True:
					erro_string = "\n**WARNING**\nIt was not possible to handle the file "+str(i)+"...\nIt will be skiped.\nPlease verify the input format.\n"
					print(erro_string)
					erro.append(erro_string)
					continue
				else:
					erro_string = "\n**WARNING**\nIt was not possible to handle the file "+str(i)+"...\nIt will be skiped.\Please verify the absolute path of the files.\n"
					print(erro_string)
					erro.append(erro_string)
					continue
	if len(k) < 10:
		erro_string = "\n**WARNING**\nThe file "+str(i)+" seems to be empty...\nIt will be skiped.\nPlease verify the input format.\n"
		print(erro_string)
		erro.append(erro_string)
		continue		
	print("The positions of the "+i+" file have been extracted")
	tempfiles.append(i)
	strain = str(i)[::-1].split('/')[0][::-1].replace(".gbk", "").replace(".gbff", "").replace(".gbf", "")
	strains.append(strain)
	files = tempfiles
	pos[i] = k
del tempfiles

if "Positions_1" not in os.listdir():
	os.mkdir("Positions_1")
else:
	shutil.rmtree("Positions_1")
	os.mkdir("Positions_1")

for i in strains:
	for j in pos.keys():
		tempName = j[::-1].split("/")[0][::-1].replace(".gbk","").replace(".gbff", "").replace(".gbf", "")
		if i == tempName:
			positions = open("Positions_1/"+i+".tab", 'w')
			for k in pos[j].keys():
				positions.write(k+'\t')
				positions.write(pos[j][k])
				positions.write('\n')
############################################Extrair posições#############################################

#################################################Extrair o faa###########################################
if "faa" not in os.listdir():
	os.mkdir("faa")
else:
	shutil.rmtree("faa")
	os.mkdir("faa")

print("\nExtracting protein sequences form GenBank files\n")

for i in files:
	for j in strains:
		tempName = i[::-1].split("/")[0][::-1].replace(".gbk","").replace(".gbff", "").replace(".gbf", "")
		if j == tempName:
			print("Extracting "+j)
			faa = open("faa/"+j+".faa", 'w')
			try:
				k = extract_faa(i)
			except:
				try:
					temp = i.replace(".gbff", ".gbk")
					k = extract_faa(temp)
				except:
					temp = i.replace(".gbf", ".gbk")
					k = extract_faa(temp)
			for l in k:
				faa.write(l)
#################################################Extrair o faa###########################################

for p in parameters:
	outputs = []
	dirinicial = os.listdir()
#################################################Alinhar#################################################
	if "Tabular_1" not in os.listdir():
		os.mkdir("Tabular_1")
		outputs.append("Tabular_1")
	else:
		shutil.rmtree("Tabular_1")
		os.mkdir("Tabular_1")
		outputs.append("Tabular_1")

	if "-card" == p:
		b = dbpath+"card_protein_homolog_model.dmnd"
	elif "-vfdb" == p:
		b = dbpath+"vfdb_core.dmnd"
	elif "-bacmet" == p:
		b = dbpath+"bacmet_2.dmnd"

	print("\nStarting alignments\n")
	for i in strains:
		a = "faa/"+i+".faa"
		c = "Tabular_1/"+i+".tab"
		align(a, b, c)
#################################################Alinhar#################################################

#################################################Mineirar#################################################
	if "Tabular_2" not in os.listdir():
		os.mkdir("Tabular_2")
		outputs.append("Tabular_2")
	else:
		shutil.rmtree("Tabular_2")
		os.mkdir("Tabular_2")
		outputs.append("Tabular_2")

	print("\nMining alingments...\n")
	for i in os.listdir("Tabular_1"):
		blastmining(i)
#################################################Mineirar#################################################

#################################################Extrair keys#################################################
	comp = {}
#################BACMET########################
	if "-bacmet" == p:
		bacmetFile = open(dbpath+'bacmet_2.fasta', 'rt')
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
		cardFile = open(dbpath+'card_protein_homolog_model.fasta', 'rt')
		card = cardFile.readlines()
		cardFile.close()
		for i in card:
			if '>' in i:
				line = i.split("|")
				ident = line[1]
				temp2 = line[-1].split("[")
				temp2 = temp2[0].split(" ")
				if len(temp2) > 2:
					gene = str(temp2[2])
				else:
					gene = line[-1][:line[-1].find(" ")].replace("\n", "")
				comp[str(ident)] = str(gene)

#################VFDB########################
	if "-vfdb" == p:
		genes_comp = {}
		vfdbFile = open(dbpath+'vfdb_core.fasta', 'rt')
		vfdb = vfdbFile.readlines()
		vfdbFile.close()
		for i in vfdb:
			if '>' in i:
				mech = re.findall("(?<=\)\s-\s)[A-z\/\-\s]*(?=\s\()", i, flags=0)
				if len(mech) == 1:
					mech = mech[0]
				else:
					print(i)
				ident = i[i.find('>')+1:i.find('(')]
				if ' ' in ident:
					ident = ident.replace(' ', '')
				gene = i[i.find(' ')+2:i.find(')',i.find(' '))].replace("\n", "")
				comp[str(ident)] = str(gene)
				genes_comp[str(gene)] = str(mech)

#################MEGARes########################
	'''if "-megares" == p:
		megaresFile = open(dbpath+'megares_v2.fasta', 'rt')
		megares = megaresFile.readlines()
		megaresFile.close()
		for i in megares:
			if '>' in i:
				ident = i[i.find('>')+1:i.find('|')]
				if ' ' in ident:
					ident = ident.replace(' ', '')
				gene = i.split('|')
				gene = gene[4].replace("\n", "")
				comp[str(ident)] = str(gene)'''
#################################################Extrair keys#################################################

#################################################Gerar a matriz#################################################
	if "-card" == p:
		titulo = "matriz_card.csv"
		outputs.append(titulo)
	elif "-vfdb" == p:
		titulo = "matriz_vfdb.csv"
		outputs.append(titulo)
	elif "-bacmet" in sys.argv:
		titulo = "matriz_bacmet.csv"
		outputs.append(titulo)
#	elif "-megares" == p:
#		titulo = "matriz_megares.csv"

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
	outputs.append("Positions")
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
	fileType = "pdf"
	if "-pdf" in sys.argv:
		fileType = "pdf"
	elif "-png" in sys.argv:
		fileType = "png"
	if p == "-card":
		out = "clustermap_card."+fileType
		outputs.append(out)
		color = "Blues"
	elif p == "-vfdb":
		out = "clustermap_vfdb."+fileType
		outputs.append(out)
		color = "Reds"
	elif p == "-bacmet":
		out = "clustermap_bacmet."+fileType
		outputs.append(out)
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
	try:
		if (x > 1.9) and (y > 1.9):
			plt.figure()
			p2 = sns.clustermap(df, figsize=(x,y), cmap=color)
			p2.savefig(out, format=fileType, dpi=200, bbox_inches="tight")
		else:
			plt.subplots(figsize=(x,y))
			p2 = sns.heatmap(df, cmap=color, vmin=0, vmax=100)
			p2.figure.savefig(out, format=fileType, dpi=200, bbox_inches="tight")
	except:
		erro_string = "\nIt was not possible to plot the "+out+" figure...\nPlease verify the GenBank files and the matrix_x.csv output."
		erro.append(erro_string)
		print(erro_string)
#################################################Gerar figura###################################################

############################################Omics#########################################
	print("\nDoing presence analysis...")
	if p == "-card":
		t1 = "card_gene_count.csv"
		outputs.append(t1)
		t2 = "card_strain_count.csv"
		outputs.append(t2)
		t3 = "card_pan.csv"
		outputs.append(t3)
		t4 = "Pan-resistome development"
		l1 = "Pan-resistome"
		l2 = "Core-resistome"
		t6 = "card_mechanisms.csv"
		outputs.append(t6)
		t7 = "card_drug_classes.csv"
		outputs.append(t7)
		t8 = "card_mechanisms_barplot."+fileType
		t9 = "card_drug_classes_barplot."+fileType
	if p == "-vfdb":
		t1 = "vfdb_gene_count.csv"
		outputs.append(t1)
		t2 = "vfdb_strain_count.csv"
		outputs.append(t2)
		t3 = "vfdb_pan.csv"
		outputs.append(t3)
		t4 = "Pan-virulome development"
		l1 = "Pan-virulome"
		l2 = "Core-virulome"
		t6 = "vfdb_mechanisms.csv"
		outputs.append(t6)
		t7 = "vfdb_mechanisms_barplot."+fileType
	if p == "-bacmet":
		t1 = "bacmet_gene_count.csv"
		outputs.append(t1)
		t2 = "bacmet_strain_count.csv"
		outputs.append(t2)
		t3 = "bacmet_pan.csv"
		outputs.append(t3)
		t4 = "Pan-resistome development"
		l1 = "Pan-resistome"
		l2 = "Core-resistome"
		t6 = "bacmet_heavy_metals.csv"
		outputs.append(t6)
		t7 = "bacmet_all_compounds.csv"
		outputs.append(t7)
		t8 = "bacmet_heavy_metals_barplot."+fileType
	t5 = t3.replace("csv", fileType)
	outputs.append(t5)
#Per genes
	df2 = df
	headers = list(df2.columns.values)
	count = open(t1, "w")
	count.write("Genes;Presence Number;Strains\n")
	for i in headers:
		if "Unnamed:" in i:
			df2 = df2.drop(columns=[i])
		for j in df2[i]:
			if j != 0:
				df2[i] = df2[i].replace(j, 1)
		temp = df2[i].tolist()
		temp2 = ""
		for j in range(0, len(temp)):
			if temp[j] == 1:
				temp2 = temp2+str(lines[j])+','
		soma = 0
		soma = sum(df2[i])
		count.write(str(i)+";"+str(int(soma))+";"+temp2[:-1]+"\n")
	count.close()
#Per strains
	strains = list(df2.index.values)
	count = open(t2, "w")
	count.write("Strains;Presence Number;Genes\n")
	dic = {}
	for i in range(0, len(strains)):
		count.write(str(strains[i])+";")
		temp2 = []
		for j in range(0, len(headers)): #j são todos os genes
			if df2[headers[j]][i] == 1:
				temp2.append(headers[j])
		dic[strains[i]] = temp2
		count.write(str(len(temp2))+";")
		for k in temp2:
			count.write(str(k)+",")
		count.write('\n')
	count.close()
#PanOmic
	panomic = open(t3, "w")
	panomic.write("Strains;Core;Pan\n")
	pan = {}
	pan_atual = []
	core_atual = []
	for i in dic.keys():
		for j in dic[i]:
			if j not in pan_atual:
				pan_atual.append(j)
			if len(core_atual) == 0:
				core_atual = dic[i]
			elif len(core_atual) > 0:
				for k in core_atual:
					if k not in dic[i]:
						core_atual.remove(k)
		pan[i] = [len(core_atual), len(pan_atual)]
	for i in pan.keys():
		temp = str(i+";"+str(pan[i][0])+";"+str(pan[i][1])+"\n")
		panomic.write(temp)
	panomic.close()
	df = pd.read_csv(t3, sep=";")
	df["Number of Genomes"] = list(range(1,(len(df["Strains"].tolist()))+1))
	df = df.rename(columns={'Core': 'Number of Genes'})
	plt.figure()
	p1 = sns.lineplot(data=df, x="Number of Genomes", y="Pan", label=l1)
	p1 = sns.lineplot(data=df, x="Number of Genomes", y="Number of Genes", label=l2).set_title(t4)
	p1.figure.savefig(t5, format=fileType, dpi=300, bbox_inches="tight")
############################################Omics#########################################

############################################Pan_distribution#########################################
	if p == "-card":
		print("\nMaking the pan-distribution...")
		aro = pd.read_csv(dbpath+"aro_index.tsv", sep="\t")
		aro_genes = aro["ARO Name"].tolist()
		aro_drug = aro["Drug Class"].tolist()
		aro_mech = aro["Resistance Mechanism"].tolist()
		aro = {}
		for i in range(0, len(aro_genes)):
			aro[aro_genes[i]] = (aro_drug[i], aro_mech[i])

		mech = []
		drug = []
		for i in aro.keys():
			if ";" in aro[i][1]:
				temp = aro[i][1].split(";")
				for j in temp:
					if j not in mech:
						mech.append(j)
			else:
				temp = aro[i][1]
				if temp not in mech:
					mech.append(temp)
			if ";" in str(aro[i][0]):
				temp = aro[i][0].split(";")
				for j in temp:
					if j not in drug:
						drug.append(j)
			else:
				temp = aro[i][0]
				if temp not in drug:
					drug.append(temp)

		genes = pd.read_csv(t1, sep=";")
		core = []
		acce = []
		exclusive = []
		g = genes["Genes"].tolist()
		n = genes["Presence Number"].tolist()
		for i in range(0, len(g)):
			if n[i] == len(strains):
				core.append(g[i])
			elif (n[i]>1) and (n[i]<len(strains)):
				acce.append(g[i])
			elif n[i] == 1:
				exclusive.append(g[i])

		out = open(t6, "w")
		out.write("Resistance Mechanism;Core;Accessory;Exclusive\n")
		for k in mech:
			coreM = 0
			for i in core:
				try:
					if str(k) in aro[i][1]:
						coreM = coreM + 1
				except:
					for j in aro.keys():
						if i in j:
							if k in aro[j]:
								coreM = coreM + 1
			accessoryM = 0
			for i in acce:
				try:
					if str(k) in aro[i][1]:
						accessoryM = accessoryM + 1
				except:
					for j in aro.keys():
						if i in j:
							if k in aro[j]:
								accessoryM = accessoryM + 1
			exclusiveM = 0
			for i in exclusive:
				try:
					if str(k) in aro[i][1]:
						exclusiveM = exclusiveM + 1
				except:
					for j in aro.keys():
						if i in j:
							if k in aro[j]:
								exclusiveM = exclusiveM + 1
			if (coreM != 0) or (accessoryM != 0) or (exclusiveM != 0):
				out.write(str(k).capitalize()+";"+str(coreM)+";"+str(accessoryM)+";"+str(exclusiveM)+"\n")
		out.close()

		out2 = open(t7, "w")
		out2.write("Drug Class;Core;Accessory;Exclusive\n")
		for k in drug:
			if str(k).capitalize() == "Nan":
				continue
			coreM = 0
			for i in core:
				try:
					if str(k) in aro[i][0]:
						coreM = coreM + 1
				except:
					for j in aro.keys():
						if i in j:
							if k in aro[j]:
								coreM = coreM + 1
			accessoryM = 0
			for i in acce:
				try:
					if str(k) in aro[i][0]:
						accessoryM = accessoryM + 1
				except:
					for j in aro.keys():
						if i in j:
							if k in aro[j]:
								accessoryM = accessoryM + 1
			exclusiveM = 0
			for i in exclusive:
				try:
					if str(k) in aro[i][0]:
						exclusiveM = exclusiveM + 1
				except:
					for j in aro.keys():
						if i in j:
							if k in aro[j]:
								exclusiveM = exclusiveM + 1
			if (coreM != 0) or (accessoryM != 0) or (exclusiveM != 0):
				out2.write(str(k).capitalize()+";"+str(coreM)+";"+str(accessoryM)+";"+str(exclusiveM)+"\n")
		out2.close()

		try:
			data = pd.read_csv(t6, sep=";", index_col=("Resistance Mechanism"))
			x = len(list(data["Core"])) * 1
			y = 0
			for i in data:
				for j in range(0, len(data["Core"])):
					if data[i][j] >= y:
						y = data[i][j] / 3.5
			ax = data.plot.bar(figsize=(x,y), fontsize=15)
			ax.figure.savefig(t8, format=fileType, dpi=300, bbox_inches="tight")

			data = pd.read_csv(t7, sep=";", index_col=("Drug Class"))
			x = len(list(data["Core"])) * 1
			y = 0
			for i in data:
				for j in range(0, len(data["Core"])):
					if data[i][j] >= y:
						y = data[i][j] / 3.5
			ax = data.plot.bar(figsize=(x,y), fontsize=15)
			ax.figure.savefig(t9, format=fileType, dpi=300, bbox_inches="tight")
			outputs.append(t8)
			outputs.append(t9)
		except:
			erro_string = "\nIt was not possible to generate the barplot.\nPlease, check the input file.\n"
			print(erro_string)
			erro.append(erro_string)

	if p == "-vfdb":
		print("\nMaking the pan-distribution...")

		genes = pd.read_csv(t1, sep=";")
		core = []
		acce = []
		exclusive = []
		g = genes["Genes"].tolist()
		n = genes["Presence Number"].tolist()
		for i in range(0, len(g)):
			if n[i] == len(strains):
				core.append(g[i])
			elif (n[i]>1) and (n[i]<len(strains)):
				acce.append(g[i])
			elif n[i] == 1:
				exclusive.append(g[i])

		out = open(t6, "w")
		out.write("Virulence Mechanism;Core;Accessory;Exclusive\n")
		mech = []
		for k in genes_comp:
			if genes_comp[k] not in mech:
				mech.append(genes_comp[k])
		for mechanism in mech:
			core_number = 0
			accessory_number = 0
			exclusive_number = 0
			for gene in core:
				if genes_comp[gene] == mechanism:
					core_number = core_number + 1
			for gene in acce:
				if genes_comp[gene] == mechanism:
					accessory_number = accessory_number + 1
			for gene in exclusive:
				if genes_comp[gene] == mechanism:
					exclusive_number = exclusive_number + 1
			if (core_number != 0) or (accessory_number != 0) or (exclusive_number != 0):
				out.write(str(mechanism).capitalize()+";"+str(core_number)+";"+str(accessory_number)+";"+str(exclusive_number)+"\n")
		out.close()

		try:
			data = pd.read_csv(t6, sep=";", index_col=("Virulence Mechanism"))
			x = len(list(data["Core"])) * 1
			y = 0
			for i in data:
				for j in range(0, len(data["Core"])):
					if data[i][j] >= y:
						y = data[i][j] / 3.5
			ax = data.plot.bar(figsize=(x,y), fontsize=15)
			ax.figure.savefig(t7, format=fileType, dpi=300, bbox_inches="tight")
			outputs.append(t7)
		except:
			erro_string = "\nIt was not possible to generate the barplot.\nPlease, check the input file.\n"
			print(erro_string)
			erro.append(erro_string)
	
	if p == "-bacmet":
		db = pd.read_csv(dbpath+"bacmet_2.txt", sep="\t")
		genes = db["Gene_name"].tolist()
		compostos = db["Compound"].tolist()
		comp = []
		for i in compostos:
			if "," in i:
				temp = i.split(",")
				for j in temp:
					if j.strip() not in comp:
						comp.append(j.strip())
			else:
				if i.strip() not in comp:
					comp.append(i.strip())
		dic = {}
		for i in range(0, len(genes)):
			if "," not in compostos[i]:
				dic[genes[i]] = compostos[i]
			else:
				temp = compostos[i].split(', ')
				for j in range(0, len(temp)):
					temp[j] = temp[j].strip()
				dic[genes[i]] = temp
		matriz = pd.read_csv("bacmet_gene_count.csv", sep=";")
		my_genes = {}
		for i in range(0, len(matriz["Genes"].tolist())):
			my_genes[matriz["Genes"][i]] = matriz["Presence Number"][i]
		genes = {}
		core_ome = []
		accessory_ome = []
		exclusive_ome = []
		for i in my_genes.keys():
			if my_genes[i] == len(strains):
				core_ome.append(i)
			elif my_genes[i] > 1:
				accessory_ome.append(i)
			else:
				exclusive_ome.append(i)
		outa = open(t6, 'w')
		outa.write("Compound;Core;Accessory;Exclusive\n")
		outb = open(t7, 'w')
		outb.write("Compound;Core;Accessory;Exclusive\n")
		for k in comp:
			core = 0
			for i in core_ome:
				try:
					if k in dic[i]:
						core = core + 1
				except:
					for j in dic.keys():
						if i in j:
							if k in dic[j]:
								core = core + 1
			accessory = 0
			for i in accessory_ome:
				try:
					if k in dic[i]:
						accessory = accessory + 1
				except:
					for j in dic.keys():
						if i in j:
							if k in dic[j]:
								accessory = accessory + 1
			exclusive = 0
			for i in exclusive_ome:
				try:
					if k in dic[i]:
						exclusive = exclusive + 1
				except:
					for j in dic.keys():
						if i in j:
							if k in dic[j]:
								exclusive = exclusive + 1
			if (core != 0) or (accessory != 0) or (exclusive != 0):
				if ("(" in k) and ("[" not in k):
					outa.write(k+";"+str(core)+";"+str(accessory)+";"+str(exclusive)+"\n")
				outb.write(k+";"+str(core)+";"+str(accessory)+";"+str(exclusive)+"\n")
		outa.close()
		outb.close()
		
		try:
			data = pd.read_csv(t6, sep=";", index_col=("Compound"))
			x = len(list(data["Core"])) * 1
			y = 1
			for i in data:
				for j in range(0, len(data["Core"])):
					if data[i][j] >= y:
						y = data[i][j] / 3.5
			ax = data.plot.bar(figsize=(x,y), fontsize=15)
			ax.figure.savefig(t8, format=fileType, dpi=300, bbox_inches="tight")
			outputs.append(t8)
		except:
			erro_string = "\nIt was not possible to generate the barplot.\nPlease, check the input file.\n"
			print(erro_string)
			erro.append(erro_string)
################################################Pan_distribution#########################################

#####################################################Organizando#################################################
	print("\nGrouping the results...\n")

	atual = "Results_"+p[1:]+"_"+datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
	os.mkdir(atual)

	for i in outputs:
		try:
			shutil.move(i, atual)
		except:
			erro_string = "It was not possible to move the file "+str(i)+" to the final directory.\nPlease, chack the output path.\n"
			print(erro_string)
			erro.append(erro_string)

if "-b" in sys.argv:
	pasta = "GenBank_"+datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
	os.mkdir(pasta)
	for i in gbff:
		try:
			shutil.move(i, pasta)
		except:
			erro_string = "It was not possible to move the file "+str(i)+" to the final directory.\nPlease, chack the output path.\n"
			erro.append(erro_string)
			print(erro_string)
#################################################Organizando#################################################
if len(erro) > 0:
	final_erro = writeERR(erro)
	print("\nNumber of errors reported: "+str(len(erro)))
	print("Please check the "+str(final_erro)+" file.\n")
else:
	print("No error reported!\n")

#################################################Removendo arquivos intermediarios#################################################
if ("-keep" not in sys.argv) and ("-k" not in sys.argv):
	shutil.rmtree("Positions_1")
	shutil.rmtree("faa")
#################################################Removendo arquivos intermediarios#################################################

print('''That's all, folks!
Thank you for using this script.
If you're going to use results from PanViTa, please cite us:

PanViTa\thttps://doi.org/10.3389/fbinf.2023.1070406\t2023

Do not forget to quote the databases used.\n''')
if "-bacmet" in sys.argv:
	print("BacMet\thttps://doi.org/10.1093/nar/gkt1252\t2014")
if "-card" in sys.argv:
	print("CARD\thttps://doi.org/10.1093/nar/gkz935\t2020")
if "-vfdb" in sys.argv:
	print("VFDB\thttps://doi.org/10.1093/nar/gky1080\t2019")
print("\nDon't forget to mention the optional software too, if you already used them.")
if "-m" in sys.argv:
	print("mlst\thttps://doi.org/10.1186/1471-2105-11-595\t2010")
	print("mlst\thttps://github.com/tseemann/mlst")
if "-a" in sys.argv:
	print("prokka\thttps://doi.org/10.1093/bioinformatics/btu153\t2014")
print('')
