import argparse, csv, os, sys

parser=argparse.ArgumentParser(description="Determines the index of individuals\
 of a given ethnicity within a 1000 Genomes VCF")
parser.add_argument("VCF", type=str, help="VCF of 1000 Genomes individuals")
parser.add_argument("Population", type=str, help="1KG identifier for population\
to be found in the index, enter \"ALL|\" to print index for all populations in \
the VCF")
arg=parser.parse_args()

a=os.path.abspath(sys.argv[0]).split("/")[:-1]
b="/".join(a)

try:
	 open(b+"/1KG_IDs.txt")
except IOError:
	print "Missing file \"1KG_IDs.txt\" containing the IDs for 1000 Genomes"

popdict={}
with open(b+"/1KG_IDs.txt") as t:
	for line in csv.reader(t,delimiter='\t'):
		popdict[line[1]]=line[0]

with open(arg.VCF) as t:
	for line in csv.reader(t,delimiter='\t'):
		if line[0]=='#CHROM':
			for j in range(len(line)):
				if 'NA' in line[j] or 'HG' in line[j]:
					vcfind=line[j:]
					break

indict={}
for i in range(len(vcfind)):
	try:
		eth=popdict[vcfind[i]]
	except:
		print"Non 1000 Genomes individuals found in VCF"
	if eth in indict.keys():
		indict[eth].append(i)
	else:
		indict[eth]=[i]

if arg.Population=="All":
	for j in indict.keys():
		l=[str(i) for i in indict[j]]
		print j,"=",",".join(l)

else:
	try:
		indict[arg.Population]
	except IOError:
		print "Population not in VCF or invalid population ID"
	l=[str(i) for i in indict[arg.Population]]
	print ",".join(l)
