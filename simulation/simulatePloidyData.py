#!/progs/users/python-2.7.8/bin//python

#Requires:
#bio-python: in Ubuntu 14.04 run: sudo apt-get install python-biopython
import argparse
import random
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os.path

###################################
# Global variables
###################################
# Create an empty array to store the dosage table
DosageTable = []

###################################
# Command line arguments
###################################
parser = argparse.ArgumentParser(description='Creates a simulated polyploid genome')
#Add license arg
parser.add_argument('--genome', '-g', metavar='haploidGenome.fasta', type=str, required=True, help='An haploid genome (fasta)', dest='genomeFile')
parser.add_argument('--ploidy', '-p', metavar='ploidy', type=int, required=True, help='The ploidy of the final genome (an integer)', dest='ploidy')
parser.add_argument('--heterozygosity', metavar='heterozygosity', type=float, required=True, help='The heterozygosity rate of resulting genome, a number betwen 0 (low) and 1 (highest)', dest='heteroz')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

args = parser.parse_args()

"""
# Checks if user requested the license to be displayed
if lic:
	print("Copyright (C) 2016 Diego Mauricio Riano Pachon\ne-mail: diriano\@gmail.com Copyright (c) 2016\nThis program is free software; you can redistribute it and/or\nmodify it under the terms of the GNU General Public License\nas published by the Free Software Foundation; either version 2\nof the License, or (at your option) any later version.")
	print("This program is distributed in the hope that it will be useful,\nbut WITHOUT ANY WARRANTY; without even the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\nGNU General Public License for more details.")
	print("You should have received a copy of the GNU General Public License\nalong with this program; if not, write to the Free Software\nFoundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.")
	exit()
"""

# Get dosage and heterozygosity from user args
ploidy = args.ploidy
heterozygosity = args.heteroz

#Check existence of output files
for k in range(ploidy):
	fname = "simulatedChroms_Ploidy" + str(ploidy) + "_Heter" + str(heterozygosity) + "_" + str(k) + ".fasta"
	#print(fname)
	if os.path.isfile(fname):
		print("Output file (%s) already exists!" % fname)
		exit()

# Generate the dosage table for a given ploidy (assuming that only two alleles are present in the individual (biallelic), which is the most common case)
# TODO: Extend to multiallelic
j = ploidy

for i in range(ploidy):
	if (j >= i) and (i != 0):
		DosageTable.append([i,j])
	j -= 1

# Get allele dosage table for further use
lenDosageTable = len(DosageTable)

def makehash():
	return defaultdict(makehash)

# Create tables that will store dosage and the second allele in each position
heterozPosChromsDosageBases = makehash()

genomeOBJ = open(args.genomeFile, "rU")

# Get random dosage (dosage table index) and second allele base (a random base, besides the original one)
for chr in SeqIO.parse(genomeOBJ, "fasta"):
	for index, base in enumerate(chr):
		randomUnif = random.uniform(0,1)
		if randomUnif <= heterozygosity:
			randomDosageTable = random.randint(0,lenDosageTable-1)
			DosageThisPosition = DosageTable[randomDosageTable]
			base = base.upper()
			GenomeAlph = ["A","T","C","G"]
			GenomeAlph.remove(base)
			baseAltAllele = GenomeAlph[random.randint(0,2)]
			ListFirstBase = DosageThisPosition[0] * base
			ListSecBase = DosageThisPosition[1] * baseAltAllele
			TotalBasesThisPosition = ListFirstBase + ListSecBase
			heterozPosChromsDosageBases[chr.id][index]=TotalBasesThisPosition

genomeOBJ.close()

# For each position, check whether it's polymorphic, then assigns a new or old (original) base, according to the randomly selected dosage
for i in range(1, ploidy+1): #Python's range function will generate a sequence of number from start (1) up to, but not including stop (ploidy+1)
	filename = "simulatedChroms_Ploidy" + str(ploidy) + "_Heter" + str(heterozygosity) + "_" + str(i) + ".fasta"
	output_handle = open(filename, "w")
	genomeOBJ = open(args.genomeFile, "rU")
	for chr in SeqIO.parse(genomeOBJ, "fasta"):
		print(chr.id)
		newChrID = chr.id + "_p" + str(ploidy) + "_h" + str(heterozygosity) + "_copy" + str(i)
		print(newChrID)
		descPloidyCopyNum = "Chr copy number " + str(i)
		newChr = SeqRecord(Seq('', IUPAC.unambiguous_dna), id=newChrID, description=descPloidyCopyNum)
		for index, base in enumerate(chr):
			if index in heterozPosChromsDosageBases[chr.id]:
				newChr = newChr + heterozPosChromsDosageBases[chr.id][index][i-1]
			else:
				newChr = newChr + base
		SeqIO.write(newChr, output_handle, "fasta")
	output_handle.close()

	genomeOBJ.close()
