#!/usr/bin/env python 

"""This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import sys

#Making sure you are running a version of python that works with this script.
if sys.version_info[0] != 2 or sys.version_info[1] < 7 or sys.version_info[2] < 8:
    print("This script requires Python version 2.7.8 or higher within major version 2")
    sys.exit(1)
    
    
import argparse
import pysam
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os.path
import ploidyNGS_utils

###################################
# Dependencies:
# - python >= 2.7.8
# - pysam  >= 0.9
# - biopython >= 1.66
###################################

###################################
# Command line arguments
###################################

version=sys.argv[0] + ' ' + ploidyNGS_utils.git_version()
parser = argparse.ArgumentParser(description='ploidyNGS: Visual exploration of ploidy levels', add_help=True)
parser.add_argument('-v','--version', action='version', version=version)
parser.add_argument('-o','--out', dest='out', metavar='file', type=str, help='Base name for a TAB file that will keep the allele counts', required=True)
parser.add_argument('-b','--bam', dest='bam', metavar='mappingGenome.bam', type=str, help='BAM file used to get allele frequencies', required=True)
parser.add_argument('-m','--max_allele_freq', dest='AllowedMaxAlleleFreq', metavar='0.95 (default)', type=float, help='Fraction of the maximum allele frequency (float betwen 0 and 1, default: 0.95)', required=False, default=0.95)
parser.add_argument('-d','--max_depth', dest='MaxDepth', metavar='100 (default)', type=int, help='Max number of reads kepth at each position in the reference genome (integer, default: 100)', required=False, default=100)

# Get information from the argparse (arguments)
args = parser.parse_args()
bamOBJ = args.bam
AllowedMaxAlleleFreq = args.AllowedMaxAlleleFreq
MaxDepth=args.MaxDepth
baseOut=args.out + "_depth" + str(MaxDepth)
fileOut=baseOut + ".tab"
outOBJ = open(fileOut,"w")

# Check if bam index is present; if not, create it
bamindexname = args.bam + ".bai"
if os.path.isfile(bamindexname):
	print("BAM index present... OK!")
else:
	print("No index available for pileup. Creating an index...")
	pysam.index(args.bam)#TODO: check that indexing worked OK

# Get number of reads mapped, using idxstats, instead of samtools view, should be much faster, the ony drawback is that
#  it only gives the number of maped reads, independant of whether they were paired or not during mapping
print("Getting the number of mapped reads from BAM")
libSize=0
for l in pysam.idxstats(args.bam).split('\n'):
	if(len(l.split('\t'))==4):
		libSize= libSize+int(l.split('\t')[2])  

print libSize

# Create a pysam object for the indexed BAM
bamfile = pysam.AlignmentFile(bamOBJ, "rb")

# Create a hash of hashes to count number of each allele at each chromosome position
def makehash():
	return defaultdict(makehash)

count = makehash()
countAlleleNormalized = makehash()

# Traversing BAM file: Count the number of reads for 
#  each observed nucleotide at each position in the chromosome/contig
#  Stores the count in a dictionary of dictionaries (count) for later processing
for contig in bamfile.references:
	for pucolumn in bamfile.pileup(contig, 0):
		pos_1 = pucolumn.pos
		countReadsPos=0
		#print("Contig: " + contig + " Chromosome Position:" + str(pos_1))
		for puread in pucolumn.pileups:
			if not puread.is_del and not puread.is_refskip and countReadsPos <= MaxDepth-1:
				countReadsPos=countReadsPos+1
				base = puread.alignment.query_sequence[puread.query_position]
				if count[contig][pos_1][base]:
					count[contig][pos_1][base]=count[contig][pos_1][base]+1
				else:
					count[contig][pos_1][base]=1
		#print("Total number of reads: " + str(countReadsPos))

#Traversing dictionary of dictionaries with number of reads for each observed nucleotide
# at each position, skips monomorphic positions and positions in which the most frequent
# nucleotide has a frequency larger than AllowedMaxAlleleFreq.
for contig, dict2 in count.iteritems():
	for pos, dict3 in dict2.iteritems():
		pos_depth = 0
		pos_bases = {}
		if count[contig][pos]['A']:
			pos_depth = pos_depth + int(count[contig][pos]['A'])
			pos_bases['A']=int(count[contig][pos]['A'])
			#print("A: %s" % (count[contig][pos]['A']))
		if count[contig][pos]['C']:
			pos_depth = pos_depth + int(count[contig][pos]['C'])
			pos_bases['C']=int(count[contig][pos]['C'])
			#print("C: %s" % (count[contig][pos]['C']))
		if count[contig][pos]['T']:
			pos_depth = pos_depth + int(count[contig][pos]['T'])
			pos_bases['T']=int(count[contig][pos]['T'])
			#print("T: %s" % (count[contig][pos]['T']))
		if count[contig][pos]['G']:
			pos_depth = pos_depth + int(count[contig][pos]['G'])
			pos_bases['G']=int(count[contig][pos]['G'])
			#print("G: %s" % (count[contig][pos]['G']))
		#print(pos_depth)
		if len(pos_bases) > 1: #Skips monomorphic positions
			#print "Max Allele", max(pos_bases, key=pos_bases.get)
			maxAlleleFreq = (float(pos_bases[max(pos_bases, key=pos_bases.get)]))/float(pos_depth)
			if maxAlleleFreq <= AllowedMaxAlleleFreq: #Skips positions in which the most frequent nucleotide has a frequency larger than AllowedMaxAlleleFreq
				for obsBase, obsCount in pos_bases.iteritems():
					percBase = (float(obsCount) / float(pos_depth)) * 100
					countAlleleNormalized[contig][pos][obsBase]=percBase
				alleleFreqDist = []
				if(countAlleleNormalized[contig][pos]['A']):
					#outOBJ.write("%s\t" % (countAlleleNormalized[contig][pos]['A']))
					alleleFreqDist = alleleFreqDist + [countAlleleNormalized[contig][pos]['A']]
				else:
					#outOBJ.write("0\t")
					alleleFreqDist = alleleFreqDist + [0]
				if(countAlleleNormalized[contig][pos]['T']):
					#outOBJ.write("%s\t" % (countAlleleNormalized[contig][pos]['T']))
					alleleFreqDist = alleleFreqDist + [countAlleleNormalized[contig][pos]['T']]
				else:
					#outOBJ.write("0\t")
					alleleFreqDist = alleleFreqDist + [0]
				if(countAlleleNormalized[contig][pos]['C']):
					#outOBJ.write("%s\t" % (countAlleleNormalized[contig][pos]['C']))
					alleleFreqDist = alleleFreqDist + [countAlleleNormalized[contig][pos]['C']]
				else:
					#outOBJ.write("0\t")
					alleleFreqDist = alleleFreqDist + [0]
				if(countAlleleNormalized[contig][pos]['G']):
					#outOBJ.write("%s\t" % (countAlleleNormalized[contig][pos]['G']))
					alleleFreqDist = alleleFreqDist + [countAlleleNormalized[contig][pos]['G']]
				else:
					#outOBJ.write("0\t")
					alleleFreqDist = alleleFreqDist + [0]
				alleleFreqDist.sort()
				outOBJ.write("%s\t%s\tFourth\t%.2f" % (contig, pos, alleleFreqDist[0]))
				outOBJ.write("\n")
				outOBJ.write("%s\t%s\tThird\t%.2f" % (contig, pos, alleleFreqDist[1]))
				outOBJ.write("\n")
				outOBJ.write("%s\t%s\tSecond\t%.2f" % (contig, pos, alleleFreqDist[2]))
				outOBJ.write("\n")
				outOBJ.write("%s\t%s\tFirst\t%.2f" % (contig, pos, alleleFreqDist[3]))
				outOBJ.write("\n")

outOBJ.close()

cmdRscript="Rscript --vanilla ploidyNGS_generateHistogram.R "+ fileOut + " " + fileOut +".PloidyNGS.pdf" + " " + str(1-AllowedMaxAlleleFreq) + " " + str(AllowedMaxAlleleFreq)
#print(cmdRscript)
os.system(cmdRscript)
