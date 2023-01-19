#!/usr/bin/env python 

#Visually exploring ploidy with Next Generation Sequencing data
__author__      = "Diego Mauricio Riano-Pachon & Renato Augusto Correa dos Santos"
__copyright__   = "Copyright 2016,2017"
__license__     = "GPL v3.0"
__maintainer__  = "Diego Mauricio Riano-Pachon"
__email__       = "diriano@gmail.com"

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
if sys.version_info[0] != 3:
    print("This script requires Python version 3")
    sys.exit(1)
    
    
import argparse
import pysam
from collections import defaultdict
import os.path
import linecache
from time import  asctime
import ploidyNGS_utils

###################################
# Dependencies:
# - python >= 3 or higher
# - pysam  >= 0.17 or higher
# - biopython >= 1.79 or higher
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
parser.add_argument('-u','--min_cov', dest='MinCov', metavar='0 (default)', type=int, help='Minimum coverage required to use a position (integer, default: 0)', required=False, default=0)
parser.add_argument('-g','--guess_ploidy', dest='guessPloidy', help='Try to guess ploidy level by comparison with simulated data', required=False, action="store_true")
parser.add_argument('-c','--coverage_guess_ploidy', dest='covGuessPloidy', choices=[15,25,50,100], help='This parameter will be automatically set based on the coverage of your BAM file, however you can overrride it.', required=False, type=int)

#Print a greeting with the date and the version number
print("###############################################################")
print("## This is ploidyNGS version " + ploidyNGS_utils.git_version())
print("## Current date and time: " + asctime())
print("###############################################################")
# Get information from the argparse (arguments)
args = parser.parse_args()
bamOBJ = args.bam
AllowedMaxAlleleFreq = args.AllowedMaxAlleleFreq
MaxDepth=args.MaxDepth
MinCov=args.MinCov
baseOut=args.out + "_MaxDepth" + str(MaxDepth) + "_MinCov" + str(MinCov)
fileHistOut=baseOut + ".tab"
fileGuessOut=baseOut + ".ks-distance.PloidyNGS.tbl"
covGuessPloidy=args.covGuessPloidy
outOBJ = open(fileHistOut,"w")
# Check if bam index is present; if not, create it
bamindexname = args.bam + ".bai"
if os.path.isfile(bamindexname):
	print("BAM index present... OK!")
else:
	print("No index available for pileup. Creating an index...")
	pysam.index(args.bam)#TODO: check that indexing worked OK

# Get number of reads mapped, using idxstats, instead of samtools view, should be much faster, the ony drawback is that
#  it only gives the number of maped reads, independant of whether they were paired or not during mapping
libSize=0
for l in pysam.idxstats(args.bam).split('\n'):
	if(len(l.split('\t'))==4):
		libSize= libSize+int(l.split('\t')[2])  

print("Number of mapped reads from BAM: %i" % libSize)

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

countTotalPositions=0
countTotalReads=0

for contig in bamfile.references:
	for pucolumn in bamfile.pileup(contig, 0):

		# skip if low coverage position
		## .get_num_aligned() returns number of reads aligned to
		## this position after a min base quality filter is applied
		## alternatively use .nsegments to ignore min base quality filter
		if pucolumn.get_num_aligned() <= MinCov:
			continue

		countTotalPositions=countTotalPositions+1
		# add 1 to position count because .pos returns 0-indexed position
		pos_1 = pucolumn.pos + 1
		countReadsPos=0
		#print("Contig: " + contig + " Chromosome Position:" + str(pos_1))
		for puread in pucolumn.pileups:
			if not puread.is_del and not puread.is_refskip and countReadsPos <= MaxDepth-1:
				countTotalReads=countTotalReads+1
				countReadsPos=countReadsPos+1
				base = puread.alignment.query_sequence[puread.query_position]
				if count[contig][pos_1][base]:
					count[contig][pos_1][base]=count[contig][pos_1][base]+1
				else:
					count[contig][pos_1][base]=1
		#print("Total number of reads: " + str(countReadsPos))

averageCoverage=countTotalReads/countTotalPositions
print("Observed average coverage: %5.2f" % averageCoverage)

#Traversing dictionary of dictionaries with number of reads for each observed nucleotide
# at each position, skips monomorphic positions and positions in which the most frequent
# nucleotide has a frequency larger than AllowedMaxAlleleFreq.

# count total number of heteromorphic positions per contig
# that is, sites with more than 1 base, and max frequency < MaxAlleFreq
countHeteroPosPerContig = {}

for contig, dict2 in count.items():

	# set heteromorphic counter to 0
	countHeteroPosPerContig[contig] = 0

	# traverse dict of dicts
	for pos, dict3 in dict2.items():
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
				
				# add heteromorphic site to counter
				countHeteroPosPerContig[contig]=countHeteroPosPerContig[contig]+1

				# calculate per-base frequencies and store in list
				for obsBase, obsCount in pos_bases.items():
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

# print total number of heteromorphic sites per contig
for contig, count in countHeteroPosPerContig.items():
	print("Number of heteromorphic positions in ", contig, ": ", count)
# and overall total
print("Total number of heteromorphic positions: ", sum( countHeteroPosPerContig.values() ))

cmdPloidyGraphRscript="Rscript --vanilla ploidyNGS_generateHistogram.R "+ fileHistOut + " " + fileHistOut + ".PloidyNGS.pdf" + " " + str(1-AllowedMaxAlleleFreq) + " " + str(AllowedMaxAlleleFreq)
#print(cmdPloidyGraphRscript)
os.system(cmdPloidyGraphRscript)

if(args.guessPloidy):
 setCoverage=0
 if(covGuessPloidy):
  setCoverage=covGuessPloidy
 else:
  if(averageCoverage < 20):
   setCoverage=15
  elif(averageCoverage >= 20 and averageCoverage < 37.5):
   setCoverage=25
  elif(averageCoverage >= 37.5 and averageCoverage < 75):
   setCoverage=50
  elif(averageCoverage >= 75):
   setCoverage=100
  else:
   setCoverage=100
 print("Coverage used for guessing ploidy: %i" % setCoverage) 
 if(os.path.isdir("./simulation/data")):
  cmdGuessPloidyRscript="Rscript --vanilla ploidyNGS_guessPloidy.R " + str(setCoverage) + " " + fileHistOut + " " + fileGuessOut
  #print(cmdGuessPloidyRscript)
  os.system(cmdGuessPloidyRscript)
  line=linecache.getline(fileGuessOut, 2)
  fields=line.replace('"','').split(" ")
  res="""
  After comparing your data with our simulated dataset
  and computing the Kolmogorov-Smirnov distance, 
  the closest ploidy to yours is %s
  """
  print(res % (fields[0]))
 else:
  print("The directoy with simulation data is missing. Guessing ploidy level cannot run!")
