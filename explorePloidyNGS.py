#!/bioinf/progs/anaconda/bin/python

import argparse
import pysam
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os.path

###################################
# Dependences:
# - samtools >= 1.2
# - python >=2.7.8
# - pysam >0.9
###################################

###################################
# Global variables
###################################

chroms_dict = defaultdict(list)
pathSamtoolsBin="/bioinf/progs/samtools-1.3/bin/samtools"

###################################
# Command line arguments
###################################

parser = argparse.ArgumentParser(description='Print allele frequencies per position', add_help=True)
parser.add_argument('--out', dest='out', metavar='file.tab', type=str, help='TAB file with allele counts', required=True)
parser.add_argument('--bam', dest='bam', metavar='mappingGenome.bam', type=str, help='BAM file used to get allele frequencies', required=True)
parser.add_argument('--genome', dest='genome', metavar='genome.fasta', type=str, help='Genome (FASTA)', required=True)
parser.add_argument('--max_allele_freq', dest='AllowedMaxAlleleFreq', metavar='0.95 (default)', type=float, help='Percentage of the maximum allele frequency', required=False, default=0.95)

# Get information from the argparse (arguments)
args = parser.parse_args()
bamOBJ = open(args.bam,"r")
outOBJ = open(args.out,"w")
genomeOBJ = open(args.genome,"r")
AllowedMaxAlleleFreq = args.AllowedMaxAlleleFreq

# get chromosome sizes from the FASTA
chromsSize = {} 

for chr in SeqIO.parse(genomeOBJ, "fasta"):
	chromsSize[chr.id]=len(chr.seq)

# Check if bam index is present; if not, create it
bamindexname = args.bam + ".bai"
if os.path.isfile(bamindexname):
	print("BAM index present... OK!")
else:
	print("No index available for pileup. Creating an index...")
	cmdSamtololsIndex=pathSamtoolsBin + " index " + args.bam
	samtololsIndexReturn=os.system(cmdSamtololsIndex)
	if (samtololsIndexReturn!=0):
		raise Exception("Problem indexing your bam file")
	else:
		print(". . . index succesfully created!")

# Get number of reads mapped, using idxstats, instead of samtools view, should be much faster, the ony draw back is that it only give the number of maped reads, independant of whether they were paired or not during mapping
print("Getting the number of mapped reads from BAM")
libSize=0
for l in pysam.idxstats(args.bam).split('\n'):
	if(len(l.split('\t'))==4):
		libSize= libSize+int(l.split('\t')[2])  

print libSize
#libSize= reduce(lambda x, y: x + y, [ int(l.split('\t')[2]) for l in pysam.idxstats(args.bam).split('\n') ])

# Create a pysam object for the indexed BAM
bamfile = pysam.AlignmentFile(bamOBJ, "rb")

# Create a hash of hashes to count number of each allele at each chromosome position
def makehash():
	return defaultdict(makehash)

count = makehash()
countAlleleNormalized = makehash()

# Count each allele at each position in chromosome
for contig, contigSize in chromsSize.iteritems():
		for pucolumn in bamfile.pileup(contig, 0):
			pos_1 = pucolumn.pos
			for puread in pucolumn.pileups:
				if not puread.is_del and not puread.is_refskip:
					base = puread.alignment.query_sequence[puread.query_position]
					if count[contig][pos_1][base]:
						count[contig][pos_1][base]=count[contig][pos_1][base]+1
					else:
						count[contig][pos_1][base]=1

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
		if len(pos_bases) > 1:
			#print "Max Allele", max(pos_bases, key=pos_bases.get)
			maxAlleleFreq = (float(pos_bases[max(pos_bases, key=pos_bases.get)]))/float(pos_depth)
			if maxAlleleFreq <= AllowedMaxAlleleFreq:
				pos_depth_cmp = (float(pos_depth) / float(libSize)) * 1000000
				for obsBase, obsCount in pos_bases.iteritems():
					normalized = (float(obsCount) / float(libSize)) * 1000000
					percBase = (normalized / pos_depth_cmp) * 100
					countAlleleNormalized[contig][pos][obsBase]=percBase
					#print pos, maxAlleleFreq, pos_depth, pos_depth_cmp, obsBase, obsCount, percBase, libSize

for contig, dict2 in countAlleleNormalized.iteritems():
        for pos, dict3 in dict2.iteritems():
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
		outOBJ.write("%s\t%s\tFourthFreq\t%.2f" % (contig, pos, alleleFreqDist[0]))
		outOBJ.write("\n")
		outOBJ.write("%s\t%s\tThirdFreq\t%.2f" % (contig, pos, alleleFreqDist[1]))
		outOBJ.write("\n")
		outOBJ.write("%s\t%s\tSecondFreq\t%.2f" % (contig, pos, alleleFreqDist[2]))
		outOBJ.write("\n")
		outOBJ.write("%s\t%s\tFirstFreq\t%.2f" % (contig, pos, alleleFreqDist[3]))
		outOBJ.write("\n")

outOBJ.close()
bamOBJ.close()

def createRscript(table):
	rfile=table + ".Rscript"
	pdfFilename=table + ".ExplorePloidy.pdf"
	title="Explore ploidy - NGS"
	rscriptOBJ = open(rfile,"w")
        rscriptOBJ.write("datain<-read.table(\"table\",header=F)")
	rscriptOBJ.write("colnames(datain)<-c('Chrom','Pos','Type','Freq')")
	rscriptOBJ.write("head(datain)")
	rscriptOBJ.write("dim(datain)")
	rscriptOBJ.write("pdf(\"pdfFilename\")")
	rscriptOBJ.write("ggplot(datain,aes(x=Freq, fill=Type)) +")
	rscriptOBJ.write(" geom_histogram(binwidth = 0.5, alpha=0.4) +")
	rscriptOBJ.write(" ggtitle(\"title\") +")
	rscriptOBJ.write(" ylab(\"Counts positions\") +")
	rscriptOBJ.write(" xlab(\"Allele Freq\") +")
	rscriptOBJ.write(" scale_x_continuous(limits=c(1,100))")
	rscriptOBJ.write("dev.off()")
	return;

createRscript(args.out)
