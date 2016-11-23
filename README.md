# ploidyNGS
Exploring ploidy levels from NGS data alone.

**Motivation:** Due to the decreasing costs of genome sequencing researchers are increasingly assessing the genome information of non-model organisms using Next Generation Sequencing technologies. An important question to answer when assessing a new genome sequence is: What is the ploidy level of the organism under study?

**Results:** We have developed \texttt{ploidyNGS}, a model-free, open source tool to visualize and explore ploidy levels in a newly sequenced genome, exploiting short read data. We tested \texttt{ploidyNGS} using both simulated and real NGS data of the model yeast \textit{Saccharomyces cerevisiae}. \texttt{ploidyNGS} allows the identification of the ploidy level of a newly sequenced genome in a visual way. We have applied \texttt{ploidyNGS} to a wide range of different genome ploidy and heterozigosity levels, as well as to a range of sequencing depths.

**Availability and implementation:** `ploidyNGS` is available under the GNU General Public License (GPL) at https://github.com/diriano/ploidyNGS. ploidyNGS is implemented in Python and R.


## Requirements

- python >=2.7.8
- pysam  >=0.9
- biopython >=1.66
- R and ggplot2

# Examples and test
-------------------

## Interpretation of results of ploidyNGS with different ploidy levels

In the directory test_data, we provide examples of ploidyNGS results.
The resulting pdf files for each ploidyNGS run are found in directory test_data/ploidyNGS_results.

- `DataTestPloidy1.tab.ExplorePloidy.pdf` - haploid organism
- `DataTestPloidy2.tab.ExplorePloidy.pdf` - diploid
- `DataTestPloidy3.tab.ExplorePloidy.pdf` - triploid
- `DataTestPloidy4.tab.ExplorePloidy.pdf` - tetraploid

The interpretation of the resulting plots is based on the frequency of the proportion of each allele, in all observed heteromorphic sites (which corresponds to the x axis in the plot) for a given sequenced organism. A table with the expected bi-allelic proportions for each ploidy level is shown below.

| Ploidy level | Genome Position | Allele 1 (expected proportion) | Allele 2 (expected proportion) |
| ------------ | --------------- |-------------------------------:|-------------------------------:|
| Haploidy | Monomorphic | 0 | 100 |
| Diploidy | Monomorphic | 0 | 100 |
| Diploidy | Heteromorphic | 50 | 50 |
| Triploidy | Monomorphic | 0 | 100 |
| Triploidy | Heteromorphic | 33.33 | 66.67 |
| Tetraploidy | Monomorphic | 0 | 100 |
| Tetraploidy | Heteromorphic | 25 | 75 |
| Tetraploidy | Heteromorphic | 50 | 50 |
| Pentaploidy | Monomorphic | 0 | 100 |
| Pentaploidy | Heteromorphic | 20 | 80 |
| Pentaploidy | Heteromorphic | 40 | 60 |
| Hexaploidy | Monomorphic | 0 | 100 |
| Hexaploidy | Heteromorphic | 16.67 | 83.33 |
| Hexaploidy | Heteromorphic | 50 | 50 |
| Hexaploidy | Heteromorphic | 33.33 | 66.67 |
| Heptaploidy | Monomorphic | 0 | 100 |
| Heptaploidy | Heteromorphic | 28.57 | 71.43 |
| Heptaploidy | Heteromorphic | 42.86 | 57.14 |
| Heptaploidy | Heteromorphic | 14.29 | 85.71 |

For example, in the diploid case, we expect only one peak with frequencies of two alleles, both representing 50% of the observations, as observed in figure `DataTestPloidy2.tab.ExplorePloidy.pdf`. Other observed proportions in this plot correspond to noise due to sequencing errors.

## Full analysis example - diploid organism

Here is a complete example of a pipeline for ploidy analysis using `ploidyNGS`.
For this example, we use a simulated diploid yeast chromosome generated using our script `simulatePloidyData.py`, which results in the plot above, `test_data/ploidyNGS_results/DataTestPloidy2.tab.ExplorePloidy.pdf`.

* NOTE: For real data, the user will have only reads from sequencing, not from simulated genome sequence.
* For real data, the BAM used in ploidyNGS should be created from mapping reads to the assembled genome sequence, or to a reference sequence (e.g. a closely related strain).

a) Generation of the simulated diploid organism sequences.

We used the *Saccharomyces cerevisiae* S288c chromosome I sequence (haploid) to generate the diploid one, as shown here:

```
$ python3 simulatePloidyData.py --genome GCA_000146045.2_R64_genomic_chromosomeI.fna --heterozygosity 0.01 --ploidy 2
```

Here, `simulatePloidyData.py` takes the input chromosome (`GCA_000146045.2_R64_genomic_chromosomeI.fna`) and creates two chromosomes (ploidy 2) with heteromorphic positions very around 100 bases (for the heterozigosity rate set, 0.01). In this case, the heteromorphic positions have the expected proportion for a diploid (50% of each allele).

b) Generation of the simulated reads for the diploid organism.

For the generation of simulated HiSeq 2500 Illumina paired-end reads, we used the command below in ART (Huang, Weichun, *et al.* 2012):

```
$ art_illumina -na -i simulatedChroms_Ploidy2_Heter0.01.fasta -p -l 100 -ss HS25 -f 100 -m 200 -s 10 -o Ploidy2_100x
```

The paired-end reads are generated with 100 bases and an average fragment size of 200.

The `.fastq` files are in the directory `test_data/simulatedDiploidGenome` :

The sofware outputs the reads in two files:
- Ploidy2_100x1.fq.gz
- Ploidy2_100x2.fq.gz

These reads are used in mapping step.

c) Mapping step using `Bowtie2`.

A BAM with reads mapping the genomic sequence is required in ploidyNGS main script.
It can be generated with Bowtie2 or other mapping algorithm.

Here, we used the haploid *S. cerevisiae* chromosome I to map our reads from the diploid. Here is how `Bowtie2` (version 2.2.3) is run:

```
$ bowtie2 --met-file align_metrics.txt -t --very-sensitive -q -x GCA_000146045.2_R64_genomic_chromosomeI.fna -1 Ploidy2_100x1.fq -2 Ploidy2_100x2.fq | samtools view -bS - | samtools sort - Ploidy2.bowtie2.sorted
```

`samtools` (Li *et al.* 2009) is also used after the pipe in order to generate the BAM file and then sort it (by genome coordinates).

If you are working on a multi-threaded cluster, `Bowtie2` (and other mappers) also allows to use multiple threads, speeding up the process.

The resulting BAM file is then used in ploidyNGS script (next step).

d) Running ploidyNGS main script (`explorePloidyNGS.py`)

`explorePloidyNGS.py` only requires a BAM and a name for the output file. Here is how to use it:

```
$ explorePloidyNGS.py --out DataTestPloidy2.tab --bam Ploidy2.bowtie2.sorted.bam 
```

The script outputs two files:
- `DataTestPloidy4.tab.ExplorePloidy.pdf`, the plot with the frequencies of allele proportions in heterozygous sites.
- `DataTestPloidy2.tab.Rscript`, the R script used to generate the plot above. 

# REFERENCES

Huang, Weichun, et al. "ART: a next-generation sequencing read simulator." *Bioinformatics* 28.4 (2012): 593-594.

Li, Heng, et al. "The sequence alignment/map format and SAMtools." *Bioinformatics* 25.16 (2009): 2078-2079.
