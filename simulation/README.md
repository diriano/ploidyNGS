# Simulate genomic short-reads for `ploidyNGS`

In order to test `ploidyNGS` with genomes of different ploidy levels we took Saccharomyces cerevisiae S288C chromosome I as basis. For this we developed `simulatePloidyData.py` that takes a chromose sequence in fasta format, an heterozigosity level and the ploidy level and produces a fasta file tht complies with the heterozugosity and ploidy level requested.

## Generation of the simulated diploid organism sequences.

```
python simulatePloidyData.py --genome GCA_000146045.2_R64_genomic_chromosomeI.fna --heterozygosity 0.01 --ploidy 2
```

This will result in the file `simulatedChroms_Ploidy2_Heter0.01.fasta`.

Here, `simulatePloidyData.py` takes the input chromosome (`GCA_000146045.2_R64_genomic_chromosomeI.fna`) and creates two chromosomes (ploidy 2) with heteromorphic positions found approx. every 100 bases (for the heterozigosity rate set, 0.01). In this case, the heteromorphic positions have the expected proportion for a diploid (50% of each allele).

## Generation of the simulated reads for the diploid organism.

For the generation of simulated HiSeq 2500 Illumina paired-end reads, we used the command below in [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/) (Huang, Weichun, *et al.* 2012):

```
art_illumina -na -i simulatedChroms_Ploidy2_Heter0.01.fasta -p -l 100 -ss HS25 -f 100 -m 200 -s 10 -o Ploidy2_100x
```

This will generate paired-end reads (`-p`) of 100 bps (`-l 100`) following the error profile of a HiSeq2500 (`-ss HS25`), with an insert size of 200bp (`-m 200) and a standard deviation of 10bp (`-s 10`). The resulting depth (coverage) should be 100x (`-f 100`). Resulting in the files with base name Ploidy2_100x.

The sofware outputs the reads in two files:
- Ploidy2_100x1.fq.gz
- Ploidy2_100x2.fq.gz

These reads should be used in mapping step.

## Mapping step using `Bowtie2`.

`ploidyNGS` requires a BAM with reads mapping to the genomic sequence. This file could be generated with Bowtie2 or other mapping software that produces BAM files.

Here, we have used the the original haploid *S. cerevisiae* S288C chromosome as reference, and the reads generated above with ART that came from a diploid sequence. We used `Bowtie2` (version 2.2.3) to map the reads to the reference as follows:

```
bowtie2 --met-file align_metrics.txt -t --very-sensitive -q -x GCA_000146045.2_R64_genomic_chromosomeI.fna -1 Ploidy2_100x1.fq -2 Ploidy2_100x2.fq | samtools view -bS - | samtools sort - Ploidy2.bowtie2.sorted
```

`samtools` (Li *et al.* 2009) is also used after the pipe in order to generate the BAM file and then sort it (by genome coordinates).

If you are working on a multi-threaded cluster, `Bowtie2` (and other mappers) also allows to use multiple threads, speeding up the process, please check its documentation.

The resulting BAM file could then be used by `ploidyNGS`.

# REFERENCES

Huang, Weichun, et al. "ART: a next-generation sequencing read simulator." *Bioinformatics* 28.4 (2012): 593-594.

Li, Heng, et al. "The sequence alignment/map format and SAMtools." *Bioinformatics* 25.16 (2009): 2078-2079.
