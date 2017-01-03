# ploidyNGS: Visually exploring ploidy with Next Generation Sequencing data

**Motivation:** Due to the decreasing costs of genome sequencing researchers are increasingly assessing the genome information of non-model organisms using Next Generation Sequencing technologies. An important question to answer when assessing a new genome sequence is: What is the ploidy level of the organism under study?

**Results:** We have developed `ploidyNGS`, a model-free, open source tool to visualize and explore ploidy levels in a newly sequenced genome, exploiting short read data. We tested `ploidyNGS` using both simulated and real NGS data of the model yeast *Saccharomyces cerevisiae*. `ploidyNGS` allows the identification of the ploidy level of a newly sequenced genome in a visual way. We have applied `ploidyNGS` to a wide range of different genome ploidy and heterozigosity levels, as well as to a range of sequencing depths.

**Availability and implementation:** `ploidyNGS` is available under the GNU General Public License (GPL) at https://github.com/diriano/ploidyNGS. ploidyNGS is implemented in Python and R.

## Requirements

- python >=2.7.8
- pysam  >=0.9
- biopython >=1.66
- R and ggplot2

## Installation

### Get ploidyNGS from github

```bash
cd ~
git clone https://github.com/diriano/ploidyNGS.git
```
### Python virtual environment and python dependencies

Make sure pip is installed. If not, then please check your distribution's documentation and install it. In Ubuntu 16.04 LTS you can do:

```bash
sudo apt install python-pip
```

Install the virtualenv module:

```bash
pip install --upgrade pip
pip install virtualenv
```

Go to the ploidyNGS folder, create and active a new virtual environment

```bash
cd ploidyNGS
virtualenv .venv
source .venv/bin/activate
```

Install dependencies:

pysam requires the zlib headers. Please make sure you have these installed. In Ubuntu 16.04 LTS you can do:

```bash
sudo apt install zlib1g-dev
```

Then install pysam:

```bash
pip install pysam
```

Install biopython dependencies, and then biopython itself (for details see: http://biopython.org/DIST/docs/install/Installation.html):

```bash
pip install numpy
pip install biopython
```

At this point you can deactive the python virtual environment:

```bash
deactivate
```
### R

In Ubuntu 16.04 LTS:

```bash
  sudo apt install r-base
  sudo apt install r-cran-ggplot2
```

### Test ploidyNGS

In order to use `ploidyNGS` please start by activating the python virtual environment that you created before:

```bash
cd ~
cd ploidyNGS
source .venv/bin/activate
```

And then run:

```bash
./ploidyNGS.py -o diploidTest -b test_data/simulatedDiploidGenome/Ploidy2.bowtie2.sorted.bam
```

This should print the following in your screen:

```bash
No index available for pileup. Creating an index...
Getting the number of mapped reads from BAM
460400
```

And generate the files:

* diploidTest_depth100.tab
* diploidTest_depth100.tab.PloidyNGS.pdf

The PDF file should have a histogram identical to this https://github.com/diriano/ploidyNGS/tree/master/images/diploidTest_depth100.tab.PloidyNGS.png

After running ploidyNGS do not forget to deactivate your python virtual environment:

```bash
deactivate
```
## ploidyNGS usage and examples

Active the python virtual environment, created above, before using `ploidyNGS`:

````bash
cd ~/ploidyNGS/
source .venv/bin/activate
```

The switch `-h` will give you access to the full help page:

```bash
./ploidyNGS.py -h
usage: ploidyNGS.py [-h] [-v] -o file -b mappingGenome.bam [-m 0.95 default)]
                    [-d 100 (default]

ploidyNGS: Visual exploration of ploidy levels

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -o file, --out file   Base name for a TAB file that will keep the allele
                        counts
  -b mappingGenome.bam, --bam mappingGenome.bam
                        BAM file used to get allele frequencies
  -m 0.95 (default), --max_allele_freq 0.95 (default)
                        Fraction of the maximum allele frequency (float betwen
                        0 and 1, default: 0.95)
  -d 100 (default), --max_depth 100 (default)
                        Max number of reads kepth at each position in the
                        reference genome (integer, default: 100)
```

There are two required parameters: the input BAM file (`-b` or `--bam`) and a string (`-o` or `--out`) that will be used to generate the output files.

Two additional parameters help you control some of the behavior of `ploidyNGS`:

* `-m` or `--max_allele_freq`: Float between 0 and 1. Default 0.95. Ignore positions where the frequencey of the most abundant allele is higher or equal than `max_allele_freq`
* `-d` or `--max_depth`: Integer. Default 100. Maximum sequencing depth to consider, e.g., if d=100, then only the first 100 mapped reads will be examined.

## Examples and test
-------------------

We have prepared four test datasets, representing ploidy levels from 1 to 4, that you can use to try `ploidyNGS`. These data are available in the folder test_data:

* HaploidGenome
* simulatedDiploidGenome
* simulatedTetraploidGenome
* simulatedTriploidGenome

### Running ploidyNGS main script (`ploidyNGS.py`)

`ploidyNGS.py` only requires a BAM and a base name for the output files. We will show you how to use it using two of the pre-generated test datasets:

First using the HaploidGenome data:

```bash
cd ~/ploidyNGS
source .venv/bin/activate
mkdir myTest
./ploidyNGS.py --out myTest/DataTestPloidy1.tab --bam test_data/HaploidGenome/Ploidy1.bowtie2.sorted.bam 
deactivate
```

The script outputs two files:
- `DataTestPloidy2.tab_depth100.tab`, the table with the frequencies of allele proportions in heterozygous sites. It has the four columns: Sequence name, Position, Type of allele (First, Second, Thirds and Fourth) ordered by abundance, and the frequencing of the allele.
- `DataTestPloidy2.tab_depth100.tab.PloidyNGS.pdf`, the histogram of allele frequencies. This is your main tool to decide on the ploidy level of yur organism.

Compare your file DataTestPloidy1.tab_depth100.tab.PloidyNGS.pdf it should be pretty similar to test_data/ploidyNGS_results/DataTestPloidy1.tab.ExplorePloidy.png

What you should see in that histogram is one peak around 95% for the series labelled First, which is the most common allele, this means that the most frequent was represented by 95% or more of the reads in most positions, which is compatible with a haploid genome. The series labelled Second, should have a peak around 5%, this means that in hereromorphic positions, the second most frequent allele have around 5% of the reads, these most likely represent sequencing errors.

Now, let's have a look into a diploid dataset:

```bash
cd ~/ploidyNGS
source .venv/bin/activate
./ploidyNGS.py --out myTest/DataTestPloidy2.tab --bam test_data/simulatedDiploidGenome/Ploidy2.bowtie2.sorted.bam
deactivate
```

Compare your file DataTestPloidy2.tab_depth100.tab.PloidyNGS.pdf it should be pretty similar to test_data/ploidyNGS_results/DataTestPloidy2.tab.ExplorePloidy.png

You should see two peaks in the histogram, both located around 50%. This means that both the most frequent and second most frequence alleles, each, have around 50% of the reads in the heteromorphic positions. This result is compatible with a diploid genome.

Try running ploidyNGS for the simulated triploid and tetraploid datasets. Read the following section for more details on how to interpret the histogram generated by `ploidyNGS`.

### General guidelines to interpret the results of ploidyNGS for different ploidy levels

Additionally you will also find the folder ploidyNGS_results that contains the graphs generated by `ploidyNGS` foreach of these datasets.
The resulting PDF (PNG versions are also provided) files for each ploidyNGS run are:

- `DataTestPloidy1.tab.ExplorePloidy.pdf` - haploid organism
- `DataTestPloidy2.tab.ExplorePloidy.pdf` - diploid
- `DataTestPloidy3.tab.ExplorePloidy.pdf` - triploid
- `DataTestPloidy4.tab.ExplorePloidy.pdf` - tetraploid

For the interpretation of the resulting plots you can either use our simulation data (https://github.com/diriano/ploidyNGS/tree/master/simulation/results) or the following table. Briefly, in `ploidyNGS`'s graphs, the X-axis position of the peaks is directly related to the ploidy level of the organism under study. For instance, if you observe peaks for the two most frequent alleles close to 50%, you have a diploid organism (see `DataTestPloidy2.tab.ExplorePloidy.pdf`). The percentages mean that for heteromorphic positions, each allele was covered by approx. that percent of total reads in the positions.


| Ploidy level | Second most frequent allele (expected proportion) | Most frequence allele (expected proportion) |
| ------------ |-------------------------------:|-------------------------------:|
| Diploidy | 50 | 50 |
| Triploidy | 33.33 | 66.67 |
| Tetraploidy | 25 | 75 |
| Tetraploidy | 50 | 50 |
| Pentaploidy | 20 | 80 |
| Pentaploidy | 40 | 60 |
| Hexaploidy | 16.67 | 83.33 |
| Hexaploidy | 50 | 50 |
| Hexaploidy | 33.33 | 66.67 |
| Heptaploidy | 28.57 | 71.43 |
| Heptaploidy | 42.86 | 57.14 |
| Heptaploidy | 14.29 | 85.71 |

## Full analysis example - diploid organism

Here is a complete example of a pipeline for ploidy analysis using `ploidyNGS`.
For this example, we use a simulated diploid yeast chromosome generated using our script `simulatePloidyData.py`, which results in the plot above, `test_data/ploidyNGS_results/DataTestPloidy2.tab.ExplorePloidy.pdf`.

# RECOMMENDATIONS

`ploidyNGS` was developed to assess ploidy level in small genomes, of a few tens Mbp. If your genome is larger than that we recommend that you:

a) look at individual chromosomes. For this you can use bamtools, check: https://www.biostars.org/p/46327/.

b) look at a region of a single chromosome.

The running time and memory usage of `ploidyNGS` depends on the number of reads in the BAM file, i.e., the sequencing depth. We have seen that a sequencing depth of 100x is enough to get a good idea of the ploidy level. So if you have sequenced your genome to a larger depth, please sub-sample you data before creating the BAM file for `ploidyNGS`. Alternatively you could use the parameter `-d` in `ploidyNGS` to look at the `d` first mapped reads, however pysam will still load the full BAM file.

# REFERENCES

Huang, Weichun, et al. "ART: a next-generation sequencing read simulator." *Bioinformatics* 28.4 (2012): 593-594.

Li, Heng, et al. "The sequence alignment/map format and SAMtools." *Bioinformatics* 25.16 (2009): 2078-2079.
