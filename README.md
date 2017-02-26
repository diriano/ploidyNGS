# ploidyNGS: Visually exploring ploidy with Next Generation Sequencing data

**Motivation:** Due to the decreasing costs of genome sequencing researchers are increasingly assessing the genome information of non-model organisms using Next Generation Sequencing technologies. An important question to answer when assessing a new genome sequence is: What is the ploidy level of the organism under study?

**Results:** We have developed `ploidyNGS`, a model-free, open source tool to visualize and explore ploidy levels in a newly sequenced genome, exploiting short read data. We tested `ploidyNGS` using both simulated and real NGS data of the model yeast *Saccharomyces cerevisiae*. `ploidyNGS` allows the identification of the ploidy level of a newly sequenced genome in a visual way. We have applied `ploidyNGS` to a wide range of different genome ploidy and heterozigosity levels, as well as to a range of sequencing depths.

**Availability and implementation:** `ploidyNGS` is available under the GNU General Public License (GPL) at https://github.com/diriano/ploidyNGS. `ploidyNGS` is implemented in Python and R.

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
###############################################################
## This is ploidyNGS version v3.0.0
## Current date and time: Sun Feb 26 12:26:05 2017
###############################################################
No index available for pileup. Creating an index...
Getting the number of mapped reads from BAM
460400
```

And should generate the files:

* diploidTest_depth100.tab
* diploidTest_depth100.tab.PloidyNGS.pdf

The PDF file should have a histogram identical to this https://github.com/diriano/ploidyNGS/tree/master/images/diploidTest_depth100.tab.PloidyNGS.png

After running `ploidyNGS` do not forget to deactivate your python virtual environment:

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
###############################################################
## This is ploidyNGS version v3.0.0
## Current date and time: Sun Feb 26 12:26:05 2017
###############################################################
usage: ploidyNGS.py [-h] [-v] -o file -b mappingGenome.bam [-m 0.95 default)]
                    [-d 100 (default] [-g] [-c {15,25,50,100}]

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
  -g, --guess_ploidy    Try to guess ploidy level by comparison with simulated
                        data
  -c {15,25,50,100}, --coverage_guess_ploidy {15,25,50,100}
                        If you selected --guess_ploidy, It should be one of
                        the coverage values for simulated data. It defaults to
                        100
```

There are two required parameters: the input BAM file (`-b` or `--bam`) and a string (`-o` or `--out`) that will be used to generate the output files.

Additional parameters help you control some of the behavior of `ploidyNGS`:

* `-m` or `--max_allele_freq`: Float between 0 and 1. Default 0.95. Ignore positions where the frequencey of the most abundant allele is higher or equal than `max_allele_freq`
* `-d` or `--max_depth`: Integer. Default 100. Maximum sequencing depth to consider, e.g., if d=100, then only the first 100 mapped reads will be examined.
* `-g` pr `--guess_ploidy`: Boolean. If present, `ploidyNGS` will try to guess the ploidy level of your data by comparing your allele frequency distribution to our simulated data, and computing the Kolmogorov-Smirnov distance, this will generate a table with the suffix: .ks-distance.PloidyNGS.tbl, also `ploidyNGS` will report on the screen the ploidy of the most similar distribution in the simulated data.
* `-c` or `--coverage_guess_ploidy`: One of 15, 25, 50, 100. If you selected `--guess_ploidy`, then you coudl use this argument. It selected the subset of simulated data to compare against in order to guess ploidy level. We provide 4 subsets of simulated data, the correspond to different sequencing depths/coverages.

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
./ploidyNGS.py --out myTest/DataTestPloidy1 --bam test_data/HaploidGenome/Ploidy1.bowtie2.sorted.bam 
deactivate
```

The script outputs two files:
- `DataTestPloidy2.tab_depth100.tab`, the table with the frequencies of allele proportions in heterozygous sites. It has the four columns: Sequence name, Position, Type of allele (First, Second, Thirds and Fourth) ordered by abundance, and the frequencing of the allele.
- `DataTestPloidy2.tab_depth100.tab.PloidyNGS.pdf`, the histogram of allele frequencies. This is your main tool to decide on the ploidy level of yur organism.

Compare your file DataTestPloidy1.tab_depth100.tab.PloidyNGS.pdf it should be pretty similar to test_data/ploidyNGS_results/DataTestPloidy1.tab.ExplorePloidy.png

What you should see in that histogram is one peak around 95% for the series labelled First, which is the most common allele, this means that the most frequent was represented by 95% or more of the reads in most positions, which is compatible with a haploid genome. The series labelled Second, should have a peak around 5%, this means that in hereromorphic positions, the second most frequent allele have around 5% of the reads, these most likely represent sequencing errors.

In case you want ploidyNGs to try to guess the ploidy level in your sample, you should use the argument --guess_ploidy, see:

```bash
cd ~/ploidyNGS
source .venv/bin/activate
mkdir myTest
./ploidyNGS.py --guess_ploidy --out myTest/DataTestPloidy1_guessPloidy --bam test_data/HaploidGenome/Ploidy1.bowtie2.sorted.bam
deactivate
```

This will print the following in your screen:

```bash
###############################################################
## This is ploidyNGS version v3.0.0
## Current date and time: Sun Feb 26 12:26:05 2017
###############################################################
BAM index present... OK!
Getting the number of mapped reads from BAM
475920
Average coverage: 99.00

  After comparing your data with our simulated dataset
  and computing the Kolmogorov-Smirnov distance, 
  the closest ploidy to yours is 4
```

You will also find the file `myTest/DataTestPloidy1_guessPloidy_depth100.ks-distance.PloidyNGS.tbl`, that has a table with the comparison of your allele frequency distribution to all the simulated datasets with 100x coverage (this is the default coverage for comparison, you can change it with the argument --coverage_guess_ploidy).

Now, let's have a look into a diploid dataset:

```bash
cd ~/ploidyNGS
source .venv/bin/activate
./ploidyNGS.py --out myTest/DataTestPloidy2.tab --bam test_data/simulatedDiploidGenome/Ploidy2.bowtie2.sorted.bam --guess_ploidy
deactivate
```

Compare your file DataTestPloidy2.tab_depth100.tab.PloidyNGS.pdf it should be pretty similar to test_data/ploidyNGS_results/DataTestPloidy2.tab.ExplorePloidy.png

You should see two peaks in the histogram, both located around 50%. This means that both the most frequent and second most frequence alleles, each, have around 50% of the reads in the heteromorphic positions. This result is compatible with a diploid genome.

Try running ploidyNGS for the simulated triploid and tetraploid datasets. Read the following section for more details on how to interpret the histogram generated by `ploidyNGS`.

### General guidelines to interpret the results of ploidyNGS for different ploidy levels

For the interpretation of the resulting plots you can either use our simulation data (https://github.com/diriano/ploidyNGS/blob/master/simulation/results/results_ploidy_simulations.html) or the following table. Briefly, in `ploidyNGS`'s histograms, the position of the peaks on the X-axis is directly related to the ploidy level of the organism under study. For instance, if you observe peaks for the two most frequent alleles close to 50%, you have a diploid organism (see `DataTestPloidy2.tab.ExplorePloidy.pdf`). The values on the X-axis represent the fraction of reads representing each of the 4 possible alleles.

| Ploidy level | Second most frequent allele (expected proportion) | Most frequence allele (expected proportion) |
| ------------ |-------------------------------:|-------------------------------:|
| Diploid | 50 | 50 |
| Triploid | 33.33 | 66.67 |
| Tetraploid | 25 | 75 |
| Tetraploid | 50 | 50 |
| Pentaploid | 20 | 80 |
| Pentaploid | 40 | 60 |
| Hexaploid | 16.67 | 83.33 |
| Hexaploid | 50 | 50 |
| Hexaploid | 33.33 | 66.67 |
| Heptaploid | 28.57 | 71.43 |
| Heptaploid | 42.86 | 57.14 |
| Heptaploid | 14.29 | 85.71 |

# RECOMMENDATIONS

`ploidyNGS` was developed to assess ploidy level in small genomes, of a few tens Mbp. If your genome is larger than that we recommend that you:

a) look at individual chromosomes. For this you can use bamtools, check: https://www.biostars.org/p/46327/.

b) look at a region of a single chromosome.

The running time and memory usage of `ploidyNGS` depends on the number of reads in the BAM file, i.e., the sequencing depth. We have seen that a sequencing depth of 100x is enough to get a good idea of the ploidy level. So if you have sequenced your genome to a larger depth, please sub-sample you data before creating the BAM file for `ploidyNGS`. Alternatively you could use the parameter `-d` in `ploidyNGS` to look at the `d` first mapped reads, however pysam will still load the full BAM file.
