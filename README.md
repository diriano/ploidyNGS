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

If you do not want, or you cannot, install these two modules globally, you could run:

```bash
pip install --user --upgrade pip
pip install --user virtualenv
```

Go to the ploidyNGS folder, create and activate a new virtual environment

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

```bash
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
                        This parameter will be automatically set based on the
                        coverage of your BAM file, however you can overrride
                        it.
```

There are two required parameters: the input BAM file (`-b` or `--bam`) and a string (`-o` or `--out`) that will be used to generate the output files.

Additional parameters help you control some of the behavior of `ploidyNGS`:

* `-m` or `--max_allele_freq`: Float between 0 and 1. Default 0.95. Ignore positions where the frequencey of the most abundant allele is higher or equal than `max_allele_freq`
* `-d` or `--max_depth`: Integer. Default 100. Maximum sequencing depth to consider, e.g., if d=100, then only the first 100 mapped reads will be examined.
* `-g` pr `--guess_ploidy`: Boolean. If present, `ploidyNGS` will try to guess the ploidy level of your data by comparing your allele frequency distribution to our simulated data, and computing the Kolmogorov-Smirnov distance, this will generate a table with the suffix: .ks-distance.PloidyNGS.tbl, also `ploidyNGS` will report on the screen the ploidy of the most similar distribution in the simulated data.
* `-c` or `--coverage_guess_ploidy`: One of 15, 25, 50, 100. When you select the option --guess_ploidy your sample will be compared with simulated datasets of a given coverage. `ploidyNGS.py` will automatically select which coverage value to use base on the actual observed average coverage in your sample. You can override this behaviour by using --coverage_guess_ploidy and one of the values 15, 25, 50 or 100, ie., the 4 subsets of simulated data that we provide, that correspond to different sequencing depths/coverages.

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
## This is ploidyNGS version v3.1.0
## Current date and time: Thu Mar  9 13:09:55 2017
###############################################################
BAM index present... OK!
Number of mapped reads from BAM: 206062
Observed average coverage: 65.00
Coverage used for guessing ploidy: 50

  After comparing your data with our simulated dataset
  and computing the Kolmogorov-Smirnov distance, 
  the closest ploidy to yours is 1
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

#### Important note on ploidy guessing with ploidyNGS

The original idea with `ploidyNGS` was to have a simple, model-free, visualization of the frequency of alternative alleles on a newly sequenced genome, using short reads next generation sequencing. Other software had tried to implement proper ploidy level estimation, but this usually requires complex models, that could be very sensitive to deviation of their assumptions (See for instance: Bao et al., 2014; Margarido & Heckerman, 2015; Yu et al., 2014). One of our main goals was to provide a model-free visualization of the polymorphisms found in a newly sequenced genome, and let the user decide about the ploidy level, based on visual comparisons with simulated data also provided by `ploidyNGS`. We understand that inferring ploidy level could be a valuable feature for end-users, and we have implemented a simple routine to compare the profile generated from the user’s data to the sets of simulated profile distributions for different ploidy and heterozygosity levels and sequencing coverages. In this way `ploidyNGS` can try to guess the ploidy level of the user provided data. However, **there are several factors impacting on the polymorphism profile distribution, and the ultimate decision about ploidy level relies on the inspection of the `ploidyNGS` generated graph by the user, and not `ploidyNGS`'s guess**, this is highlighted in the output of `ploidyNGS`. The routine that we have implemented in ploidyNGS computes the Kolmogorov-Smirnov distance between pairs of distributions (user vs simulated) and suggests the ploidy of the distribution with the smallest distance.

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

When using the parameter --guess_ploidy, ``ploidyNGS.py` will try to guess the ploidy in your sample by comparing the pattern in which first and second most frequent alleles appear to the same patterns of simulated data of known ploidy. For this we have generated simulated data at different coverages, i.e., 15, 25, 50 and 100. `ploidyNGS.py` will select the closest coverage to yours, in the following way: 
* If the average coverage in your sample is smaller than 20, your data will be compared against datasets of coverage 15. 
* If you coverage is >=20 and < 37.5, then it will compared against datasets of coverage 25. 
* If you coverage is >=37.5 and < 75, then it will compared against datasets of coverage 50. 
* If you coverage is >=75, then it will compared against datasets of coverage 100. 

The Kolmogorov-Smirnov will be computed between pairs of datasets, i.e., you dataset vs one of the simulated dataset of a given coverage and ploidy. The comparison with the smallest distance will be reported as the possible ploidy of your sample. We advise that you always look at the PDF generated by `ploidyNGS.py`, and do not rely only on the guess that we provide. If you find any inconsistencies please contact us.


# RECOMMENDATIONS

`ploidyNGS` was developed to assess ploidy level in small genomes, of a few tens Mbp. If your genome is larger than that we recommend that you:

a) look at individual chromosomes. For this you can use bamtools, check: https://www.biostars.org/p/46327/.

b) look at a region of a single chromosome.

c) Look at regions of at least 1Mbp.

The running time and memory usage of `ploidyNGS` depends on the number of reads in the BAM file, i.e., the sequencing depth. We have seen that a sequencing depth of 100x is enough to get a good idea of the ploidy level. So if you have sequenced your genome to a larger depth, please sub-sample you data before creating the BAM file for `ploidyNGS`. Alternatively you could use the parameter `-d` in `ploidyNGS` to look at the `d` first mapped reads, however pysam will still load the full BAM file.

# REFERENCES

Margarido, G. R. A. and Heckerman, D. (2015). ConPADE: Genome assembly ploidy estimation from next-generation sequencing data. PLoS Comput. Biol., 11(4), e1004229.

Bao, L., Pu, M., and Messer, K. (2014). AbsCN-seq: a statistical method to estimate tumor purity, ploidy and absolute copy numbers from nextgeneration sequencing data. Bioinformatics, 30(8), 1056–1063.

Yu, Z., Liu, Y., Shen, Y.,Wang, M., and Li, A. (2014). CLImAT: accurate detection of copy number alteration and loss of heterozygosity in impure and aneuploid tumor samples using whole-genome sequencing data. Bioinformatics, 30(18), 2576–2583.
