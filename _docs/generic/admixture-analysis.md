---
title: Admixture
permalink: admixture-analysis.html
sidebar: generic
#tags: [PopGen1]
#product: Generic
---

## Data

In this section, we are going to use the [1000 genomes](https://www.internationalgenome.org/) data maintained by EMBL-EBI. As writen in the [NIH](https://www.genome.gov/27528684/1000-genomes-project) (National Human Genome Research Institute).

```
The 1000 Genomes Project is a collaboration among research groups in the US, UK, and China and Germany to produce an extensive catalog of human genetic variation that will support future medical research studies. It will extend the data from the International HapMap Project [...] The genomes of over 1000 unidentified individuals from around the world will be sequenced using next generation sequencing technologies. The results of the study will be publicly accessible to researchers worldwide.
```

The data use in this tutorial is a ***modified VCF file*** of the chromosome 22. The data was downloaded from the [ftp site](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) of 1000genomes and reduced to a smaller dataset to speed up the computation as much as possible. In case you have a bigger or better computer (e.g., a cluster) you can select and analyze the chromosome you prefer from the mention ftp site (be aware of adapting your code to your neccesities). Both the annotation of the samples, in a ordered csv file, and the modified _VCF_ file can be download it from [here](https://github.com/adriangeerre/popgen.github.io/tree/master/analysis/admixture). The annotation is the same one as for the chrosome 22 complete _VCF_ file but it does not affect our analisis.

## File format

The ***Variant Call Format*** or _VCF_ is defined by ***Samtools*** as a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome. The format also has the ability to contain genotype information on samples for each position. I need to add that the _VCF_ file is created from a _BAM_ file which is created from a _SAM_ file. The _SAM_ file is created by aligning the genome of interest against the reference genome, in this case an individual genome from a person inside of a population against the human reference genome. For practical purposes, we are not doing all the process but directly using the _VCF_. The _fastq_ sequences and _BAM_ files can be found in the [ftp site](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/). 

## Admixture

The definition of the word _admixture_ is "A mixture" or "something mixed with something else". In genetics, the word references to the mixture of individuals, indeed, the mixture of the genetic information, from different groups. We understand the mixture as the mating and the generation of descendants. Basically, the estimate of admixture between differentiated populations is done to determine how mixed they are and characterize their evolutionary history.  For that purpose, we use the variants/SNPs in the populations to compare their similarity.

The admixture can be estimated using variants in respect to a reference genome. In this section, we will use the software [ADMIXTURE](https://dalexander.github.io/admixture/download.html) (which used to be available in the _University of California_). Other packages are able to estimate admixture within a population, for example, _ipADMIXTURE_ (2020) which uses Q matrices. Also, inside the software _ADMIXTURE_ of the last program you can get sample data files from the [HapMap3](https://www.sanger.ac.uk/resources/downloads/human/hapmap3.html). However, I find more insteresting to use [1000 genomes](https://www.internationalgenome.org/) since [NCBI](https://www.ncbi.nlm.nih.gov/variation/news/NCBI_retiring_HapMap/) took HapMap 3 down and it is more recent.

The [manual](https://dalexander.github.io/admixture/admixture-manual.pdf) of _ADMIXTURE_ explains how to run the program. The software can work with three different inputs: **.bed**, **.ped** and **.geno** formats. The output is a space-delimited file know as the _Q estimates_. Also, we need to define the number of population expected (_K_). In our case, we know there are 26 different human-based populations but it the genome may be divided in a different number.

<p>&nbsp;</p>

**Choosing K**

K is the number of populations inside our dataset and it is defined by the value of the cross-validation. A low cross-validation resembles the best K for our dataset. An easy way to define K is by running the model for the different values and compare their error. As the manual says:

{% highlight Bash %}
for K in <range>; do admixture --cv <input> $K | tee log$K.out; done
{% endhighlight %}

Seems easy, right? Run a For loop and select the minimum error value. But, what if we do not know the maximum an minimum values? For example, our data is made out of 26 populations potentially divided in 5 regions, maybe the real number of populations is 3 but it could also be 26. We can not use intuition to define the range but we cannot run 26 models for masive data in a common computer (also, computation in a cluster cost money!!). In order to get a hint on the possible range we will use **Principal Components Analysis** (PCA).

<p>&nbsp;</p>

**Install BiocManager in R**

{% highlight R %}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install() 	# Push an update
{% endhighlight %}

After having the package manager for biological libraries, BiocManager, we can install the packages required for our analysis.

{% highlight R %}
BiocManager::install("SNPRelate")
{% endhighlight %}

This command should install the library SNPRelate and its dependecies.

<p>&nbsp;</p>


**PCA**

Once we have the package _SNPRelate_ installed in R, we can start analysing our genotype data of the chromosome 22. The first step is to load the libraries and transform our _VCF_ file into a binary version of it, a _GDS_ file. This process only needs to be run one time because it will produce a file that we will use to compute the PCA.

{% highlight R %}
# Load libraries
library(tidyverse)
library(SNPRelate)

# Set working directory
setwd("<Working Directory with the csv and vcf/vcf.gz file>")

# Transform VCF to GDS (Run one time!!)
vcf_file <- "chr22.phase3.vcf.gz"
snpgdsVCF2GDS(vcf_file, "chr22.phase3.gds", method="biallelic.only")
{% endhighlight %}

Following, we will load the metadata file and the recently created _GDS_ file. After that, we will run a first PCA and check the content of the output. 

{% highlight R %}
# Load data
metadata <- read.csv("igsr_samples_phase3.csv", sep = '\t') # Metadata
variants <- snpgdsOpen("chr22.phase3.gds") # Open gds file (Variants)

# PCA
pca <- snpgdsPCA(variants)
summary(pca)
{% endhighlight %}

The output of the PCA is a list that contains the eigenvectors (or PCs) and eigenvalues, also, the id of the sample and the variance per PC. In order to plot the PCs, first, we need to check the variance explain by each PC individually. The easies way to check the variance is to plot them.

{% highlight R %}
# PCA variance
pca_var <- pca$varprop[!is.nan(pca$varprop)] * 100 # It contains NAN because the variance is below the representable limit
pcs <- seq(1, length(pca_var))
plot(pcs, pca_var, type = "b", col = "red", pch = c(16))
{% endhighlight %}

<a href="http://adriangeerre.github.io/popgen.github.io/analysis/admixture/images/PCA_variance.png">
	<img src="http://adriangeerre.github.io/popgen.github.io/analysis/admixture/images/PCA_variance.png" alt="PCA variance" style="width:100%">
</a>

The basic plot shows that the first two PCs contains the most variance (5.49%) while the rest contains a small percent of the variance. Given the results, we will plot the two first PCs. First, we match the sample id with the population and region. Then, we plot the two first PCs colored by region.

<a href="http://adriangeerre.github.io/popgen.github.io/analysis/admixture/images/PCA.png">
	<img src="http://adriangeerre.github.io/popgen.github.io/analysis/admixture/images/PCA.png" alt="PCA" style="width:100%">
</a>

In order to run _ADMIXTURE_ we need to generate the correct input file. In this case we are going to use a _.bed_ file and the supporting file (_bim_ and _fam_). We can use the following code to generate the file in R.

{% highlight R %}
# GDS to BED to use in ADMIXTURE
snpgdsGDS2BED(variants, "chr22.phase3.gds", sample.id = NULL)
{% endhighlight %}

<p>&nbsp;</p>


**Running Admixture**

To compute the admixture between the population I followed the recomendation of the manual and compute a range of different K values. We expect to have 4/5 populations, meaning K equal to 4/5, were the population for America is not independent but a mixture of all the different populations. 

{% highlight Bash %}
for K in 1 2 3 4 5 6 7 8 9 10; do admixture_linux-1.3.0/admixture --cv chr22.phase3.gds.bed $K | tee log$K.out; done
{% endhighlight %}

After running the admixture for the different values of _K_ (it tooks a 3/4 hours in my computer), here are out results of the cross-validation:

```
cat log*.out | grep CV
CV error (K=1): 0.07605
CV error (K=2): 0.06851
CV error (K=3): 0.06639
CV error (K=4): 0.06621
CV error (K=5): 0.06569
CV error (K=6): 0.06631
CV error (K=7): 0.06770
CV error (K=8): 0.06877
```

As expected by the PCA, the lowest error (0.06569) is given for K equal to 5 (5 populations). However, this value is close to K equal to 4 or 6 but it may be an artifact of the small sample used for the tutorial. The values of _Fst_ for K equal to 5 are:

```
Fst divergences between estimated populations: 

	    Pop0	Pop1	Pop2	Pop3	
Pop0	
Pop1	0.120	
Pop2	0.125	0.191	
Pop3	0.048	0.137	0.129	
Pop4	0.080	0.131	0.161	0.110	
```

The _F-statistic_ or _Fst_, is the proportion of the total genetic variance contained in a subpopulation relative to the total genetic variance. Values range from 0 to 1. High _Fst_ implies a considerable degree of differentiation among populations. It can be estimated with:

<p>&nbsp;</p>

<center><img src="https://latex.codecogs.com/svg.latex?F_{st}=\frac{H_T - H_S}{H_T}" title="F-statistic estimation"/></center>

<p>&nbsp;</p>

where:

<p>&nbsp;</p>

<center><img src="https://latex.codecogs.com/svg.latex?H_T=\frac{n_T}{n_T - 1}(1 - \sum \bar{p}^2_i)&space;\quad&space;H_S=\frac{n_S}{n_S - 1}(1 - \sum p^2_i)" title="Allele frequency (di-allelic model)"/></center>

<p>&nbsp;</p>

Given the formulas, the first term _Ht_ is the expected heterozygosity in the total sample population, and _Hs_ is the average of the expected heterozygosity calculated from each sampled subpopulation (_nt_ is the sample size and _ns_ is the own sample size). Basically, the T stands for _total population_ and the S for _subpopulation_.

Then, the _Fst_ matrix from the admixture of _K_ equal to 5 shows that the Pop0 and Pop3 are the closest followed by Pop0 and Pop4. On the other hand, the populations with the largest difference are Pop1 and Pop2. However, the values are close to 0 so we do not see a large diferentiation between populations.

<p>&nbsp;</p>


**Plotting Admixture**

After running the admixture software and select the value of K. We can load the Q matrix in R and plot the results (as the Admixture manual defines). The Q-matrix I used are available [here](https://github.com/adriangeerre/popgen.github.io/tree/master/analysis/admixture).

{% highlight R %}
# Q-matrix
q_matrix <- read.table("chr22.phase3.gds.5.Q")
ord <- q_matrix[order(q_matrix$V1, q_matrix$V2, q_matrix$V3, q_matrix$V4, q_matrix$V5),] # Order min to max value
admixture_plot <- barplot(t(as.matrix(ord)), space=c(0.2), col=rainbow(5), xlab="Individual #", ylab="Ancestry", border=NA, las=2)
{% endhighlight %}

<a href="http://adriangeerre.github.io/popgen.github.io/analysis/admixture/images/Admixture_plot.png">
	<img src="http://adriangeerre.github.io/popgen.github.io/analysis/admixture/images/Admixture_plot.png" alt="Admixture for K=5" style="width:100%">
</a>

The admixture plot shows 5 different sections. The first three sections (purple, blue and green) contain almost pure individuals. The last two (yellow and red) contain a remarkable amounf of individuals with variants from another regions. The yellow region, include the most admixture from other regions. The different regions cannot be defined from the admixture plot but we can expect that the yellow region is America, the red Africa and the other three the rest, similar to the PCA clustering.

Finally, the package [starmie](https://github.com/sa-lee/starmie), from 2016, plot admixture with advanced options (similar to tidyverse).

<p>&nbsp;</p>
