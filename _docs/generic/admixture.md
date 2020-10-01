---
title: Admixture
permalink: admixture.html
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

***Include definition of admixture***

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

<p>&nbsp;</p>

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

<p>&nbsp;</p>

Following, we will load the metadata file and the recently created _GDS_ file. After that, we will run a first PCA and check the content of the output. 

<p>&nbsp;</p>

{% highlight R %}
# Load data
metadata <- read.csv("igsr_samples_phase3.csv", sep = '\t') # Metadata
variants <- snpgdsOpen("chr22.phase3.gds") # Open gds file (Variants)

# PCA
pca <- snpgdsPCA(variants)
summary(pca)
{% endhighlight %}

<p>&nbsp;</p>

The output of the PCA is a list that contains the eigenvectors (or PCs) and eigenvalues, also, the id of the sample and the variance per PC. In order to plot the PCs, first, we need to check the variance explain by each PC individually.

<p>&nbsp;</p>

{% highlight R %}
# PCA variance
pca_var <- pca$varprop[!is.nan(pca$varprop)] * 100 # It contains NAN because the variance is below the representable limit
pcs <- seq(1, length(pca_var))
plot(pcs, pca_var, type = "b", col = "red", pch = c(16)) # The first two PCs contains the most variance (2.32%)
{% endhighlight %}

<p>&nbsp;</p>





