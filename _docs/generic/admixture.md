---
title: Admixture
permalink: admixture.html
sidebar: generic
#tags: [PopGen1]
#product: Generic
---

## Data

In this section, we are going to use the [1000 genomes](https://www.internationalgenome.org/) data maintained by EMBL-EBI. As writen in the [NIH](https://www.genome.gov/27528684/1000-genomes-project) (National Human Genome Research Institute).

<p>&nbsp;</p>

```
The 1000 Genomes Project is a collaboration among research groups in the US, UK, and China and Germany to produce an extensive catalog of human genetic variation that will support future medical research studies. It will extend the data from the International HapMap Project [...] The genomes of over 1000 unidentified individuals from around the world will be sequenced using next generation sequencing technologies. The results of the study will be publicly accessible to researchers worldwide.
```

<p>&nbsp;</p>

We will use the dataset ***Phase 3*** (released in 2013) which consists of 26 human populations with a total 2504 individuals. The data can be downloaded in _VCF_ format, it is divided by chromosomes and stored in the [ftp site](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). For didactical purpose, we will use the smallest human autosomal chromosome which is the **chr22** (files **_vcf.gz_** and **_vcf.gz.tbi_**). Please, acces the [ftp site](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) and download the data. In case you have a bigger or better computer (e.g., a cluster) you can select and analyze the chromosome you prefer (be aware of adapting your code to your neccesities).

The ***Variant Call Format*** or _VCF_ is defined by ***Samtools*** as a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome. The format also has the ability to contain genotype information on samples for each position. I need to add that the _VCF_ file is created from a _BAM_ file which is created from a _SAM_ file. The _SAM_ file is created by aligning the genome of interest against the reference genome, in this case an individual genome from a person inside of a population against the human reference genome. For practical purposes, we are not doing all the process but directly using the _VCF_. The _fastq_ sequences and _BAM_ files can be found in the [ftp site](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).

## Admixture

***Include definition of admixture***

The admixture can be estimated using variants in respect to a reference genome. In this section, we will use the software [ADMIXTURE](https://dalexander.github.io/admixture/download.html) (which used to be available in the _University of California_). Other packages are able to estimate admixture within a population, for example, _ipADMIXTURE_ (2020) which uses Q matrices. Also, inside the software _ADMIXTURE_ of the last program you can get sample data files from the [HapMap3](https://www.sanger.ac.uk/resources/downloads/human/hapmap3.html). However, I find more insteresting to use [1000 genomes](https://www.internationalgenome.org/) since [NCBI](https://www.ncbi.nlm.nih.gov/variation/news/NCBI_retiring_HapMap/) took HapMap 3 down and it is more recent.

The [manual](https://dalexander.github.io/admixture/admixture-manual.pdf) of _ADMIXTURE_ explains how to run the program. The software can work with three different inputs: **.bed**, **.ped** and **.geno** formats. The output is a space-delimited file know as the _Q estimates_. Also, we need to define the number of population expected (_K_). In our case, we know there are 26 different human-based populations but it the genome may be divided in a different number.



**Choosing K**

K is the number of populations inside our dataset and it is defined by the value of the cross-validation. A low cross-validation resembles the best K for our dataset. An easy way to define K is by running the model for the different values and compare their error. As the manual says:

<p>&nbsp;</p>

{% highlight Bash %}
for K in <range>; do admixture --cv <input> $K | tee log$K.out; done
{% endhighlight %}

<p>&nbsp;</p>

Seems easy, right? Run a For loop and select the minimum error value. But, what if we do not know the maximum an minimum values? For example, our data is made out of 26 populations, maybe the real number of populations is 3 but it could also be 50. We can not use intuition to define the range but we cannot run 50 models for masive data. In order to get a hint on the possible range we will use **Principal Components Analysis** (PCA).

**Install BiocManager in R**

After having the package manager for biological libraries, BiocManager, we can install the packages required for our analysis.

This command should install the library SNPRelate and its dependecies.