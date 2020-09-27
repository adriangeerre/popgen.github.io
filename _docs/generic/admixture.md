---
title: Admixture
permalink: admixture.html
sidebar: generic
#tags: [PopGen1]
#product: Generic
---

### Data

In this section, we are going to use the [1000 genomes](https://www.internationalgenome.org/) data maintained by EMBL-EBI. As writen in the [NIH](https://www.genome.gov/27528684/1000-genomes-project) (National Human Genome Research Institute).

```
The 1000 Genomes Project is a collaboration among research groups in the US, UK, and China and Germany to produce an extensive catalog of human genetic variation that will support future medical research studies. It will extend the data from the International HapMap Project [...] The genomes of over 1000 unidentified individuals from around the world will be sequenced using next generation sequencing technologies. The results of the study will be publicly accessible to researchers worldwide.
```

We will use the dataset _Phase 3_ which consists of 26 human populations with a total 2504 individuals. The data can be downloaded in VCF format, it is divided by chromosomes and stored in the [ftp site](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). 
For didactical purpose, we will use the smallest human autosomal chromosome which is the **chromosome 22** (2 files: **_vcf.gz_** and **_vcf.gz.tbi_**). Please, acces the ftp site and download the data. In case you have a bigger or better computer (e.g., a cluster) you can select and analyze the chromosome you prefer (be aware of adapting your code to your neccesities).

The Variant Call Format or VCF is defined by _Samtools_ as a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome. The format also has the ability to contain genotype information on samples for each position. I need to add that the VCF file is created from a BAM file which is created from a SAM file. The SAM file is created by aligning the genome of interest against the reference genome, in this case an individual genome from a person inside of a population against the human reference genome. For practical purposes, we are not doing all the process but directly using the VCF.

### Admixture




### References:

Patterson et al., 2012
