---
title: Variant Call
permalink: variant_call.html
sidebar: generic
tags: [PopGen1]
#product: Generic
---

In order to prepare the data for the Admixture analysis, we need to get all variants in the genomes of each population.

Files:
	* Variants: a vcf file (Variant Call Format)
	* Meta data: a csv file

Find more information of the different file types required at [Samtools](https://samtools.github.io/hts-specs/).

### VCF file:

Samtools define a _vcf_ file as a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome. The format also has the ability to contain genotype information on samples for each position.

### Create a VCF file (Platypus):

We need a reference genome in order to generate the variant call. Also a haplotype file in order to generate it.

We will use the tool _Platypus_ from the University of Oxford ([Platypus](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/platypus-a-haplotype-based-variant-caller-for-next-generation-sequence-data)). Platypus requieres a list of BAM files and the reference, both must be properly indexed. Samtools is recommended for that purpose. BAM (Binary Alignment/Map format) files are a binary version of SAM files. SAM files (Sequence Alignment/Map format)are a tab-delimited text file containing sequence alignment data.

In other words, the steps to generate a VCF file from scratch are:

<p>&nbsp;</p>

	* Get sequences to align and reference sequence for aligment.
	* Generate the SAM file with the alignment information of sequences against reference.
	* Transform SAM to BAM file.
	* Call the variants using the BAM file/s and the reference.

<p>&nbsp;</p>

For the second step, once we have the data, we will use _Samtools_. The source code of the software can be download from [Samtools](http://www.htslib.org/download/). The website shows how to build and install the software. For a pre-build software, the one we will use, we will download the packages from [Samtools Github page](https://github.com/samtools). More exactly, we will clone the repositories _htslib_ and _samtools_.

<p>&nbsp;</p>

In a linux environment:
```
mkdir /opt/tools
cd /opt/tools
#git clone https://github.com/samtools/htslib.git
git clone https://github.com/samtools/samtools.git
export PATH=/opt/tools/samtools/bin:$PATH
```
<p>&nbsp;</p>

In order to download _Platypus_ we can follow the instructions in the [main website](https://www.rdm.ox.ac.uk/research/lunter-group/lunter-group/platypus-a-haplotype-based-variant-caller-for-next-generation-sequence-data). In the section _Quick links_ we can find the link to download the latest stable version. After downloading the file and into a linux enviroment, I install _Platypus_.

<p>&nbsp;</p>

```
tar -xvzf platypus-latest.tgz
cd Platypus_0.8.1
```
<p>&nbsp;</p>

Basic usage:

	python Platypus.py callVariants --bamFiles=data.bam --refFile=ref.fa --output=out.vcf





**Install BiocManager in R**

<p>&nbsp;</p>

{% highlight R %}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install() 	# Push an update
{% endhighlight %}

<p>&nbsp;</p>

After having the package manager for biological libraries, BiocManager, we can install the packages required for our analysis.

<p>&nbsp;</p>

{% highlight R %}
BiocManager::install("SNPRelate")
{% endhighlight %}

<p>&nbsp;</p>

This command should install the library _SNPRelate_ and its dependecies.



DOWNLOAD:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/