---
title: Data for Admixture analysis
permalink: admixture-data.html
sidebar: generic
tags: [PopGen1]
#product: Generic
---

## Data download

We will use the dataset ***Phase 3*** (released in 2013) which consists of around 26 human populations, 2504 individuals divided, and 5 Superpopulations (Africa, America, East Asia, Europe, and South Asia). The data can be downloaded in _VCF_ format, it is divided by chromosomes and stored in the [ftp site](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). We will use the smallest human autosomal chromosome which is the **chr22** (files **_vcf.gz_** and **_vcf.gz.tbi_**). The data was downloaded from the [ftp site](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) of 1000genomes.

## Genomic Data preparation

Decompress the download. I run the command `awk 'END{print FNR}'` to get the total number of lines. The total number of lines was 1103800. Using `grep '#' <file> | wc -l` I measured the number of lines for the information in the vcf file. A total of 253 lines, including the header, defined the file. Following, I split the file into information and genomic data.

{% highlight Bash %}
head -n 252 ALL.chr22.phase3.genotypes.vcf > header         # Index information
tail -n 1103548 ALL.chr22.phase3.genotypes.vcf > variants   # Genomic data with header
{% endhighlight %}

Then, I remove columns inside of the genomic data to reduce the load of computation. This done for teaching purposes. It is obvious that the amount of data affects the results (the more the better) but in this case, we want to have an example file that can run easily (fast) on an average computer. To remove the columns (individuals in our sample) we first count the total number.

{% highlight Bash %}
awk 'NR==1{print NF}' variants
{% endhighlight %}

The total number of columns is 2513 but the 9 columns are information not related to the individuals. That gives us a total of 2504 columns to select the data. I use a random generator to select randomly the individuals. Then, I select the columns and create the final vcf file.

{% highlight Bash %}
pre_cols=`shuf -i 10-2513 -n 300 | sort -n`
cols=`echo $pre_cols | sed 's/ /,/g'`
cut -f 1-9,$cols variants > pre_vcf
cat header pre_vcf > chr22.phase3.vcf
{% endhighlight %}

We can always check that the number of columns is correct after creating the file _pre_vcf_. The test example is less heavy (2.3GB) than the original file (11GB). We can gzip the vcf file (`gzip chr22.phase3.vcf`) to reduce the weight (61MB). Finally, I remove the intermediate files (i.e., header, variants, and pre_vcf).

## Metadata Preparation

To download the metadata access this [link](https://www.internationalgenome.org/data-portal/sample). The link should show a table inside of 1000 genomes website. Click on _Filter by data collection_, select _1000 Genomes phase 3 release_ and download the data by clicking on _Download the list_. Once I have downloaded the file and renamed it to _igsr_samples_phase3.tsv_ with the metadata for the downloaded _VCF_ files, I modify it by hand using [libreoffice](https://libreoffice.org/). I select and reorder the columns to contain:

```
Biosample ID
Sample name
Sex
Population name
Superpopulation code
Superpopulation name
```

However, it is not necessary since the data is going to be used by columns. I decide to clean and tidy the data f so it is easier to understand.