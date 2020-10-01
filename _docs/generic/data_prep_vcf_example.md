---
title: Data preparation: VCF example for Admixture
permalink: data_prep_vcf_example.html
sidebar: generic
tags: [PopGen1]
#product: Generic
---

## Data download

**Describe the data download**

## Genomic Data preparation

Decompress the download. I run the command `awk 'END{print FNR}'` to get the total number of lines. The total number of lines was 1103800. Using `grep '#' <file> | wc -l` I measured the number of lines for the information in the vcf file. A total of 253 lines, including the header, defined the file. Following, I split the file in information and genomic data.

{% highlight Bash %}
head -n 252 ALL.chr22.phase3.genotypes.vcf > header         # Index information
tail -n 1103548 ALL.chr22.phase3.genotypes.vcf > variants   # Genomic data with header
{% endhighlight %}

Then, I remove columns inside of the genomic data to reduce the load of computation. This done for teaching purposes. It is obvious that the amount of data affects the results (the more the better) but in this case we want to have a example file that can run easily (fast) in an average computer. To remove the columns (individuals in our sample) we first count the total number.

{% highlight Bash %}
awk 'NR==1{print NF}' variants
{% endhighlight %}

The total number of columns is 2513 but the 9 columns are information not related to the individuals. That give us a total of 2504 columns to select the data. I use a random generator to select randomly the individuals. Then, I select the columns and create final vcf file.

{% highlight Bash %}
pre_cols=`shuf -i 10-2513 -n 350 | sort -n`
cols=`echo $pre_cols | sed 's/ /,/g'`
cut -f 1-9,$cols variants > pre_vcf
cat header pre_vcf > chr22.phase3.vcf
{% endhighlight %}

We can always check that the number of columns is correct after creating the file _pre_vcf_. The test example is less heavier (2.3GB) than the original file (11GB). We can gzip the vcf file (`gzip chr22.phase3.vcf`) to reduce the weight (61MB). Finally, I remove the intermediate files (i.e., header, variants and pre_vcf).

## Metadata Preparation
