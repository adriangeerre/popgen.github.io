---
title: Mapping reads
permalink: mapping_reads.html
sidebar: generic
#tags: [PopGen1]
#product: Generic
---

## Data

**Explanation**

Everybody knows about Covid-19, Coronavirus or SARS-CoV-2 (Severe acute respiratory syndrome coronavirus 2). The information of the consensus sequence can be found [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). The Wuhan virus (Accession: NC_045512) contains around 29903 bp in a ss-RNA sequence. The reference genome (the first isolate virus) can be found in [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). To download the reference DNA sequence click on _FASTA_ under the title. Once inside, you can click on _Send to_, select _Complete Record_, _File_ and format _FASTA_. Finally, click on create file and save as _SARS-CoV-2-reference.fasta_. All the information in NCBI for the virus can be access from [here](https://www.ncbi.nlm.nih.gov/sars-cov-2/).

<p>&nbsp;</p>

**Download data**

To get _Illumina_ data can be found in [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) searching for _SARS-CoV-2_. In this case, I choose the accesion [SRX9197062](https://www.ncbi.nlm.nih.gov/sra/SRX9197062[accn]) because of didactical purposes. The short genome of the virus allows a wider vision and an easier computation than, for example, the human genome. The following image shows a screenshot of the accesion:

<p>&nbsp;</p>

<img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/project_selected.png" alt="SARS-CoV-2 selected project" style="width:100%">

<p>&nbsp;</p>

The picture shows the information related to the accession. For example, the number of bases, the size of the file, the day it was publish and other data, such as, sequencing machine used. By clicking on the highlighted area of the picture you will access the run:

<p>&nbsp;</p>

<img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/run_of_project.png" alt="Project run" style="width:100%">

<p>&nbsp;</p>

The run contains different tabs (Metadata, Analysis, Reads and Data Access). We are interested in _Reads_. Inside of it, we can find each read individually and we can also download the data by clicking on the highlighted area.

<p>&nbsp;</p>

<img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/reads_of_run.png" alt="Project reads" style="width:100%">

<p>&nbsp;</p>

Finally, to download the data click on _Download_. The link is highlighted in the following picture:

<p>&nbsp;</p>

<img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/download_illumina.png" alt="Download Illumina reads" style="width:100%">

<p>&nbsp;</p>

A second option to download the data is to go to the tab _Data access_ and use the external links. For example, in this case the data is stored inside _Amazon S3_. Following the url we can direct download the data to our computer. In relation to the data, the size of the downloaded file, compress with gunzip (.gz) is of 204.6 MB. The decompress file weights 479 MB. The downloaded file contain a total of 6.479.048 Illumina reads where all of them have a length of 36 bp. 

<p>&nbsp;</p>

**Available data**

The data, together with the pictures, is available in the github of the website. Find the resources [here](https://github.com/adriangeerre/popgen.github.io/tree/master/analysis/mapping_reads).

The data was decompress and divide in smaller pieces to reduce the size of the fasta file and allow Github to storage it.

{% highlight Bash %}
gzip -d SARS-CoV-2_exper-SRX9197062.fasta.gz
wc -l SARS-CoV-2_exper-SRX9197062.fasta # Count the number of lines
split -a 2 -d -l 1000000 SARS-CoV-2_exper-SRX9197062.fasta 'SARS-CoV-2_exper-SRX9197062.fasta_part'
{% endhighlight %}

The output were 13 files (part00 to part12, all under 50 MB) containing 1 million lines (500.000 reads) each. In order to re-create the original file, download all the files and inside of a linux/Mac terminal run:

{% highlight Bash %}
cat SARS-CoV-2_exper-SRX9197062.fasta_part* > SARS-CoV-2_exper-SRX9197062.fasta
{% endhighlight %}

## Analysis

**Explanation**

.................

<p>&nbsp;</p>

**Software**

Download the software _BWA_ from [here](bio-bwa.sourceforge.net). Place the file in the folder you prefer and run the following to decompress:

{% highlight Bash %}
bzip2 -d bwa-<version>.tar.bz2
tar -xvf bwa-<version>.tar
{% endhighlight %}

That should create a folder called bwa-_version_ were the program is contained. Now, we should compile the program by running the following (extracted from [Github](https://github.com/lh3/bwa)):

{% highlight Bash %}
cd bwa-<version>
bwa make
{% endhighlight %}

If everything works, a executable file called bwa would be created. Then, we can add the program folder into the path to quick access bwa by running `export PATH=$PATH:<path>/<to>/<bwa>` (temporal) or modifying the path by defining it inside the _.bashrc_ file.

<p>&nbsp;</p>

**Mapping**

Once the program is installed, we need to do three steps:

	1. Index the reference fasta sequence.
	2. Find the coordinates of the input reads.
	3. Map the reads against the reference sequence.

For the first step we need to run the following:

{% highlight Bash %}
bwa index SARS-CoV-2-reference.fasta
{% endhighlight %}

This should takes a few seconds and should generate the following output:

```
SARS-CoV-2-reference.fasta.amb
SARS-CoV-2-reference.fasta.ann
SARS-CoV-2-reference.fasta.bwt
SARS-CoV-2-reference.fasta.pac
SARS-CoV-2-reference.fasta.sa
```

Once we have generated the database from the reference genome, we can compute the second step. Here we can define a few parameters in order to locate the reads, for example, number of gaps, extensions, deletions, among others. **We will keep all values in default** but I recommend to play around and check the different outputs. This process will generate a binary file with format _.sai_ that will be used in the final step, in my case it took 56 seconds.

{% highlight Bash %}
bwa aln SARS-CoV-2-reference.fasta SARS-CoV-2_exper-SRX9197062.fasta > SARS-CoV-2_exper-SRX9197062.sai
{% endhighlight %}

Once we have located coordinates of the input reads, we can compute the final step. In this case, our reads are **single-end** and smaller than 100 bp (average: 36bp). Given the manual, we should run the _BWA-backtrack_ software with the subcomand _samse_ to generate a _SAM_ (Sequence Alignment/Map) file output.

{% highlight Bash %}
bwa samse SARS-CoV-2-reference.fasta SARS-CoV-2_exper-SRX9197062.sai SARS-CoV-2_exper-SRX9197062.fasta
{% endhighlight %}

The output is a SAM file with around 919 MB weight and 6.5 million lines. It is recomendable to pipe the output into gzip to compress the file and avoid consuming extra storage (`bwa samse <files> | gzip -3`).

<p>&nbsp;</p>
