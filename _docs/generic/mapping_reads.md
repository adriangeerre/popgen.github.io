---
title: Mapping reads
permalink: mapping_reads.html
sidebar: generic
#tags: [PopGen1]
#product: Generic
---

## Data

**Explanation**

Everybody knows about Covid-19, Coronavirus or SARS-CoV-2 (Severe acute respiratory syndrome coronavirus 2). The information of the consensus sequence can be found [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). The Wuhan virus (Accession: NC_045512) contains around 29903 bp in a ss-RNA sequence. The reference genome (the first isolate virus) can be found in [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). To download the reference DNA sequence click on _FASTA_ under the title. Once inside, you can click on _Send to_, select _Complete Record_, _File_ and format _FASTA_. Finally, click on create file and save as _SARS-CoV-2-reference.fasta_. Moreover, following the same procedure we can download the gene annotation for the reference genome. Click on _Send to_, select _Complete Record_, _File_ and format _GFF3_. Save the file with the name _SARS-CoV-2-reference.gff3_ All the information in NCBI for the virus can be access from [here](https://www.ncbi.nlm.nih.gov/sars-cov-2/).

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

In this tutorial we are going to **map** _Illumina_ reads against a reference genome (without a reference we will _de novo_ assembly the genome). We are using reads comming from a Polymerase Chain Reaction (PCR) so we will have several times the same reads and also, reads that overlaps between them. The reference genome will be the template that will help us to place each reads in the proper location. Once we have our sample genome we can compare and define the differences (variants). The output of the mapping is a _SAM_ (Sequence Alignment/Map) file, read more about the format [here](http://samtools.github.io/hts-specs/SAMv1.pdf).

<p>&nbsp;</p>

**Software**

The software we are going to use in the first section is ***BWA***. There are many available softwares for mapping reads, for example, TopHat, MAQ or Bowtie. This [article](https://academic.oup.com/bioinformatics/article/28/24/3169/245777) list a large number of them. Different software may have different qualities or specializations depending the input. I do not have any especial interest in using _BWA_ for the tutorial, as long as it has a fast algorithm.

Download the software ***BWA*** from [here](bio-bwa.sourceforge.net). Place the file in the folder you prefer and run the following to decompress:

{% highlight Bash %}
bzip2 -d bwa-<version>.tar.bz2
tar -xvf bwa-<version>.tar
{% endhighlight %}

That should create a folder called bwa-_version_ were the files to compile the program are contained. Now, we should compile the program by running the following (extracted from [Github](https://github.com/lh3/bwa)):

{% highlight Bash %}
cd bwa-<version>
bwa make
{% endhighlight %}

If everything works, a executable file called bwa would be created. Then, we can add the program folder into the path to quick access bwa by running `export PATH=$PATH:<path>/<to>/<bwa>` (temporal) or modifying the path by defining it inside the _.bashrc_ file.

The software we are going to use in the second section is ***Picard***. The software was develop by the Broad Institute and has an open-source license. _Picard_ contains a huge number of tools to manipulate high-throughput sequencing (HTS) data. In order to use the software we need _java-1.8_ installed. Read more and download from [here](https://broadinstitute.github.io/picard/) by clicking the first box of the right upper corner (Latest Jar Release).

The software we are going to use in the thirs section are ***Samtools*** and ***IGV***. Download the _Samtools_ from [here](http://www.htslib.org/). To install the software, move it to the folder of your interest and run: 

{% highlight Bash %}
bzip2 -d samtools-<version>.tar.bz2
tar -xvf samtools-<version>.tar
rm samtool-<version>.tar
cd samtools-<version>    # and similarly for bcftools and htslib
./configure --prefix=<path>
make
make install
{% endhighlight %}

After compiling the software, if everything runs without errors, the executable _samtools_ will be inside the main folder. Inside of the _htslib_ folder we will found other tools. Add both to the path as explained before for quick access. Also from the Broad Institute, the Integrative Genomic Viewer (_IGV_) software is one of the most common tools to visualize our _SAM_/_BAM_ file. It is a desktop tool Access the link to download [IGV](https://software.broadinstitute.org/software/igv/download) in a zip file for Linux. 

{% highlight Bash %}
unzip IGV_Linux_2.8.10_WithJava.zip
{% endhighlight %}

The generated folder contains the pre-build software so there is no need to install. We can run _IGV_ executing `./igv.sh` inside of the folder using the terminal. The terminal will be capture by the software and a window would appear.

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

Once we have located coordinates of the input reads, we can compute the final step. In this case, our reads are **single-end** and smaller than 100 bp (average: 36bp). Given the manual, we should run the _BWA-backtrack_ software with the subcomand _samse_ to generate a _SAM_ file output.

{% highlight Bash %}
bwa samse SARS-CoV-2-reference.fasta SARS-CoV-2_exper-SRX9197062.sai SARS-CoV-2_exper-SRX9197062.fasta
{% endhighlight %}

The output is a SAM file with around 919 MB weight and 6.5 million lines. It is recomendable to pipe the output into gzip to compress the file and avoid consuming extra storage (`bwa samse <files> | gzip -3`).

<p>&nbsp;</p>

**Quality control**

Once we have our _Illumina_ reads mapped against the reference, we need to check the quality of the mapping. The percentage of mapped reads is a global indicator of the overall sequencing accuracy and of the presence of contaminating DNA. We can see all the different available programs by running:

{% highlight Bash %}
java -jar picard.jar -h
{% endhighlight %}

In our case, we will use ***CollectWgsMetrics*** inside of the _Base Calling_ tools. We can run the flag _-h_ to see the required and optional arguments.

Before the quality control, we need to sort the _SAM_ file.

{% highlight Bash %}
java -jar picard.jar SortSam -h
java -jar picard.jar SortSam -I SARS-CoV-2_exper-SRX9197062.sam -O SARS-CoV-2_exper-SRX9197062_sorted.sam -SORT_ORDER coordinate
{% endhighlight %}

The process took less than a minute. Now, we can run the quality control:

{% highlight Bash %}
java -jar picard.jar CollectWgsMetrics -h
java -jar picard.jar CollectWgsMetrics -I SARS-CoV-2_exper-SRX9197062_sorted.sam -R SARS-CoV-2-reference.fasta -O SARS-CoV-2_exper-SRX9197062_quality-control.txt --INCLUDE_BQ_HISTOGRAM --READ_LENGTH 36
{% endhighlight %}

It took a bit more than a minute to run. The ouput includes three sections: run options, metrics and histogram. ***The downloaded data from NCBI SRA has been curated previously so we have a mean coverage of 0 and a genome territory of 29903 bp***. The values in the section _##METRICS CLASS_ shows what I said.

<p>&nbsp;</p>

**Visualization**

To visualize the mapped reads against the reference genome, we need to do two things:

	1. Load our reference genome and annotation.
	2. Load our SAM file

In order to do the first step, we need to go to _Genomes_ > _Load Genome from File_. Then, select the file _SARS-CoV-2-reference.fasta_ and load it. The upper part of the window should show the full genome with a value of 29 kb in the middle of the line. In the upper right corner we can modify the zoom. If you move the blue line to the maximum you would be able to see the nucleotides per position. To load the gene annotation we can go to _File_ > _Load from File_ and select the file _SARS-CoV-2-reference.gff3_.

For the second step, we first need to create an index file to load the _SAM_ file. We will use the software _Samtools_. The _SAM_ input have to be a gzip compress file with _bgzip_. After compressed, we will create the index.

{% highlight Bash %}
bgzip SARS-CoV-2_exper-SRX9197062_sorted.sam
samtools index SARS-CoV-2_exper-SRX9197062_sorted.sam.gz SARS-CoV-2_exper-SRX9197062_sorted.sam.bai
bgzip -d SARS-CoV-2_exper-SRX9197062_sorted.sam.gz # Decompress to read in IGV
{% endhighlight %}

After running the three commands, we can load the files in _IGV_ to see the alignments. With _IGV_ started, go to _File_ > _Load from File_ and select _SARS-CoV-2_exper-SRX9197062_sorted.sam_. The wizard will tell you that the index file was not found and you can click _Go_ to create one. Click _Go_ and wait until finished and loaded, it may take several minutes. It will create a _.sai_ file that is the index for our _SAM_ file. Once loaded, 

<img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/IGV_full.png" alt="IGV full" style="width:100%">

<p>&nbsp;</p>

<figure>
    <a href="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/IGV_variant_pos3037"><img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/IGV_variant_pos3037.png"></a>
    <figcaption>Figure 6: SARS-CoV-2 mapped reads against reference. Variant (T/C) located at position 3037.</figcaption>
</figure>

CollectMultipleMetrics
CollectInsertSizeMetrics
CollectAlignmentSummaryMetrics



