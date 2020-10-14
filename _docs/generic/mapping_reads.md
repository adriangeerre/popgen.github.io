---
title: Mapping reads
permalink: mapping_reads.html
sidebar: generic
#tags: [PopGen1]
#product: Generic
---

## Data

**Explanation & Download data**

Everybody knows about Covid-19, Coronavirus, or SARS-CoV-2 (Severe acute respiratory syndrome coronavirus 2). Nowadays, given the situation all over the world, there is a lot of data available for this virus. We will download the reference genome from Wuhan and Illumina data to map. The Wuhan virus (Accession: NC_045512) contains around 29903 bp in an ss-RNA sequence. The reference genome (the first isolate virus) can be found in [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). To download the reference DNA sequence click on _FASTA_ under the title. Once inside, you can click on _Send to_, select _Complete Record_, _File_, and format _FASTA_. Finally, click on _create file_ and save as _SARS-CoV-2-reference.fasta_. Moreover, following the same procedure, we can download the gene annotation for the reference genome. Click on _Send to_, select _Complete Record_, _File_, and format _GFF3_. Save the file with the name _SARS-CoV-2-reference.gff3_ All the information in NCBI for the virus can be accessed from [here](https://www.ncbi.nlm.nih.gov/sars-cov-2/).

Illumina data (reads) can be found in [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) searching for _SARS-CoV-2_. In this case, I choose the accession [SRX9197062](https://www.ncbi.nlm.nih.gov/sra/SRX9197062[accn]) because of didactical purposes. The short genome of the virus allows a wider vision and an easier computation than, for example, the human genome. The following image shows a screenshot of the accession:

<p>&nbsp;</p>

<a href="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/project_selected.png">
	<img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/project_selected.png" alt="SARS-CoV-2 selected project" style="width:100%">
</a>

<p>&nbsp;</p>

The picture shows the information related to the accession. For example, the number of bases, the size of the file, the day it was publish and other data, such as, sequencing machine used. By clicking on the highlighted area of the picture you will access the sequencing run:

<p>&nbsp;</p>

<a href="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/run_of_project.png">
	<img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/run_of_project.png" alt="Project run" style="width:100%">
</a>

<p>&nbsp;</p>

Inside of the run, we have different sections: _Metadata_, _Analysis_, _Reads_, and _Data Access_. The first two are a recap of the data pre and post-analysis. We are interested in downloading the data so we access the section _Reads_ (highlighted in the previous picture).

<p>&nbsp;</p>

<a href="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/reads_of_run.png">
	<img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/reads_of_run.png" alt="Project reads" style="width:100%">
</a>

<p>&nbsp;</p>

Once inside of _Reads_, we can see a interactive table that contains all reads. To download the reads, click on _Filtered Download_. Finally, select both, the fasta and fastq files and download them (or download only the fastq and transform it into fasta). The download link is highlighted in the following picture:

<p>&nbsp;</p>

<a href="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/download_illumina.png">
	<img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/download_illumina.png" alt="Download Illumina reads" style="width:100%">
</a>

<p>&nbsp;</p>

A second option to download the data is to go to the tab _Data access_ and use the external links. For example, in this case, the data is stored inside _Amazon S3_. Following the URL we can directly download the data to our computer. In relation to the data, the size of the downloaded file, compress with gunzip (.gz) is of 204.6 MB. The decompressed file weighs 479 MB. The downloaded file contains a total of 6.479.048 Illumina reads where all of them have a length of 36 bp.

<p>&nbsp;</p>

**Available data**

The data, together with the pictures, is available in the GitHub of the website. Find the resources [here](https://github.com/adriangeerre/popgen.github.io/tree/master/analysis/mapping_reads).

The data was decompressed and divide into smaller pieces to reduce the size of the fasta file and allow Github to storage it.

{% highlight Bash %}
gzip -d SARS-CoV-2_exper-SRX9197062.fasta.gz
wc -l SARS-CoV-2_exper-SRX9197062.fasta # Count the number of lines
split -a 2 -d -l 1000000 SARS-CoV-2_exper-SRX9197062.fasta 'SARS-CoV-2_exper-SRX9197062.fasta_part'
{% endhighlight %}

The output was 13 files (part00 to part12, all under 50 MB) containing 1 million lines (500.000 reads) each. In order to re-create the original file, download all the files and inside of a Linux/Mac terminal run:

{% highlight Bash %}
cat SARS-CoV-2_exper-SRX9197062.fasta_part* > SARS-CoV-2_exper-SRX9197062.fasta
{% endhighlight %}

## Analysis

**Explanation**

In this tutorial, we are going to **map** Illumina reads against a reference genome (without a reference we will _de novo_ assembly of the genome). We are using reads coming from a Polymerase Chain Reaction (PCR) so we will have several times the same reads and also, reads that overlap between them. The reference genome will be the template that will help us to place each read in the proper location. Once we have our sample genome we can compare and define the differences (variants). The output of the mapping is a _SAM_ (Sequence Alignment/Map) file, read more about the format [here](http://samtools.github.io/hts-specs/SAMv1.pdf).

<p>&nbsp;</p>

**Software**

The software we are going to use in the first section is ***FastQC***. The software was developed by the Babraham Bioinformatics and has an open-source license. Read more and download from [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Once downloaded, unzip the file `unzip fastqc_<version>.zip`. The files inside are pre-build and the software is ready to use. The software can run as a desktop interface or a terminal program. Add the path for quick access.

The software we are going to use in the second section is ***BWA***. There are many available softwares for mapping reads, for example, TopHat, MAQ, or Bowtie2. This [article](https://academic.oup.com/bioinformatics/article/28/24/3169/245777) lists a large number of them. Different software may have different qualities or specializations depending on the input. I do not have any special interest in using _BWA_ for the tutorial, as long as it does a fast computation.

Download the software ***BWA*** from [here](http://bio-bwa.sourceforge.net/). Place the file in the folder you prefer and run the following to decompress:

{% highlight Bash %}
bzip2 -d bwa-<version>.tar.bz2
tar -xvf bwa-<version>.tar
{% endhighlight %}

That should create a folder called _bwa-version_ where the files to compile the program are contained. Now, we should compile the program by running the following (extracted from [Github](https://github.com/lh3/bwa)):

{% highlight Bash %}
cd bwa-<version>
bwa make
{% endhighlight %}

If everything works, an executable file called _bwa_ would be created. Then, we can add the program folder into the path to quick access bwa by running `export PATH=$PATH:/<path>/<to>/<bwa>` (temporal) or modifying the path by defining it inside the _.bashrc_ file.

The software we are going to use from the third section to the end are ***Samtools*** and ***IGV***. Download the _Samtools_ from [here](http://www.htslib.org/). To install the software, move it to the folder of your interest and run: 

{% highlight Bash %}
bzip2 -d samtools-<version>.tar.bz2
tar -xvf samtools-<version>.tar
rm samtool-<version>.tar
cd samtools-<version>    # and similarly for bcftools and htslib
./configure --prefix=<path>
make
make install
{% endhighlight %}

After compiling the software, if everything runs without errors, the executable _samtools_ will be inside the main folder (_samtools-version_). Inside of the _htslib_ folder, we will found other tools. Add both to the path as explained before for quick access. Also from the Broad Institute, the Integrative Genomic Viewer (_IGV_) software is one of the most common desktop tools to visualize our _SAM_/_BAM_ files. Access the link to download [IGV](https://software.broadinstitute.org/software/igv/download) in a zip file for Linux. 

{% highlight Bash %}
unzip IGV_Linux_2.8.10_WithJava.zip
{% endhighlight %}

The generated folder contains the pre-build software so there is no need to install. We can run _IGV_ executing `./igv.sh` inside of the folder using the terminal. The terminal will be capture by the software and a window would appear.

<p>&nbsp;</p>

**Sequencing Quality control**

Before we map our Illumina reads against the reference, we need to check the quality of the reads. Doing quality controls is something we should always do to avoid increasing the storage, the running time, and to avoid errors in our analysis. We can run the quality control for the compress fastq file directly.

{% highlight Bash %}
fastqc -f fastq SARS-CoV-2_exper-SRX9197062.fastq.gz
{% endhighlight %}

It took a bit more than a minute to run and it displays the progress. The outputs are two, a zip file, and an HTML file (find them [here](https://github.com/adriangeerre/popgen.github.io/tree/master/analysis/mapping_reads/QC)). The first one contains the raw metrics and all the output data while the HTML contains a group of plots summarising the quality tests. If you open the HTML file with a browser, you will see that there are 3 different icons in the left column: fail (red), warning (orange), and pass (green). In our case, the data has been curated before uploading so almost everything seems correct except for the _Per base sequence content_ and the _Sequence duplication levels_. From my point of view, given that I have not done the lab work, the data seems right with minor issues. We should try to correct those issues and, after, proceed to the mapping of the reads.

<p>&nbsp;</p>

**Correct reads errors**

If we found some issues with the reads we can try to solve them. For example, as the information of the [run](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR12718124) says in section _Analysis_ there are only 11.23% of the reads identified for SARS-CoV-2. The others are contaminants that could cause a distortion of our results. In this case, ***I am not running the correction because it takes time, storage and the data was already cleaned for us***. For this purpose, we can use different software: _BLAST_ (online or installed), _Kraken2_, or R scripting with the packages related to the databases _Silva_, _NCBI_ (_taxize_) or _Bold_ (_bold_). 

Also, other mistakes can be corrected after mapping/assembling the reads. For example, the software _Blobtools_ (previously known as _Blobology_) is an option to visualize the content of genome assembly datasets.

<p>&nbsp;</p>

**Mapping**

Once the program is installed, we need to do three steps:

	1. Index the reference fasta sequence.
	2. Find the coordinates of the input reads.
	3. Map reads against the reference sequence.

For the first step we need to run the following:

{% highlight Bash %}
bwa index SARS-CoV-2-reference.fasta
{% endhighlight %}

This should take a few seconds and should generate the following output:

```
SARS-CoV-2-reference.fasta.amb
SARS-CoV-2-reference.fasta.ann
SARS-CoV-2-reference.fasta.bwt
SARS-CoV-2-reference.fasta.pac
SARS-CoV-2-reference.fasta.sa
```

Once we have generated the database from the reference genome, we can compute the second step. Here we can define a few parameters in order to locate the reads, for example, the number of gaps, extensions, deletions, among others. **We will keep all values in default** but I recommend playing around and check the different outputs. This process will generate a binary file with format _.sai_ that will be used in the final step, in my case it took 56 seconds.

{% highlight Bash %}
bwa aln SARS-CoV-2-reference.fasta SARS-CoV-2_exper-SRX9197062.fasta > SARS-CoV-2_exper-SRX9197062.sai
{% endhighlight %}

Once we have located the coordinates of the input reads, we can compute the final step. In this case, our reads are **single-read** (SE) and smaller than 100 bp (average: 36bp). Single-read means that the reads have been sequenced from only one end, this method is fast, economical, and common in RNAseq studies. Given the manual, we should run the _BWA-backtrack_ software with the subcommand _samse_ to generate a _SAM_ file output.

{% highlight Bash %}
bwa samse SARS-CoV-2-reference.fasta SARS-CoV-2_exper-SRX9197062.sai SARS-CoV-2_exper-SRX9197062.fasta > SARS-CoV-2_exper-SRX9197062.sam
{% endhighlight %}

The output is a SAM file with around 919 MB weight and 6.5 million lines. It is recommendable to pipe the output into gzip to compress the file and avoid consuming extra storage (`bwa samse <files> | gzip -3`). Moreover, we can use the subcommand _bamse_ to generate a _BAM_ file (binary SAM file). However, we will transform our _SAM_ file to _BAM_ later in the tutorial.

<p>&nbsp;</p>

**Remove duplicates**

Duplicate reads can distort our findings. For example, can introduce PCR errors that will be treated as SNPs if they are not corrected. In order to correct those we need to do a few steps:

{% highlight Bash %}
samtools fixmate -m SARS-CoV-2_exper-SRX9197062.sam SARS-CoV-2_exper-SRX9197062_fixmate.sam
samtools sort SARS-CoV-2_exper-SRX9197062_fixmate.sam -o SARS-CoV-2_exper-SRX9197062_fixmate_sorted.sam
{% endhighlight %}

After adding the mate score tag, we are ready to remove the duplicates:

{% highlight Bash %}
samtools markdup SARS-CoV-2_exper-SRX9197062_fixmate_sorted.sam SARS-CoV-2_exper-SRX9197062_fixmate_sorted_markdup.sam
{% endhighlight %}

Since we are using single-read data the process would not reduce the number of sequences but it will improve the posterior visualization. The removal of duplicates is more likely to improve paired-end reads datasets. Now, we can remove the intermediate files and keep the last one computed (markdup).

<p>&nbsp;</p>

**Transform _SAM_ to _BAM_ format**

In order to run obtain the alignment metrics using _Samtools flagstat_ we need a _BAM_ file. _BAM_ format is just a binary version of the _SAM_ format.

{% highlight Bash %}
samtools view -bS -T SARS-CoV-2-reference.fasta SARS-CoV-2_exper-SRX9197062_fixmate_sorted_markdup.sam > SARS-CoV-2_exper-SRX9197062.bam
{% endhighlight %}

<p>&nbsp;</p>

**Alignment Quality control**

The percentage of mapped reads is a global indicator of the overall sequencing accuracy and of the presence of contaminating DNA. We will run the program in the terminal but feel free to launch the desktop window.

{% highlight Bash %}
samtools index SARS-CoV-2_exper-SRX9197062.bam # Result = .bai
samtools flagstat SARS-CoV-2_exper-SRX9197062.bam > SARS-CoV-2_exper-SRX9197062_Mapping_Metrics.txt
{% endhighlight %}

The output is a text file including a few lines. Since our data is single-read and was curated before uploading to NCBI, the results shown are poor. However, the procedure can be done for other data and generate different metrics.

```
6479048 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
6239545 + 0 duplicates
6275945 + 0 mapped (96.87% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

The output shows that close to 97% of the reads have been mapped. Obviously, the lines related to the pairing of reads are uninformative for us.

<p>&nbsp;</p>

**Visualization**

To visualize the mapped reads against the reference genome, we need to do two things:

	1. Load our reference genome and annotation.
	2. Load our BAM file

In order to do the first step, we need to go to _Genomes_ > _Load Genome from File_. Then, select the file _SARS-CoV-2-reference.fasta_ and load it. The upper part of the window should show the full genome with a value of 29 kb in the middle of the line. In the upper right corner, we can modify the zoom. If you move the blue line to the maximum you would be able to see the nucleotides per position. To load the gene annotation we can go to _File_ > _Load from File_ and select the file _SARS-CoV-2-reference.gff3_.

For the second step, we already created the required index for the _BAM_ file using _samtools index_ (_SARS-CoV-2_exper-SRX9197062.bam.bai_). If we want to load a _SAM_ file we need to create a _.sai_ format file. Also, if we load it in _IGV_ it will ask to create an index file with a pop-up wizard. To load the files in _IGV_ to see the alignments. With _IGV_ started, go to _File_ > _Load from File_ and select _SARS-CoV-2_exper-SRX9197062.bam_. It may take 1-2 minutes to load. Once loaded you may be able to see something similar to the pictures. I said similar because the first picture shows the mapping of reads of the full genome with all the depth and coverage. The second picture shows the region to the left (position 3037 bp), which is marked in red color and where we can find a variant (T instead of C). Also, I modify the appearance a bit, such as the color of the coverage track.

<a href="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/IGV_full.png">
	<img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/IGV_full.png" alt="IGV full" style="width:100%">
</a>

<p>&nbsp;</p>

<a href="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/IGV_variant_pos3037.png">
	<img src="http://adriangeerre.github.io/popgen.github.io/analysis/mapping_reads/images/IGV_variant_pos3037.png" alt="IGV variant pos3037" style="width:100%">
</a>

After the visualization, the errors found in the quality control seems not to be interfering with the alignment. Finally, you can obtain the consensus sequence of your mapping using _IGV_. Right-click on the box containing the reads (Alignment track), a submenu would be displayed. This menu allows you to define parameters or sort the reads, among others. To get the consensus sequence you can click on _Copy consensus sequence_ and paste into a file, adapt, and save to fasta format. If you want to visualize the consensus sequence click on _Quick consensus mode_. I recommend playing with the different options displayed by the menu.

<p>&nbsp;</p>
