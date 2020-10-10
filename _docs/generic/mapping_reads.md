---
title: Mapping reads
permalink: mapping_reads.html
sidebar: generic
tags: [PopGen1]
#product: Generic
---

## Data

**Explanation**

Everybody knows about Covid-19, Coronavirus or SARS-CoV-2 (Severe acute respiratory syndrome coronavirus 2). The information of the consensus sequence can be found [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). The Wuhan virus (Accession: NC_045512) contains around 29903 bp in a ss-RNA sequence. The reference genome (the first isolate virus) can be found in [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). To download the reference DNA sequence click on _FASTA_ under the title. Once inside, you can click on _Send to_, select _Complete Record_, _File_ and format _FASTA_. Finally, click on create file and save as _SARS-CoV-2-reference.fasta_. All the information in NCBI for the virus can be access from [here](https://www.ncbi.nlm.nih.gov/sars-cov-2/).

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

A second option to download the data is to go to the tab _Data access_ and use the external links. For example, in this case the data is stored inside _Amazon S3_. Following the url we can direct download the data to our computer. In relation to the data, the size of the downloaded file, compress with gunzip (.gz) is of 206 MB. The Illumina reads downloaded have a length of 36. 

## Analysis

**Explanation**

