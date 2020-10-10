---
title: Mapping reads
permalink: mapping_reads.html
sidebar: generic
tags: [PopGen1]
#product: Generic
---

## Data

Everybody knows about Covid-19, Coronavirus or SARS-CoV-2 (Severe acute respiratory syndrome coronavirus 2). The information of the consensus sequence can be found [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). The Wuhan virus (Accession: NC_045512) contains around 29903 bp in a ss-RNA sequence. The reference genome (the first isolate virus) can be found in [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). To download the reference DNA sequence click on _FASTA_ under the title. Once inside, you can click on _Send to_, select _Complete Record_, _File_ and format _FASTA_. Finally, click on create file and save as _SARS-CoV-2-reference.fasta_.

To get _Illumina_ data can be found in [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) searching for _SARS-CoV-2_. In this case, I choose the accesion _SRX9197062_ because of didactical purposes. The short genome of the virus allows a wider vision and an easier computation than, for example, the human genome.




For the analysis, I downloaded the metadata from this [link](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRX9197062&o=acc_s%3Aa). Other information about the data, can be found [here](https://www.ncbi.nlm.nih.gov/sra/SRX9197062[accn]).

The size of the downloaded file, compress with gunzip (.gz) is of 206 MB. The Illumina reads downloaded have a length of 36. 

## Analysis

