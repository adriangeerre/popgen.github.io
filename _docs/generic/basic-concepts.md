---
title: Basic concepts
permalink: basic-concepts.html
sidebar: generic
#tags: [PopGen1]
#product: Generic
---

As defined by _Rasmus Nielsen & Montgomery Slatkin_ in the book _An Introduction to Population Genetics_, Population Genetics is the study of alleles in populations. The definition included two concepts that are fundamental to understand, _population_ and _alleles_. A population, as a simple explanation, is a group of individual living in a defined area and capable of giving rise to new individuals by means of reproduction. In other words, individuals of a given species that share a common area (the scale goes from a small area to a full country or continent). The term allele refers to the genetic variants within a population. To understand the second term I need to be more specific. A locus (loci in plural), in Population Genetics, is a region of the genome where there are one or more alleles segregating (in other fields it referes to a _coding gene_). A locus can be a gene, a region of a gene or a single nucleotide. The nucleotide variation inside of a locus is defined as a segregating site or polymorphism or mutation or single nucleotide polymorphism (SNP). For example, if a given position of the genome the population shows that there are individuals carrying a T where others are carrying a C (position _e_ in Diagram 1), that site is an SNP. Something importat is to say that _indels_ (insertions or deletions) are not considered SNPs. In brief, this field of Biology uses the variation within the genomes of the different individuals of a population to estimate how it was, it is and it will be. This type of studies depend mainly on lab work, PCR or NGS, and computation, bioinformatic software but also in mathematical and statistical models.

<p>&nbsp;</p>	# Blank line

```
        a b c d e f g h i j k l m
ind 1:  A C C A T G A T A C T A C 
ind 2:  A C C A C G A T A C T A C 
ind 3:  A C C A C G A T G C T A C
ind 4:  A C C A C G A T G C T A C
ind 5:  A C C A C G A T A C T A C
```
**Diagram 1:** Comparison of a genetic region between 5 individuals of the same population of an _haploid_ species. The positions _e_ and _i_ are segregating sites or SNPs.

<p>&nbsp;</p>

**How do SNPs occur?**
SNPs affecting evolution, affecting the genetic material of gametes, occur during the division and recombination of the chromosomes. The replication of the DNA is not perfect and sometimes the process fails including errors. Most of the errors are corrected but a few of them remain. The errors that raise in a new individual are called mutations and the effect of them could be deletereous, neutral or advantageous (More information in the section _Selection_).

To start, a basic analysis applied to the populations is to measure the allele frequencies. To make it simple, first we are going to work with a di-allelic model. In other words, SNPs are going to be defined by a dual option even though there are 3 more potential nucleotides for a position. Moreover, mutations are not likely to occur and the most common situation is to determine a di-allelic behaviour of SNPs.

<p>&nbsp;</p>

<center><img src="https://latex.codecogs.com/svg.latex?F_A=\frac{N_A}{2N}&space;and&space;F_a=\frac{N_a}{2N}" title="Allele frequency (di-allelic model)"/></center>

<p>&nbsp;</p>

The formula contain various terms:

	- F: Is the frequency, the subindex relate to the allele.
	- A: First allele. It could define the major or dominant allele.
	- a: Second allele. It could define the minor or recesive allele.
	- N: Number of individuals in the population. It uses 2N because is defined for a diploid population, such as, humans. The subindex related to the number of individuals of a given allele.
	

Moreover, the sum of both frequencies is equal to 1. For example:

The position _e_ in the example above (Diagram 1) contains 1 T's and 4 C's. The frequencies in the haploid population are 0.20 (=1/5) and 0.80 (=4/5). In the position _i_ we have 0.60 and 0.40 for A and G, respectively.

