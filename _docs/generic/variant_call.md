---
title: Variant Call
permalink: variant_call.html
sidebar: generic
# tags: [PopGen1]
#product: Generic
---

## Data

In this tutorial on how to do a variant call over genetic data we are going to continue with the data for _SARS-CoV-2_ used in the tutorial of [mapping reads against a reference genome](https://adriangeerre.github.io/popgen.github.io/mapping_reads.html). That means that the data (_BAM_ file) is obtained after running the mapping.

## Analysis

**Explanation**

Variant call analysis is a process to extract the variants or SNPs from our data in relation to a reference genome. These regions contain information that can be used to estimate different parameters in population genetics. The current input is a _SAM_/_BAM_ file and the current output a _VCF_ file. Samtools define a _VCF_ file as a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome. The format also has the ability to contain genotype information on samples for each position. We also need the reference genome file in fasta/fastq format. Find all the information about _VCF_ files [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

**Software**

The software we are going to use is ***Samtools***. In this tutorial we are going to use ***bcftools***.  Download the _bcftools_ from [here](http://www.htslib.org/download/). To install the software, move it to the folder of your interest and run: 

{% highlight Bash %}
bzip2 -d bcftools-<version>.tar.bz2
tar -xvf bcftools-<version>.tar
rm bcftool-<version>.tar
cd bcftools-<version>    # and similarly for bcftools and htslib
./configure --prefix=<path>
make
make install
{% endhighlight %}

Remove the source code folder. After compiling the software, if everything runs without errors, the executable _bcftools_ will be inside the main folder (_bcftools-version_). Finally, add the folder to the path for quick access.

**Create a VCF file**

The steps to generate a VCF file are:

	1. Get sequences to align and reference sequence for aligment.
	2. Generate the SAM file with the alignment information of sequences against reference.
	3. Transform SAM to BAM file.
	4. Calculate the read coverage.
	5. Call the variants using the BAM file/s and the reference.

In this case we have already define our sorted _BAM_ file and we have the reference sequence. We have the indexes of both files. Remember to generate the data following the mapping of Illumina reads explained [here](https://adriangeerre.github.io/popgen.github.io/mapping_reads.html). There is also an explanation on how to install the software on Linux. We will start in the **step 4**. Read about the meaning of the argument by running `bcftools mpileup`.

{% highlight Bash %}
bcftools mpileup -f SARS-CoV-2-reference.fasta SARS-CoV-2_exper-SRX9197062.bam -O b -o SARS-CoV-2_exper-SRX9197062.bcf
{% endhighlight %}

After running the command, you may see the message **maximum number of reads per input file set to -d 250**. This message means that over 250 sequence depth the rest of them would not be taken into account to avoid excess memory usage. You can redefine this argument if you want. Now, with the _.bcf_ file, we can call the variants or SNPs. We need to define the ploidy, in our case, _SARS-CoV-2_ has a ploidy of 1 (indeed its genomes is a RNA string). Read about the meaning of the argument by running `bcftools call`.

{% highlight Bash %}
bcftools call --ploidy 1 -m -v SARS-CoV-2_exper-SRX9197062.bcf -o SARS-CoV-2_exper-SRX9197062.vcf
{% endhighlight %}

We can have a look at the file using the command `less`. The file contains a few variants (13) as defined by _IGV_ during the visualization of the mapping. Indeed, it contains the SNP in the position 3037 as the example image showed. However, for didactical purposes, we will filter the _VCF_ file to generate the final variant call output.

{% highlight Bash %}
vcfutils.pl varFilter SARS-CoV-2_exper-SRX9197062.vcf > SARS-CoV-2_exper-SRX9197062_final.vcf
{% endhighlight %}

After running the comand, and applying the default values ,the number of variants (SNPs) have reduce to 10. We can visualize the result using _IGV_ (as in the read mapping tutorial). We should load the reference, _.gff3_, _BAM_ and _VCF_ files.

<p>&nbsp;</p>

<a href="http://adriangeerre.github.io/popgen.github.io/analysis/variant_call/images/IGV_VCF.png">
	<img src="http://adriangeerre.github.io/popgen.github.io/analysis/variant_call/images/IGV_VCF.png" alt="Project reads" style="width:100%">
</a>

<p>&nbsp;</p>

Finally, the variant call (_VCF_ file) shows 10 red lines on top of the alignment. If we compare with the ones defined in the coverage we can identify the 3 SNPs that were removed. Moreover, the variant marked at the end of the genome (position 29872) seems to be a false positive because it contains a small number of reads defining it.

<p>&nbsp;</p>
