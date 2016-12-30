# Chip-Seq analysis for histone modification H3K27ac Tutorial

## Analysis preparation
The [European Nucleotide Archive] (http://www.ebi.ac.uk/ena) provides many types of raw sequencingdata, sequence assembly information and functional annotation. 

###FASTQ data download
Download the [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) data corresponding toChIP-seq experiment mapping the H3K27ac histone modification in mouse Embryonic Stem cells (mES cells) along with the input control sample (see [here](http://epigenie.com/guide-getting-started-with-chip-seq/) for further informations about chip-seq sequencing).

From bash, connect to the ftp server `ftp://ftp.sra.ebi.ac.uk/` via `ftp` command and download these files

```
get vol1/fastq/SRR066/SRR066787/SRR066787.fastq.gzget vol1/fastq/SRR066/SRR066766/SRR066766.fastq.gzget vol1/fastq/SRR066/SRR066767/SRR066767.fastq.gz
```

In the folder you have downloaded it, run

```
gunzip SRR066787.fastq.gzgunzip SRR066766.fastq.gzgunzip SRR066767.fastq.gz
```

### Install bowtie2, samtools and bedtools
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/bowtie2/README.html)

First install [conda](http://conda.pydata.org/docs/intro.html) package manager using [miniconda](http://conda.pydata.org/docs/install/quick.html#os-x-miniconda-install) installer. Then add the following channels

```
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
```

Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.
Install [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

```
conda install bowtie2
conda update bowtie2
```

SAM (Sequence Alignment/Map) format is a generic format for storing large nucleotide sequence alignments. 
[SAM Tools](http://samtools.sourceforge.net/) provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.

Install samtools using conda

```
conda install samtools
conda update samtools
```

The [bedtools](http://bedtools.readthedocs.io/en/latest/) utilities support a wide range of operations for interrogating and manipulating genomic features. 

Install bedtools using conda

```
conda install bedtools
conda update bedtools
```

### Install MACS
[MACS](https://github.com/taoliu/MACS/) uses python2. In order to install it using conda, first create a python2 environment, activate it, then install macs2

```
conda create --name py27env python=2
source activate py27env
conda install macs2

```
Remember to activate the proper environment before using macs.

### Reference Genome Indexes
Use Illumina's [iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html) collection: each iGenomes archive contains pre-built Bowtie2 indexes.
Download the Ensembl Mouse Reference Genome from [here](ftp://ussd-ftp.illumina.com/Mus_musculus/Ensembl/NCBIM37/) and unpack it.
The Bowtie2 index files are located under the path 

```
Mus_musculus_Ensembl_NCBIM37/Mus_musculus/Ensembl/NCBIM37/Sequence/Bowtie2Index/
```

### Quality Assessment
Sequenced reads are saved in .fastq files. The very first step in the analyses of sequencing results consistsin quality assessment. The R package *ShortRead* provides a *qa* function to perform this analysis.

In the following R code a vector with fastq file names is generated and then for each of these files and the custom qas function is applied in order to assess the quality of the reads in each file; then a Quality Assessment report is generated

```
source("https://bioconductor.org/biocLite.R")fls = list.files(dataDirectory, ".fastq$", full=TRUE)names(fls) = sub(".fastq", "", basename(fls))biocLite("ShortRead")library(ShortRead)qas = lapply(seq_along(fls), function(i, fls) qa(readFastq(fls[i]), names(fls)[i]), fls)qa = do.call(rbind, qas)rpt = report(qa,dest = 'QA_report.html')
```

### Read alignment
The next step is to align the reads to mm9 mouse genome assembly. This is done using Bowtie2 tool. 

```
bowtie2 -p 8 -x {path_to_index_files}/{index_prefix} -U {fastq_filename} -S {sam_filename}
```
In our cases:

```
bowtie2 -p 8 -x Mus_musculus_Ensembl_NCBIM37/Mus_musculus/Ensembl/NCBIM37/Sequence/Bowtie2Index/genome -U SRR066787.fastq -S ES_input.sam
bowtie2 -p 8 -x Mus_musculus_Ensembl_NCBIM37/Mus_musculus/Ensembl/NCBIM37/Sequence/Bowtie2Index/genome -U SRR066766.fastq -S H3K27ac_rep1.sam
bowtie2 -p 8 -x Mus_musculus_Ensembl_NCBIM37/Mus_musculus/Ensembl/NCBIM37/Sequence/Bowtie2Index/genome -U SRR066767.fastq -S H3K27ac_rep2.sam

```

Output example:

```
13662226 reads; of these:
  13662226 (100.00%) were unpaired; of these:
    1852897 (13.56%) aligned 0 times
    8365623 (61.23%) aligned exactly 1 time
    3443706 (25.21%) aligned >1 times
86.44% overall alignment rate
```

The resulting [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) files are next transformed to [BAM](http://genome.sph.umich.edu/wiki/BAM) files and filtered for best aligned reads using samtools. In this case only reads with mapping quality equal to or more than 40 are considered ```
samtools view -bS -q 40 ES_input.sam > ES_input_bestAlignment.bamsamtools view -bS -q 40 H3K27ac_rep1.sam > H3K27ac_rep1_bestAlignment.bamsamtools view -bS -q 40 H3K27ac_rep2.sam > H3K27ac_rep2_bestAlignment.bam
```### Read FilteringRemove PCR duplicates

```
samtools rmdup -s ES_input_bestAlignment.bam ES_input_filtered.bamsamtools rmdup -s H3K27ac_rep1_bestAlignment.bam H3K27ac_rep1_filtered.bamsamtools rmdup -s H3K27ac_rep2_bestAlignment.bam H3K27ac_rep2_filtered.bam
```

Output example:

```
[bam_rmdupse_core] 23 / 9986171 = 0.0000 in library '	'
```
It means that 23 aligned reads of 99861671 were duplicates and got removed.

Transform BAM files to [BED](http://www.ensembl.org/info/website/upload/bed.html) files using bedtools

```
bedtools bamtobed -i ES_input_filtered.bam > ES_input_filtered.bedbedtools bamtobed -i H3K27ac_rep1_filtered.bam > H3K27ac_rep1_filtered.bedbedtools bamtobed -i H3K27ac_rep2_filtered.bam > H3K27ac_rep2_filtered.bed
```

In order to be consistent with other tools, in the final step of data preprocessing add a 'chr' prefix to the chromosome names using[awk](https://en.wikipedia.org/wiki/AWK)

```
awk '$0="chr"$0' ES_input_filtered.bed > ES_input_filtered_ucsc.bedawk '$0="chr"$0' H3K27ac_rep1_filtered.bed > H3K27ac_rep1_filtered_ucsc.bedawk '$0="chr"$0' H3K27ac_rep2_filtered.bed > H3K27ac_rep2_filtered_ucsc.bed
```

The resulting files contains the following columns:

- *chrom*: name of the chromosome or scaffold. Chromosome names can be given with or without the 'chr' prefix.
- *chromStart*: start position of the feature in standard chromosomal coordinates (i.e. first base is 0)
- *chromEnd*: end position of the feature in standard chromosomal coordinates
populated if higher-numbered ones are used
- *name*: label to be displayed under the feature, if turned on in "Configure this page".
- *score*: a score between 0 and 1000 (quality scores)
- *strand*: defined as + (forward) or - (reverse)

For example:

```
head ES_input_filtered_ucsc.bed >>
chr2	147239339	147239375	SRR066787.66	40	+
chr14	59395610	59395646	SRR066787.69	42	-
chr6	83400347	83400383	SRR066787.70	42	+
chr2	145449796	145449832	SRR066787.64	40	-
chr16	23534349	23534385	SRR066787.74	42	+
chr18	61314957	61314993	SRR066787.78	40	-
chr4	45962178	45962214	SRR066787.63	40	+
chr9	3594089	3594125	SRR066787.62	40	-
chrX	140285356	140285392	SRR066787.82	42	-
chr12	106006035	106006071	SRR066787.77	42	+
```

[Here](http://bedtools.readthedocs.io/en/latest/content/general-usage.html) there is a further explanation of the general usage of the BED format.

Finally, isolate data for chromosome 6

```
awk '{if($1=="chr6") print $0}' ES_input_filtered_ucsc.bed > ES_input_filtered_ucsc_chr6.bedawk '{if($1=="chr6") print $0}' H3K27ac_rep1_filtered_ucsc.bed > H3K27ac_rep1_filtered_ucsc_chr6.bedawk '{if($1=="chr6") print $0}' H3K27ac_rep2_filtered_ucsc.bed > H3K27ac_rep2_filtered_ucsc_chr6.bed
```

### Peak Finding
**Narrow Peak Finding**

Using [MACS](https://github.com/taoliu/MACS/), [find peaks](https://en.wikipedia.org/wiki/Peak_calling) in the genome using the control data file 

```
macs2 callpeak -t H3K27ac_rep1_filtered_ucsc.bed -c ES_input_filtered_ucsc.bed -f BED -g mm --nomodel -n Rep1macs2 callpeak -t H3K27ac_rep2_filtered_ucsc.bed -c ES_input_filtered_ucsc.bed -f BED -g mm --nomodel -n Rep2
```

The resulting output files are the following:

- *NAME_peaks.xls* is a tabular file which contains information about called peaks. You can open it in excel and sort/filter using excel functions. Information include:

	- chromosome name
	- start position of peak
	- end position of peak
	- length of peak region
	- absolute peak summit position
	- pileup height at peak summit 
	- -log10(pvalue) for the peak summit (e.g. pvalue = 1e-10, 	then this value should be 10)
	fold enrichment for this peak summit against random Poisson distribution with local 	lambda, 
	- -log10(qvalue) at peak summit	
	
	Coordinates in XLS is 1-based which is different with BED format.

- *NAME_peaks.narrowPeak* is BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue. You can load it to UCSC genome browser. The BED6+4 format is used to provide called peaks of signal enrichment based on pooled, normalized (interpreted) data. Columns are:

	- *chrom*: name of the chromosome (or contig, scaffold, etc.).
	chromStart - The starting position of the feature in the chromosome or 	scaffold. The first base in a chromosome is numbered 0.
	- *chromEnd*: the ending position of the feature in the chromosome or 	scaffold. The chromEnd base is not included in the display of the feature. 	For example, the first 100 bases of a chromosome are defined as 	chromStart=0, chromEnd=100, and span the bases numbered 0-99.
	- *name*: name given to a region (preferably unique). Use '.' if no name is 	assigned.
	- *score*: indicates how dark the peak will be displayed in the browser 	(0-1000). If all scores were '0' when the data were submitted to the DCC, 	the DCC assigned scores 1-1000 based on signal value. Ideally the average 	signalValue per base spread is between 100-1000.
	- *strand*: +/- to denote strand or orientation (whenever applicable). Use 	'.' if no orientation is assigned.
	- *signalValue*: measurement of overall (usually, average) enrichment for 	the region.
	- *pValue*: measurement of statistical significance (-log10). Use -1 if no 	pValue is assigned.
	- *qValue*: measurement of statistical significance using false discovery 	rate (-log10). Use -1 if no qValue is assigned.
	- *peak*: point-source called for this peak (relative summit position to 	peak start); 0-based offset from chromStart. 	Use -1 if no point-source called.

	The file can be loaded directly to UCSC genome browser. Remove the beginning track line 	if you want to analyze it by other tools.

- *NAME_summits.bed* is in BED format, which contains the peak summits locations for every peaks. The 5th column in this file is -log10pvalue the same as *NAME_peaks.bed*. If you want to find the motifs at the binding sites, this file is recommended. The file can be loaded directly to UCSC genome browser. Remove the beginning track line if you want to analyze it by other tools.


Finally, isolate data for chromosome 6

```awk '{if($1=="chr6") print $0}' Rep1_peaks.narrowPeak > Rep1_peaks_ucsc_chr6.narrowPeakawk '{if($1=="chr6") print $0}' Rep2_peaks.narrowPeak > Rep2_peaks_ucsc_chr6.narroPeak
```


**Broad Peak Finding**

Following the [ENCODE approach for histone modification ChIP-seq peak calling](https://sites.google.com/site/anshulkundaje/projects/encodehistonemods), use [MACS2](https://github.com/taoliu/MACS/) to identify broader regions of enrichment (broadPeaks) that pass a Poisson p-value threshold of 0.1 (using MACS2’s broad peak mode) and narrow peaks of contiguous enrichment (narrowPeaks) that pass a Poisson p-value threshold of 0.01

```
macs2 callpeak -t H3K27ac_rep1_filtered_ucsc.bed -c ES_input_filtered_ucsc.bed -f BED --broad -g mm -p 0.01 --broad-cutoff 0.1 --nomodel -n Rep1
macs2 callpeak -t H3K27ac_rep2_filtered_ucsc.bed -c ES_input_filtered_ucsc.bed -f BED --broad -g mm -p 0.01 --broad-cutoff 0.1 --nomodel -n Rep2

```

The following is an example of the arguments list output

```
# ARGUMENTS LIST:
# name = Rep2
# format = BED
# ChIP-seq file = ['H3K27ac_rep2_filtered_ucsc.bed']
# control file = ['ES_input_filtered_ucsc.bed']
# effective genome size = 1.87e+09
# band width = 300
# model fold = [5, 50]
# pvalue cutoff for narrow/strong regions = 1.00e-02
# pvalue cutoff for broad/weak regions = 1.00e-01
# qvalue will not be calculated and reported as -1 in the final output.
# Larger dataset will be scaled towards smaller dataset.
# Range for calculating regional lambda is: 1000 bps and 10000 bps
# Broad region calling is on
# Paired-End mode is off
```

Above the resulting output files we are interesting on:

- *NAME_peaks.broadPeak* is in BED6+3 format which is similar to narrowPeak file, except for missing the 10th column for annotating peak summits. Columns are:
	- *chrom*: name of the chromosome (or contig, scaffold, etc.).
	- *chromStart*: the starting position of the feature in the chromosome or 	scaffold. The first base in a chromosome is numbered 0.
	- *chromEnd*: the ending position of the feature in the chromosome or 	scaffold. The chromEnd base is not included in the display of the feature. 	For example, the first 100 bases of a chromosome are defined as 	chromStart=0, chromEnd=100, and span the bases numbered 0-99. If all scores 	were '0' when the data were submitted to the DCC, the DCC assigned scores 	1-1000 based on signal value. Ideally the average signalValue per base 	spread is between 100-1000.
	- *name*: name given to a region (preferably unique). Use '.' if no name is 	assigned.
	- *score*: indicates how dark the peak will be displayed in the browser 	(0-1000).
	- *strand* - +/- to denote strand or orientation (whenever applicable). Use 	'.' if no orientation is assigned.
	- *signalValue*: measurement of overall (usually, average) enrichment for 	the region.
	- *pValue*: measurement of statistical significance (-log10). Use -1 if no 	pValue is assigned.
	- *qValue*: measurement of statistical significance using false discovery 	rate (-log10). Use -1 if no qValue is assigned.


Isolate then data for chromosome 6

```awk '{if($1=="chr6") print $0}' Rep1_peaks.broadPeak > Rep1_peaks_ucsc_chr6.broadPeakawk '{if($1=="chr6") print $0}' Rep2_peaks.broadPeak > Rep2_peaks_ucsc_chr6.broadPeak
```


## Bibliography
- Creyghton, M. P., Cheng, A. W., Welstead, G. G., Kooistra, T., Carey, B. W., Steine, E. J., … Jaenisch, R. (2010). Histone H3K27ac separates active from poised enhancers and predicts developmental state. Proceedings of the National Academy of Sciences of the United States of America, 107(50), 21931–21936. http://doi.org/10.1073/pnas.1016071107



