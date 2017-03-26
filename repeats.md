
#OverlapsDB And RepeatMasker
#####Manuel Vanzetti
___
## Table of Contents  

[Genome Repeats](#rep)

[Bibliography](#bib)  

### Note on Reference Genome Assembly
ENCODE, FANTOM and dbSUPER enhancer candidates have been annotated with the hg19 genome build. 
For the sake of consistency, all the performed analysis have been made considering the same hg19 reference genome from [here](http://www.repeatmasker.org/species/hg.html)

<a name="null"/>
## RepeatMasker results on hg19
[]()


See [here](http://www.repeatmasker.org/webrepeatmaskerhelp.html#reading)

The annotation file contains the cross_match output lines. It lists all best matches (above a set minimum score) between the query sequence and any of the sequences in the repeat database or with low complexity DNA. The term "best matches" reflects that a match is not shown if its domain is over 80% contained within the domain of a higher scoring match, where the "domain" of a match is the region in the query sequence that is defined by the alignment start and stop. These domains have been masked in the returned masked sequence file. In the output, matches are ordered by query name, and for each query by position of the start of the alignment. 

Example:

```
 1306 15.6  6.2  0.0 HSU08988  6563  6781  (22462) C  MER7A    DNA/MER2_type    (0)   336   103
12204 10.0  2.4  1.8 HSU08988  6782  7714  (21529) C  TIGGER1  DNA/MER2_type    (0)  2418  1493
  279  3.0  0.0  0.0 HSU08988  7719  7751  (21492) +  (TTTTA)n Simple_repeat      1    33   (0)
 1765 13.4  6.5  1.8 HSU08988  7752  8022  (21221) C  AluSx    SINE/Alu        (23)   289     1
12204 10.0  2.4  1.8 HSU08988  8023  8694  (20549) C  TIGGER1  DNA/MER2_type  (925)  1493   827
 1984 11.1  0.3  0.7 HSU08988  8695  9000  (20243) C  AluSg    SINE/Alu         (5)   305     1
12204 10.0  2.4  1.8 HSU08988  9001  9695  (19548) C  TIGGER1  DNA/MER2_type (1591)   827     2
  711 21.2  1.4  0.0 HSU08988  9696  9816  (19427) C  MER7A    DNA/MER2_type  (224)   122     2
```

This is a sequence in which a Tigger1 DNA transposon has integrated into a MER7 DNA transposon copy. Subsequently two Alus integrated in the Tigger1 sequence. The simple repeat is derived from the poly A of the Alu element. The first line is interpreted like this:

```
  1306    = Smith-Waterman score of the match, usually complexity adjusted
        The SW scores are not always directly comparable. Sometimes
        the complexity adjustment has been turned off, and a variety of
        scoring-matrices are used.
  15.6    = % substitutions in matching region compared to the consensus
  6.2     = % of bases opposite a gap in the query sequence (deleted bp)
  0.0     = % of bases opposite a gap in the repeat consensus (inserted bp)
  HSU08988 = name of query sequence
  6563    = starting position of match in query sequence
  7714    = ending position of match in query sequence
  (22462) = no. of bases in query sequence past the ending position of match
  C       = match is with the Complement of the consensus sequence in the database
  MER7A   = name of the matching interspersed repeat
  DNA/MER2_type = the class of the repeat, in this case a DNA transposon 
            fossil of the MER2 group (see below for list and references)
  (0)     = no. of bases in (complement of) the repeat consensus sequence 
            prior to beginning of the match (so 0 means that the match extended 
            all the way to the end of the repeat consensus sequence)
  2418    = starting position of match in database sequence (using top-strand numbering)
  1465    = ending position of match in database sequence
```
An asterisk (*) in the final column (no example shown) indicates that there is a higher-scoring match whose domain partly (<80%) includes the domain of this match. 

Note that the SW score and divergence numbers for the three Tigger1 lines are identical. This is because the information is derived from a single alignment (the Alus were deleted from the query before the alignment with the Tigger element was performed). The program makes educated guesses about many fragments if they are derived from the same element (e.g. it knows that the MER7A fragments represent one insert). In a next version I can identify each element with a unique ID, if interest exists (this could help to represent repeats cleaner in graphic displays). 

## Processing
From the assembly chromosome naming scheme we know that in addition to the "regular" chromosomes, the hg19 browser contains nine haplotype chromosomes, 39 unplaced contigs, and 20 unlocalized contigs. For unlocalized contigs, the contig name is appended to the regular chromosome name, as in `chr1_gl000191_random`. If the chromosome is unknown, the contig is represented with the name "chrUn" followed by the contig identifier, as in `chrUn_gl000211`. Note that the chrUn contigs are no longer placed in a single, artificial chromosome as they have been in previous UCSC assemblies.

The processing removes all the annotations resulting from unplaced or unlocalized sequences on reference chromosomes

- The `chr*_random` sequences are unplaced sequence on those reference
chromosomes.

- The `chrUn_*` sequences are unlocalized sequences where the corresponding
reference chromosome has not been determined.


Then the full dataset is divided by repeat class family obtaining 66 different samples, for example:

- Satellite/telo
- LINE/L1
- DNA/hAT-Charlie
- SINE/MIR
- LINE/L2
- LINE/CR1
-  ...

## No enough overlaps, moving to closeness
The overlap analysis tool built is not able to find evidences of an enrichment of transposable elements inside enhancers. 
Roughly speaking, results are less than random overlap results.

The new need is to build a set of analyses considering the closeness of sequences. Qualitatively, both MIR and Alu elements are found to be enriched in some ENCODE samples (eg: neuronal stem cell) respect to random sequences. There is in fact an indication of an enrichment of repeated sequences within a +- 10 kb range around the inferred center (end-start)/2 of ENCODE enhancer candidates. 

It is important now to focus on the development of a null-model based approached able to measure z values based on closeness-related features and then perform the analyses on Encode-only, FANTOM-only and then "Encode and FANTOM" enhancers list.

## Fisher
[See here for hypergeom](http://stackoverflow.com/questions/6594840/what-are-equivalents-to-rs-phyper-function-in-python)

[See here for survival function](https://en.wikipedia.org/wiki/Survival_function)


<a name="bib"/>
## Bibliography
[]()

See [Smith-Waterman Algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) for aligning reads


