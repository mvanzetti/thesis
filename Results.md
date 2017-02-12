
#OverlapsDB Test results
#####Manuel Vanzetti
___
## Table of Contents  
[Null Models](#null)  
[Fisher's Exact Test](#fisher)  
[Jaccard Similarity Index](#jaccard)  
[Relative Distance](#reldist)  
[Relative Distance](#bib)  


 
<a name="null"/>
## Null Models
[]()


Two main types of null overlap model have been built: 

- Null overlap model by **shuffle** resampling
- Null overlap model by **random** interval generation

### Shuffle resampling
Given a list of intervals A and a second list B:

1. The the genomic locations in B are randomly permuted among the given genome assembly (eg: hg19), obtaining say B'
2. The overlaps between A and B' are computed

This method ensures that the distribution of sizes along the genome in the B' list is the same of B, i.e. the method ensures the *conservation* of the size distribution.

The resulting overlap database A &#8745; B' can be used as a null model for the overlap databae A &#8745; B.

### Random Interval 
Given a list of intervarls A and a second list B:

1. The size of B (number of intervals) and the mean size of intervals inside B is computed
2. A random list of intervals among the specified genome assembly with size = mean(size of B) and with the same length of B is generated, say B'
3. The overlaps between A and B' are computed

This method does not conserve the size distribution of the B candidates but ensures the randomness.

<a name="fisher"/>
## Fisher's Exact Test
[]()

Here we perform Fisher’s exact test on the number of overlaps/unique intervals between 2 enhancer candidates lists. 

We may wish to know if the amount of overlap between the 2 sets of intervals is more than we would expect given their coverage and the size of the genome.

Considering the T-cell biosample from ENCODE and the FANTOM permissive list, performing a Fisher's Exact Test requesting a minimal 10% of overlap we obtain

` encode_bed_sorted.fisher(fantom_bed_sorted, f=0.50, genome='hg19')`

```
# Number of query intervals: 42001
# Number of db intervals: 43010
# Number of overlaps: 1000
# Number of possible intervals (estimated): 1813917
# phyper(1000 - 1, 42001, 1813917 - 42001, 43010, lower.tail=F)
# Contingency Table Of Counts
#_________________________________________
#           |  in -b       | not in -b    |
#     in -a | 1000         | 41001        |
# not in -a | 42010        | 1729906      |
#_________________________________________
# p-values for fisher's exact test
left		right		two-tail	ratio
0.56134		0.45147		0.88391		1.004

```
We can see the constructed contingency table and the pvalues for left, right and two-tail tests. From here, we can say that given the hg19 genome, it is unlikely that we would see as many overlaps as we do if the intervals from a and b were not related.

The main idea is that the lower is the right tail, the more are the overlap found than expected (randomly). The ratio could be interpreted as a measure of the enrichment found (odds ratio).

Notice that the total number of possible intervals in the above example was estimated to be 1813917. This is based on a heuristic that uses the mean sizes of intervals in the a and b sets and the size of the genome. The reported p-value will depend greatly on this. 

<a name="jaccard"/>
## Jaccard Similarity Index
[]()

The Jaccard statistic is used in set theory to represent the ratio of the intersection of two sets to the union of the two sets. Similarly, Favorov et al [1] reported the use of the Jaccard statistic for genome intervals: specifically, it measures the ratio of the number of intersecting base pairs between two sets to the number of base pairs in the union of the two sets. 


The used Jaccard statistics relies on the bedtools jaccard tool, which implements this statistic, yet modifies the statistic such that the length of the intersection is subtracted from the length of the union. As a result, the final statistic ranges from 0.0 to 1.0, where 0.0 represents no overlap and 1.0 represent complete overlap.

<a name="reldist"/>
## Relative Distance
[]()

Traditional approaches to summarizing the similarity between two sets of genomic intervals are based upon the number or proportion of intersecting intervals. However, such measures are largely blind to spatial correlations between the two sets where, dpesite consistent spacing or proximity, intersections are rare (for example, enhancers and transcription start sites rarely overlap, yet they are much closer to one another than two sets of random intervals). Favorov et al [1] proposed a relative distance metric that describes distribution of relative distances between each interval in one set nd the two closest intervals in another set (see figure above). If there is no spatial correlation between the two sets, one would expect the relative distances to be uniformaly distributed among the relative distances ranging from 0 to 0.5. If, however, the intervals tend to be much closer than expected by chance, the distribution of observed relative distances would be shifted towards low relative distance values

<a name="bib"/>
## Bibliography
[]()

```
[1] Exploring Massive, Genome Scale Datasets with the GenometriCorr Package.
Favorov A, Mularoni L, Cope LM, Medvedeva Y, Mironov AA, et al. (2012)
PLoS Comput Biol 8(5): e1002529. doi:10.1371/journal.pcbi.1002529
```
