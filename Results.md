
#OverlapsDB results
#####mvanzetti


### Table of Contents  
[Testing Results](#test)  
 
<a name="test"/>
##Testing OverlapsDB results
[]()

### Null Models

### Fisher's Exact Test
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




### Jaccard Similarity
The Jaccard statistic is used in set theory to represent the ratio of the intersection of two sets to the union of the two sets. Similarly, Favorov et al [1] reported the use of the Jaccard statistic for genome intervals: specifically, it measures the ratio of the number of intersecting base pairs between two sets to the number of base pairs in the union of the two sets. 

```
[1] Exploring Massive, Genome Scale Datasets with the GenometriCorr Package.
Favorov A, Mularoni L, Cope LM, Medvedeva Y, Mironov AA, et al. (2012)
PLoS Comput Biol 8(5): e1002529. doi:10.1371/journal.pcbi.1002529
```

The used Jaccard statistics relies on the bedtools jaccard tool, which implements this statistic, yet modifies the statistic such that the length of the intersection is subtracted from the length of the union. As a result, the final statistic ranges from 0.0 to 1.0, where 0.0 represents no overlap and 1.0 represent complete overlap.

### Relative Distance
