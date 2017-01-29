
# OverlapsDB Tool
#####mvanzetti


### Table of Contents  
[Candidate Lists Pipeline](#candidate)  
[Finding Overlaps](#overlaps)  
[Building OverlapsDB](#db)  
[Environments](#envs)
 
<a name="candidate"/>
## Candidate Lists Pipeline
[]()
TODO: Provide a picture of the full pipeline


### ENCODE REST API Downloader tool
Describe the tool built to download and reclassify the 437 samples containing putative enhancers from ENCODE database.
Because of a bug in the ENCODE download tool, an implementation of a data scraping tool using the ENCODE REST APIs was necessary.


### ENCODE Resource Processor
TODO: Describe the tool built to process downloaded ENCODE files and the possibile obtainable outputs

### DNase and H3K27ac ENCODE candidate enhancers
The full 437 samples list was then filtered following [this approach](http://zlab-annotations.umassmed.edu/enhancers/methods), so the list was restricted to putative enhancers of 47 human cell types showing both DNase and H3K27ac signals. The H3K27ac histone mark is associated with enhancer activity.

### dbSUPER Downloader tool
[TODO] Describe the tool built to download the full list of super-enhancers in two fashions:

- A csv file with columns ['ID', 'Chrom', 'Start', 'End', 'Size', 'Associated Gene', 'Method', 'Rank', 'Cell/Tissue', 'Genome'] useful for enhancers classification
- The full BED file, crucial for overlap analysis



### FANTOM5 permissive enhancers
A preliminary analysis was performed in order to understand the putative enhancers and formats available from FANTOM5 experiments.

Some references:

- At [this link](http://pressto.binf.ku.dk/about.php) a brief description of the Human Transcribed Enhancer Atlas is available
- At [this link](http://enhancer.binf.ku.dk/enhancers.php) is possibile to download enhancer lists by primary cell or organ/tissue expression rate 

It seemed more suitable to start from [this resource](https://anderssonlab.org/tag/fantom5/) where the enhancer prediction process is briefly explained, according to Andersson et al. article *An atlas of active enhancers across human cell types and tissues*.

Moreover, [here](http://enhancer.binf.ku.dk/presets/) the resources used in the Andersson et al. article are available, the resources are divided in these sections:

1. Extensive enhancers, divided in ubiquitous (cells, organs), tss-enhancer associations, extensive, permissive
2. Enhancers specifically expressed in cells
3. Enhancers specifically expressed in organs/tissues
4. Enhancer Expression (expressed in all facets,  positively differentially expressed in each facet against any other facet, ...)
5. Enhancer - FANTOM Robust Promoter associations

Citing Andersson et al. "In total, 38,554 enhancers were transcribed at a significant expression level in at least one primary cell or tissue sample. Below, we refer to this set as the ‘robust set’ of enhancers and indicate whenever it was used. For all analyses, we use the whole (‘permissive’) set of 43,011 enhancers if not otherwise mentioned", the starting point choosen list was the **permissive enhancers list**.

### FANTOM Resource processor
[TODO] Describe the tool built to process downloaded FANTOM5 files and the possibile obtainable outputs

### Epigenomics Roadmap Downloader Tool
The built downloader is able to download all the DNaseI putative enhancers BED files with coordinates for regions of each region type for each epigenome separately (see [here](http://egg2.wustl.edu/roadmap/web_portal/DNase_reg.html#delieation) and [here](http://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_enh/)).


### Epigenomics Roadmap Resource Processor
The processor merges together all the epigenomes candidate enhancers bed file sproducing one big bed file using the metadata file. The final processed bed file is ready for overlapping analysis. 

Output example:

```
./process.py ROADMAP
2017-01-28 19:06:55,820 : MainProcess : INFO : ./process.py : OverlapsDB Process 🚀 🚀 🚀 🚀 🚀
2017-01-28 19:06:55,822 : MainProcess : INFO : ./process.py : Processing data for ROADMAP
2017-01-28 19:06:55,822 : MainProcess : INFO : ./process.py : Optional parameters: Force: False
2017-01-28 19:06:55,822 : MainProcess : INFO : ./process.py : Processing Epigenomics Roadmap data
Copying metadata to staging...
127 epigenomes found
Building unified bed file in staging...
Skipping: no bed file for E114
Skipping: no bed file for E115
Skipping: no bed file for E116
Skipping: no bed file for E117
Skipping: no bed file for E118
Skipping: no bed file for E119
Skipping: no bed file for E120
Skipping: no bed file for E121
Skipping: no bed file for E122
Skipping: no bed file for E123
Skipping: no bed file for E124
Skipping: no bed file for E125
Skipping: no bed file for E126
Skipping: no bed file for E127
Skipping: no bed file for E128
Skipping: no bed file for E129
Exporting processed file in bed format (substituting names with candidate_ids) to: /Users/manuel/development/thesis/staging/EpigenomicsRoadmap/processed.bed
Completed
2017-01-28 19:17:09,326 : MainProcess : INFO : ./process.py : Process completed
```

<a name="overlaps"/>
## Finding Overlaps
[]()
The following are some preliminary analysis useful to seek the best strategy to build a final database

### Overlapping ENCODE hg19 DNase and H3K27ac candidates and FANTOM5 permissive enhancers

**Tools**

R, GRanges, subsetByOverlaps

Overlaps between genomic ranges were built using the default overlap type (any), i.e. one of start, end, within or equal.

Notice that the ENCODE3 DNase+H3K27ac candidate enhancers are divided in 47 different samples of the following types:

- primary cell
- tissue
- in vitro differentiated cells
- immortalized cell line
- stem cell
- induced pluripotent stem cell line
 
while FANTOM5 permissive enhancers are derived from original 432 primary cell, 135 tissue and 241 cell line samples and listed without redundancy

**Full overlaps**

Used lists:

- 43,011 enhancers from FANTOM5 permissive list
- 1,785,580 enhancers from ENCODE3 DNase+H3K27ac list 

Overlaps:

- 266,881 enhancers from ENCODE overlapping FANTOM5
- 31,792 enhancers from FANTOM5 overlapping ENCODE

Here a sample R script used to find the overlaps

```
# try to open a big filtered file with DNase+H3K27ac candidate enhancers 
filtered_dir <- "/Users/manuel/development/thesis/staging/ENCODE/filtered/"
filtered_file <- "filtered_20170106.csv"
big.df <- read.csv(paste0(filtered_dir, filtered_file), sep="\t")
big.ranges <- makeGRangesFromDataFrame(big.df, keep.extra.columns = TRUE)

# load FANTOM5 permissive enhancers and put in range obj
permissive_dir <- "/Users/manuel/development/thesis/staging/FANTOM/permissive/"
permissive_file <- "PERMISSIVE.csv"
permissive.df <- read.csv(paste0(permissive_dir, permissive_file), sep="\t")
permissive.ranges <-  makeGRangesFromDataFrame(permissive.df, keep.extra.columns = TRUE)

# find ENCODE candidates overlapping FANTOM5 permissive enhancers
encode.overlapping.fantom <- subsetByOverlaps(big.ranges, permissive.ranges)
fantom.overlapping.encode <- subsetByOverlaps(permissive.ranges, big.ranges)

cat(sprintf("regions of ENCODE overlapping regions from FANTOM are %d", length(encode.overlapping.fantom)))
cat(sprintf("regions of FANTOM overlapping regions from ENCODE are %d", length(fantom.overlapping.encode)))

```
Results:

- regions from ENCODE3 DNase+H3K27ac candidate enhancers overlapping regions from FANTOM5 permissive enhancers are 266881
- regions of FANTOM5 permissive enhancers overlapping regions from ENCODE3 DNase+H3K27ac candidate enhancers are 31792

**Overlaps by sample type (eg: primary cell)**

Used lists:

- 43,011 enhancers from FANTOM5 permissive list
- 659,128 primary cell enhancers from ENCODE3 DNase+H3K27ac list 

Overlaps:

- 119,852 enhancers from ENCODE overlapping FANTOM5

A possibile approach could be the one that provides an overlap database for each sample type, for example:

```
big.df.primarycell <-  subset(big.df, biosample_type == "primary cell")
unique(big.df.primarycell['biosample_term_name'])
unique(big.df.primarycell['organ_slims'])

encode.pc.ranges <- makeGRangesFromDataFrame(big.df.primarycell, keep.extra.columns = TRUE)

# ENCODE3 primary cells enhancers overlapping FANTOM5 permissive enhancers
encode.pc.overlapping.fantom <- subsetByOverlaps(encode.pc.ranges, permissive.ranges)
cat(sprintf("regions of ENCODE primary cells overlapping regions from FANTOM are %d", 
            length(encode.pc.overlapping.fantom)))

```
Result: 

- regions of ENCODE primary cells overlapping regions from FANTOM are 119852

**Overlaps by sample (eg: T-cell)**

Used lists:

- 43,011 enhancers from FANTOM5 permissive list
- 42,001 T-cell sample primary cell type enhancers from ENCODE3 DNase+H3K27ac list 

Overlaps:

- 6,482 enhancers from ENCODE overlapping FANTOM5

Another possibile approach could provide for each single ENCODE3 DNase+H3K27ac sample the overlaps with the FANTOM5 permissive enhancers, for example:

```
# filter by single ENCODE sample and find overlaps: primary cell type, T-cell sample
big.df.pc.Tcell <-  subset(big.df.primarycell, biosample_term_name == "T-cell")
encode.pc.Tcell.ranges <- makeGRangesFromDataFrame(big.df.pc.Tcell, keep.extra.columns = TRUE)
encode.pc.Tcell.overlapping.fantom <- subsetByOverlaps(encode.pc.Tcell.ranges, permissive.ranges)
cat(sprintf("regions of ENCODE T-cell primary cells overlapping regions from FANTOM are %d", 
            length(encode.pc.Tcell.overlapping.fantom)))
```
Result: 

- regions of ENCODE T-cell primary cells overlapping regions from FANTOM are 6482


### Overlaps DB pipeline: ENCODE in FANTOM

This pipeline ensures quickly to preserve informations provided by ENCODE DNase+H3K27ac samples without the need of sample specific enhancer list from FANTOM experiments.

- Select a single sample ENCODE3 DNase+H3K27ac enhancer list
- Find enhancers of this list overlapping FANTOM5 permissive list
- Build an R dataframe consisting of the resulting ranges, adding the additional columns:
	- *encode*: set to True
	- *fantom*: set to True
	- *overlap_type*: the type of the overlap (e.g.: any, within, equals, ...)
	- *max_gap*: the eventually max gap applied (e.g.: if equals, max_gap could be 10, 20, ...)
	- *overlap_min*: the minimum overlap (num of positions) that the algorithm must ensure

	In this way is possibile to understand clearly the criteria applied to build the database and choose to produce different databases changing the criteria or the sources.
- Extract the non-overlapping enhancers in an R dataframe and add the additional columns with encode=True, fantom=False
- Unify the the overlapping and non-overlapping in a new R dataframe

This pipeline could be executed for all the 47 samples simply unifing the results for each sample iteration.

The defined function is the following

```
buildEncodeOverlappingFantomBySampleDataframe <- 
  function(encode.df, fantom.df, encode_sample_type, encode_sample_term_name, overlap_type='any', overlap_maxgap=, overlap_min){
    
    # subset by sample type, term name
    query.subset.df <- encode.df
    query.subset.df <- subset(query.subset.df, biosample_type == encode_sample_type)
    query.subset.df <- subset(query.subset.df, biosample_term_name == encode_sample_term_name)
    
    # build ranges
    query.subset.ranges <- makeGRangesFromDataFrame(query.subset.df, keep.extra.columns = TRUE)
    subject.ranges <- makeGRangesFromDataFrame(fantom.df, keep.extra.columns = TRUE)
    
    query.overlapping.subject <- 
      subsetByOverlaps(query.subset.ranges, subject.ranges, type = overlap_type, maxgap = overlap_maxgap, minoverlap = overlap_min)
      
    overlap.df <- as.data.frame(query.overlapping.subject)
    overlap.df['encode'] = TRUE
    overlap.df['fantom'] = TRUE
    overlap.df['overlap_type'] = overlap_type
    overlap.df['overlap_maxgap'] = overlap_maxgap
    overlap.df['overlap_min'] = overlap_min
    
    overlap.df.subsetcols <- 
      subset(overlap.df, select = c('candidate_id', 'encode', 'fantom', 'overlap_type', 'overlap_maxgap', 'overlap_min'))
    
    merged.df <- merge(overlap.df.subsetcols, query.subset.df, by='candidate_id', all = TRUE)
    
    merged.df['encode'][is.na(merged.df['encode'])] <- TRUE
    merged.df['fantom'][is.na(merged.df['fantom'])] <- FALSE
    
    cat(
      sprintf(
        "Found %d candidates from ENCODE %s %s overlapping candidates from FANTOM with overlap type '%s', maxgap %d, min overlap %d positions",
        length(query.overlapping.subject), encode_sample_type, encode_sample_term_name, overlap_type, overlap_maxgap, overlap_min
        )
      )
    
    return(merged.df)
  }


```
The min overlap criteria function is

```
computeMinOverlapByMeanLen <-  function(encode.df, fantom.df, factor) {
  mean.encode = mean(encode.df$end - encode.df$start)
  mean.fantom = mean(permissive.df$end - permissive.df$start)
  return(round(min(mean.encode, mean.fantom) * factor))
}

```

Two examples of usage (on primary cell/T-cell and tissue/placenta) with 70% of min overlap requested:

```
> minoverlap.Tcell <- computeMinOverlapByMeanLen(big.df, permissive.df, 0.7)
> minoverlap.placenta <- computeMinOverlapByMeanLen(big.df, permissive.df, 0.7)
> test.primaryCell.Tcell.df <- 
+   buildEncodeOverlappingFantomBySampleDataframe(
+     big.df, permissive.df, 'primary cell', 'T-cell', 'any', 0L, minoverlap.Tcell)
Found 4380 candidates from ENCODE primary cell T-cell overlapping candidates from FANTOM with overlap type 'any', maxgap 0, min overlap 202 positions
> 
> test.tissue.placenta.df <- 
+   buildEncodeOverlappingFantomBySampleDataframe(
+     big.df, permissive.df, 'tissue', 'placenta', 'any', 0L, minoverlap.placenta)
Found 4199 candidates from ENCODE tissue placenta overlapping candidates from FANTOM with overlap type 'any', maxgap 0, min overlap 202 positions
```

Two examples of usage (on primary cell/T-cell and tissue/placenta) with no min overlap are (1 position is the default):

```
> test.primaryCell.Tcell.df <- 
+   buildEncodeOverlappingFantomBySampleDataframe(big.df, permissive.df, 'primary cell', 'T-cell', 'any', 0L)
Found 6482 candidates from ENCODE primary cell T-cell overlapping candidates from FANTOM with overlap type 'any' and maxgap 0
> test.tissue.placenta.df <- 
+   buildEncodeOverlappingFantomBySampleDataframe(big.df, permissive.df, 'tissue', 'placenta', 'any', 0L)
Found 5883 candidates from ENCODE tissue placenta overlapping candidates from FANTOM with overlap type 'any' and maxgap 0
```

where the dataframes are

```
filtered_dir <- "/Users/manuel/development/thesis/staging/ENCODE/filtered/"
filtered_file <- "filtered_20170106.csv"
big.df <- read.csv(paste0(filtered_dir, filtered_file), sep="\t")

permissive_dir <- "/Users/manuel/development/thesis/staging/FANTOM/permissive/"
permissive_file <- "PERMISSIVE.csv"
permissive.df <- read.csv(paste0(permissive_dir, permissive_file), sep="\t")
```
<a name="db"/>
## Building OverlapsDB
[]()

### Command Line Interface 
The following explains the utility 

```
usage: overlap.py [-h] [--assembly ASSEMBLY] [--method METHOD]
                  [--minoverlap MIN_OVERLAP] [--temp EXPORT_TEMP]
                  {ENCODE} {FANTOM,dbSUPER,All}

positional arguments:
  {ENCODE}              Select one of the available sources
  {FANTOM,dbSUPER,All}  Select the encyclopedia to consider

optional arguments:
  -h, --help            show this help message and exit
  --assembly ASSEMBLY   The assembly to use to build the overlaps files. For
                        the human genome assembly, type hg19.
  --method METHOD       Filter the source enhancer list by the provided method
  --minoverlap MIN_OVERLAP
                        The minimum overlap requested while overlapping.
  --temp EXPORT_TEMP    If specified, export temporary files representing
                        different overlapping phases

```

###ENCODE in dbSUPER

The built component is able to overlap the full files of ENCODE and dbSUPER

In the example showed below, the full list of 47 candidate enhancers from ENCODE DNase+H3K27ac is overlapped with the full human dbSUPER super-enhancer list and it's enriched with additional informations from both ENCODE and dbSUPER.
The min overlap requested is 10%.


Default options are `assembly='hg19', method='DNase_H3K27ac', min_overlap=0.1`

The resulting output is:

```
./overlap.py ENCODE dbSUPER
2017-01-15 20:44:50,273 : MainProcess : INFO : ./overlap.py : OverlapsDB Overlap 🚀 🚀 🚀 🚀 🚀
2017-01-15 20:44:50,275 : MainProcess : INFO : ./overlap.py : Initializing EncodeOverlapper
2017-01-15 20:44:50,275 : MainProcess : INFO : ./overlap.py : Overlapping with dbSUPER
ENCODE bed file: /Users/manuel/development/thesis/staging/ENCODE/filtered/filtered_hg19DNase_H3K27ac.bed
dbSUPER bed file: /Users/manuel/development/thesis/download/dbSUPER/all_hg19_bed.bed
Starting overlap...
8526833 ENCODE.intersect(dbSUPER) results
dbSUPER details file: /Users/manuel/development/thesis/download/dbSUPER/super-enhancers-annotations.csv
Merging details from dbSUPER...
ENCODE details file: /Users/manuel/development/thesis/staging/ENCODE/filtered/filtered_hg19DNase_H3K27ac.csv
Merging details from ENCODE...
Rearranging columns...
Exporting overlapped file to: /Users/manuel/development/thesis/overlap/filtered_hg19DNase_H3K27ac_dbSUPER_overlapped.csv
Completed
2017-01-15 21:00:53,869 : MainProcess : INFO : ./overlap.py : Overlapping completed
```

The resulting overlapping file size is about 2.66 GB 

The first three (transposed) records are

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>chrom</th>
      <td>chr3</td>
      <td>chr3</td>
      <td>chr1</td>
    </tr>
    <tr>
      <th>start</th>
      <td>152855118</td>
      <td>152855118</td>
      <td>214611302</td>
    </tr>
    <tr>
      <th>end</th>
      <td>152861069</td>
      <td>152861069</td>
      <td>214622352</td>
    </tr>
    <tr>
      <th>name</th>
      <td>ENCODE.3.ENCFF778PVS.0</td>
      <td>ENCODE.3.ENCFF778PVS.0</td>
      <td>ENCODE.3.ENCFF778PVS.1</td>
    </tr>
    <tr>
      <th>score</th>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>strand</th>
      <td>.</td>
      <td>.</td>
      <td>.</td>
    </tr>
    <tr>
      <th>size</th>
      <td>5951</td>
      <td>5951</td>
      <td>11050</td>
    </tr>
    <tr>
      <th>method</th>
      <td>DNase_H3K27ac</td>
      <td>DNase_H3K27ac</td>
      <td>DNase_H3K27ac</td>
    </tr>
    <tr>
      <th>description</th>
      <td>Enhancer-like regions using DNase and H3K27ac ...</td>
      <td>Enhancer-like regions using DNase and H3K27ac ...</td>
      <td>Enhancer-like regions using DNase and H3K27ac ...</td>
    </tr>
    <tr>
      <th>assembly</th>
      <td>hg19</td>
      <td>hg19</td>
      <td>hg19</td>
    </tr>
    <tr>
      <th>biosample_type</th>
      <td>primary cell</td>
      <td>primary cell</td>
      <td>primary cell</td>
    </tr>
    <tr>
      <th>biosample_term_id</th>
      <td>CL:0002327</td>
      <td>CL:0002327</td>
      <td>CL:0002327</td>
    </tr>
    <tr>
      <th>biosample_term_name</th>
      <td>mammary epithelial cell</td>
      <td>mammary epithelial cell</td>
      <td>mammary epithelial cell</td>
    </tr>
    <tr>
      <th>developmental_slims</th>
      <td>['ectoderm']</td>
      <td>['ectoderm']</td>
      <td>['ectoderm']</td>
    </tr>
    <tr>
      <th>system_slims</th>
      <td>['integumental system']</td>
      <td>['integumental system']</td>
      <td>['integumental system']</td>
    </tr>
    <tr>
      <th>organ_slims</th>
      <td>['mammary gland']</td>
      <td>['mammary gland']</td>
      <td>['mammary gland']</td>
    </tr>
    <tr>
      <th>encyclopedia</th>
      <td>ENCODE</td>
      <td>ENCODE</td>
      <td>ENCODE</td>
    </tr>
    <tr>
      <th>SE_chrom</th>
      <td>chr3</td>
      <td>chr3</td>
      <td>chr1</td>
    </tr>
    <tr>
      <th>SE_start</th>
      <td>152847915</td>
      <td>152847755</td>
      <td>214594593</td>
    </tr>
    <tr>
      <th>SE_end</th>
      <td>152884602</td>
      <td>152864205</td>
      <td>214629357</td>
    </tr>
    <tr>
      <th>SE_name</th>
      <td>SE_35899</td>
      <td>SE_64659</td>
      <td>SE_02460</td>
    </tr>
    <tr>
      <th>SE_score</th>
      <td>88</td>
      <td>441</td>
      <td>225</td>
    </tr>
    <tr>
      <th>SE_size</th>
      <td>36687.0</td>
      <td>16450.0</td>
      <td>34764.0</td>
    </tr>
    <tr>
      <th>SE_associated_gene</th>
      <td>RAP2B</td>
      <td>RAP2B</td>
      <td>PTPN14</td>
    </tr>
    <tr>
      <th>SE_method</th>
      <td>H3K27ac</td>
      <td>H3K27ac</td>
      <td>H3K27ac</td>
    </tr>
    <tr>
      <th>SE_biosample</th>
      <td>HMEC</td>
      <td>NHEK</td>
      <td>Astrocytes</td>
    </tr>
    <tr>
      <th>SE_ovlp_len</th>
      <td>5951</td>
      <td>5951</td>
      <td>11050</td>
    </tr>
    <tr>
      <th>SE_ovlp_pct</th>
      <td>100.0</td>
      <td>100.0</td>
      <td>100.0</td>
    </tr>
    <tr>
      <th>SE_encyclopedia</th>
      <td>dbSUPER</td>
      <td>dbSUPER</td>
      <td>dbSUPER</td>
    </tr>
  </tbody>
</table>
</div>




###ENCODE in FANTOM
The built component is able to overlap the full files of ENCODE and FANTOM permissive enhancers.

In the example showed below, the full list of 47 candidate enhancers from ENCODE DNase+H3K27ac is overlapped with the full list of FANTOM permissive enhancers and  enriched with additional informations from both ENCODE and FANTOM. The min overlap requested is 10%.

Default options are `assembly='hg19', method='DNase_H3K27ac', min_overlap=0.1`

The resulting output is:

```
./overlap.py ENCODE FANTOM
2017-01-15 20:18:41,039 : MainProcess : INFO : ./overlap.py : OverlapsDB Overlap 🚀 🚀 🚀 🚀 🚀
ENCODE bed file: /Users/manuel/development/thesis/staging/ENCODE/filtered/filtered_hg19DNase_H3K27ac.bed
FANTOM bed file: /Users/manuel/development/thesis/staging/FANTOM/permissive/PERMISSIVE.bed
Starting overlap...
1801781 ENCODE.intersect(FANTOM) results
Adding details from FANTOM...
ENCODE details file: /Users/manuel/development/thesis/staging/ENCODE/filtered/filtered_hg19DNase_H3K27ac.csv
Merging details from ENCODE...
Rearranging columns...
Exporting overlapped file to: /Users/manuel/development/thesis/overlap/filtered_hg19DNase_H3K27ac_FANTOM_overlapped.csv
Completed

```
The resulting overlapping file size is about 493.5 MB

The first three (transposed) overlapping records are


<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>chrom</th>
      <td>chr10</td>
      <td>chr3</td>
      <td>chr8</td>
    </tr>
    <tr>
      <th>start</th>
      <td>3892558</td>
      <td>5062817</td>
      <td>126230865</td>
    </tr>
    <tr>
      <th>end</th>
      <td>3895911</td>
      <td>5068862</td>
      <td>126234434</td>
    </tr>
    <tr>
      <th>name</th>
      <td>ENCODE.3.ENCFF778PVS.6</td>
      <td>ENCODE.3.ENCFF778PVS.8</td>
      <td>ENCODE.3.ENCFF778PVS.9</td>
    </tr>
    <tr>
      <th>score</th>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>strand</th>
      <td>.</td>
      <td>.</td>
      <td>.</td>
    </tr>
    <tr>
      <th>size</th>
      <td>3353</td>
      <td>6045</td>
      <td>3569</td>
    </tr>
    <tr>
      <th>method</th>
      <td>DNase_H3K27ac</td>
      <td>DNase_H3K27ac</td>
      <td>DNase_H3K27ac</td>
    </tr>
    <tr>
      <th>description</th>
      <td>Enhancer-like regions using DNase and H3K27ac ...</td>
      <td>Enhancer-like regions using DNase and H3K27ac ...</td>
      <td>Enhancer-like regions using DNase and H3K27ac ...</td>
    </tr>
    <tr>
      <th>assembly</th>
      <td>hg19</td>
      <td>hg19</td>
      <td>hg19</td>
    </tr>
    <tr>
      <th>biosample_type</th>
      <td>primary cell</td>
      <td>primary cell</td>
      <td>primary cell</td>
    </tr>
    <tr>
      <th>biosample_term_id</th>
      <td>CL:0002327</td>
      <td>CL:0002327</td>
      <td>CL:0002327</td>
    </tr>
    <tr>
      <th>biosample_term_name</th>
      <td>mammary epithelial cell</td>
      <td>mammary epithelial cell</td>
      <td>mammary epithelial cell</td>
    </tr>
    <tr>
      <th>developmental_slims</th>
      <td>['ectoderm']</td>
      <td>['ectoderm']</td>
      <td>['ectoderm']</td>
    </tr>
    <tr>
      <th>system_slims</th>
      <td>['integumental system']</td>
      <td>['integumental system']</td>
      <td>['integumental system']</td>
    </tr>
    <tr>
      <th>organ_slims</th>
      <td>['mammary gland']</td>
      <td>['mammary gland']</td>
      <td>['mammary gland']</td>
    </tr>
    <tr>
      <th>encyclopedia</th>
      <td>ENCODE</td>
      <td>ENCODE</td>
      <td>ENCODE</td>
    </tr>
    <tr>
      <th>FA_chrom</th>
      <td>chr10</td>
      <td>chr3</td>
      <td>chr8</td>
    </tr>
    <tr>
      <th>FA_start</th>
      <td>3893365</td>
      <td>5067974</td>
      <td>126231490</td>
    </tr>
    <tr>
      <th>FA_end</th>
      <td>3894190</td>
      <td>5068590</td>
      <td>126231859</td>
    </tr>
    <tr>
      <th>FA_name</th>
      <td>FANTOM.5.PERMISSIVE.3965</td>
      <td>FANTOM.5.PERMISSIVE.26422</td>
      <td>FANTOM.5.PERMISSIVE.39851</td>
    </tr>
    <tr>
      <th>FA_score</th>
      <td>693</td>
      <td>256</td>
      <td>140</td>
    </tr>
    <tr>
      <th>FA_size</th>
      <td>825</td>
      <td>616</td>
      <td>369</td>
    </tr>
    <tr>
      <th>FA_method</th>
      <td>CAGE_TCs</td>
      <td>CAGE_TCs</td>
      <td>CAGE_TCs</td>
    </tr>
    <tr>
      <th>FA_ovlp_len</th>
      <td>825</td>
      <td>616</td>
      <td>369</td>
    </tr>
    <tr>
      <th>FA_ovlp_pct</th>
      <td>24.6048</td>
      <td>10.1902</td>
      <td>10.339</td>
    </tr>
    <tr>
      <th>FA_encyclopedia</th>
      <td>FANTOM</td>
      <td>FANTOM</td>
      <td>FANTOM</td>
    </tr>
  </tbody>
</table>
</div>

###ENCODE in Epigenomics Roadmap

The built component is able to overlap the full files of ENCODE and Epigenomics Roadmap putative enhancers.

In the example showed below, the full list aggregating the 47 experiments of candidate enhancers from ENCODE DNase+H3K27ac is overlapped with the full list of the 111 Epigenomes DNaseI putative enhancers and enriched with additional informations from both ENCODE and Epigenomics Roadmap. The min overlap requested in this case is 40% (the sets are very large).

Options are `assembly='hg19', method='DNase_H3K27ac', min_overlap=0.4`

The resulting output is:

```
./overlap.py ENCODE ROADMAP --minoverlap 0.4 
2017-01-29 13:19:08,923 : MainProcess : INFO : ./overlap.py : OverlapsDB Overlap 🚀 🚀 🚀 🚀 🚀
2017-01-29 13:19:08,924 : MainProcess : INFO : ./overlap.py : Building OverlapsDB for ENCODE in ROADMAP
2017-01-29 13:19:08,925 : MainProcess : INFO : ./overlap.py : Optional parameters: Assembly: hg19, Method: DNase_H3K27ac, Min Overlap: 0.4, Export temp: False
2017-01-29 13:19:08,925 : MainProcess : INFO : ./overlap.py : Initializing EncodeOverlapper
2017-01-29 13:19:08,925 : MainProcess : INFO : ./overlap.py : Overlapping with Epigenomics Roadmap
ENCODE bed file: /Users/manuel/development/thesis/staging/ENCODE/filtered/filtered_hg19DNase_H3K27ac.bed
ROADMAP bed file: /Users/manuel/development/thesis/staging/EpigenomicsRoadmap/processed.bed
Starting overlap...
25604650 ENCODE.intersect(ROADMAP) results
ROADMAP details file: /Users/manuel/development/thesis/staging/EpigenomicsRoadmap/roadmap_metadata.csv
Merging details from ROADMAP...
ENCODE details file: /Users/manuel/development/thesis/staging/ENCODE/filtered/filtered_hg19DNase_H3K27ac.csv
Merging details from ENCODE...
Rearranging columns...
Exporting overlapped file to: /Users/manuel/development/thesis/overlap/filtered_hg19DNase_H3K27ac_ROADMAP_overlapped.csv
Completed
2017-01-29 15:57:32,090 : MainProcess : INFO : ./overlap.py : Overlapping completed
```

The resulting overlapping file size is about 9.4 GB

An example of three (transposed) overlapping records is the following

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>chrom</th>
      <td>chr5</td>
      <td>chr5</td>
      <td>chr5</td>
    </tr>
    <tr>
      <th>start</th>
      <td>148864272</td>
      <td>148864272</td>
      <td>148864272</td>
    </tr>
    <tr>
      <th>end</th>
      <td>148869801</td>
      <td>148869801</td>
      <td>148869801</td>
    </tr>
    <tr>
      <th>name</th>
      <td>ENCODE.3.ENCFF778PVS.5</td>
      <td>ENCODE.3.ENCFF778PVS.5</td>
      <td>ENCODE.3.ENCFF778PVS.5</td>
    </tr>
    <tr>
      <th>score</th>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>strand</th>
      <td>.</td>
      <td>.</td>
      <td>.</td>
    </tr>
    <tr>
      <th>size</th>
      <td>5529</td>
      <td>5529</td>
      <td>5529</td>
    </tr>
    <tr>
      <th>method</th>
      <td>DNase_H3K27ac</td>
      <td>DNase_H3K27ac</td>
      <td>DNase_H3K27ac</td>
    </tr>
    <tr>
      <th>description</th>
      <td>Enhancer-like regions using DNase and H3K27ac ...</td>
      <td>Enhancer-like regions using DNase and H3K27ac ...</td>
      <td>Enhancer-like regions using DNase and H3K27ac ...</td>
    </tr>
    <tr>
      <th>assembly</th>
      <td>hg19</td>
      <td>hg19</td>
      <td>hg19</td>
    </tr>
    <tr>
      <th>biosample_type</th>
      <td>primary cell</td>
      <td>primary cell</td>
      <td>primary cell</td>
    </tr>
    <tr>
      <th>biosample_term_id</th>
      <td>CL:0002327</td>
      <td>CL:0002327</td>
      <td>CL:0002327</td>
    </tr>
    <tr>
      <th>biosample_term_name</th>
      <td>mammary epithelial cell</td>
      <td>mammary epithelial cell</td>
      <td>mammary epithelial cell</td>
    </tr>
    <tr>
      <th>developmental_slims</th>
      <td>['ectoderm']</td>
      <td>['ectoderm']</td>
      <td>['ectoderm']</td>
    </tr>
    <tr>
      <th>system_slims</th>
      <td>['integumental system']</td>
      <td>['integumental system']</td>
      <td>['integumental system']</td>
    </tr>
    <tr>
      <th>organ_slims</th>
      <td>['mammary gland']</td>
      <td>['mammary gland']</td>
      <td>['mammary gland']</td>
    </tr>
    <tr>
      <th>encyclopedia</th>
      <td>ENCODE</td>
      <td>ENCODE</td>
      <td>ENCODE</td>
    </tr>
    <tr>
      <th>RO_chrom</th>
      <td>chr5</td>
      <td>chr5</td>
      <td>chr5</td>
    </tr>
    <tr>
      <th>RO_start</th>
      <td>148865036</td>
      <td>148865036</td>
      <td>148865036</td>
    </tr>
    <tr>
      <th>RO_end</th>
      <td>148867514</td>
      <td>148867514</td>
      <td>148867514</td>
    </tr>
    <tr>
      <th>RO_name</th>
      <td>ROADMAP.E001.84443</td>
      <td>ROADMAP.E002.28334</td>
      <td>ROADMAP.E005.61152</td>
    </tr>
    <tr>
      <th>RO_score</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>RO_strand</th>
      <td>.</td>
      <td>.</td>
      <td>.</td>
    </tr>
    <tr>
      <th>RO_size</th>
      <td>2478</td>
      <td>2478</td>
      <td>2478</td>
    </tr>
    <tr>
      <th>RO_method</th>
      <td>DNase</td>
      <td>DNase</td>
      <td>DNase</td>
    </tr>
    <tr>
      <th>RO_biosample_type</th>
      <td>PrimaryCulture</td>
      <td>PrimaryCulture</td>
      <td>ESCDerived</td>
    </tr>
    <tr>
      <th>RO_biosample_group</th>
      <td>ESC</td>
      <td>ESC</td>
      <td>ES-deriv</td>
    </tr>
    <tr>
      <th>RO_biosample_name</th>
      <td>ES-I3 Cells</td>
      <td>ES-WA7 Cells</td>
      <td>H1 BMP4 Derived Trophoblast Cultured Cells</td>
    </tr>
    <tr>
      <th>RO_biosample_anatomy</th>
      <td>ESC</td>
      <td>ESC</td>
      <td>ESC_DERIVED</td>
    </tr>
    <tr>
      <th>RO_ovlp_len</th>
      <td>2478</td>
      <td>2478</td>
      <td>2478</td>
    </tr>
    <tr>
      <th>RO_ovlp_pct</th>
      <td>44.81823114487249</td>
      <td>44.81823114487249</td>
      <td>44.81823114487249</td>
    </tr>
    <tr>
      <th>RO_encyclopedia</th>
      <td>ROADMAP</td>
      <td>ROADMAP</td>
      <td>ROADMAP</td>
    </tr>
  </tbody>
</table>
</div>



<a name="envs"/>
##Enviroments
[]()
### Conda environment
```
> conda list
# packages in environment at /Users/manuel/anaconda:
#
_license                  1.1                      py35_1  
_nb_ext_conf              0.3.0                    py35_0  
alabaster                 0.7.9                    py35_0  
anaconda                  custom                   py35_0  
anaconda-clean            1.0.0                    py35_0  
anaconda-client           1.5.1                    py35_0  
anaconda-navigator        1.3.1                    py35_0  
appnope                   0.1.0                    py35_0  
appscript                 1.0.1                    py35_0  
argcomplete               1.0.0                    py35_1  
astroid                   1.4.7                    py35_0  
astropy                   1.2.1               np111py35_0  
babel                     2.3.4                    py35_0  
backports                 1.0                      py35_0  
bcftools                  1.3.1                         1    bioconda
beautifulsoup4            4.5.3                    py35_0  
bedtools                  2.26.0                        0    bioconda
bitarray                  0.8.1                    py35_0  
blaze                     0.10.1                   py35_0  
bokeh                     0.12.2                   py35_0  
boto                      2.42.0                   py35_0  
bottleneck                1.1.0               np111py35_0  
bowtie2                   2.2.8                    py35_2    bioconda
cffi                      1.7.0                    py35_0  
chest                     0.2.3                    py35_0  
click                     6.6                      py35_0  
cloudpickle               0.2.1                    py35_0  
clyent                    1.2.2                    py35_0  
colorama                  0.3.7                    py35_0  
conda                     4.2.13                   py35_0  
conda-build               2.0.2                    py35_0  
conda-env                 2.6.0                         0  
configobj                 5.0.6                    py35_0  
contextlib2               0.5.3                    py35_0  
cryptography              1.5                      py35_0  
curl                      7.49.0                        1  
cycler                    0.10.0                   py35_0  
cython                    0.24.1                   py35_0  
cytoolz                   0.8.0                    py35_0  
dask                      0.11.0                   py35_0  
datashape                 0.5.2                    py35_0  
decorator                 4.0.10                   py35_0  
dill                      0.2.5                    py35_0  
docutils                  0.12                     py35_2  
dynd-python               0.7.2                    py35_0  
entrypoints               0.2.2                    py35_0  
et_xmlfile                1.0.1                    py35_0  
fastcache                 1.0.2                    py35_1  
filelock                  2.0.6                    py35_0  
flask                     0.11.1                   py35_0  
flask-cors                2.1.2                    py35_0  
freetype                  2.5.5                         1  
get_terminal_size         1.0.0                    py35_0  
gevent                    1.1.2                    py35_0  
greenlet                  0.4.10                   py35_0  
h5py                      2.6.0               np111py35_2  
hdf5                      1.8.17                        1  
heapdict                  1.0.0                    py35_1  
htslib                    1.3.2                         0    bioconda
icu                       54.1                          0  
idna                      2.1                      py35_0  
imagesize                 0.7.1                    py35_0  
ipykernel                 4.5.0                    py35_0  
ipython                   5.1.0                    py35_0  
ipython_genutils          0.1.0                    py35_0  
ipywidgets                5.2.2                    py35_0  
itsdangerous              0.24                     py35_0  
jbig                      2.1                           0  
jdcal                     1.2                      py35_1  
jedi                      0.9.0                    py35_1  
jinja2                    2.8                      py35_1  
jpeg                      8d                            2  
jsonschema                2.5.1                    py35_0  
jupyter                   1.0.0                    py35_3  
jupyter_client            4.4.0                    py35_0  
jupyter_console           5.0.0                    py35_0  
jupyter_core              4.2.0                    py35_0  
lazy-object-proxy         1.2.1                    py35_0  
libdynd                   0.7.2                         0  
libpng                    1.6.22                        0  
libtiff                   4.0.6                         2  
libxml2                   2.9.2                         0  
libxslt                   1.1.28                        2  
llvmlite                  0.13.0                   py35_0  
locket                    0.2.0                    py35_1  
lxml                      3.6.4                    py35_0  
markupsafe                0.23                     py35_2  
matplotlib                1.5.3               np111py35_0  
mistune                   0.7.3                    py35_1  
mkl                       11.3.3                        0  
mkl-service               1.1.2                    py35_2  
mpmath                    0.19                     py35_1  
multipledispatch          0.4.8                    py35_0  
nb_anacondacloud          1.2.0                    py35_0  
nb_conda                  2.0.0                    py35_0  
nb_conda_kernels          2.0.0                    py35_0  
nbconvert                 4.2.0                    py35_0  
nbformat                  4.1.0                    py35_0  
nbpresent                 3.0.2                    py35_0  
networkx                  1.11                     py35_0  
nltk                      3.2.1                    py35_0  
nose                      1.3.7                    py35_1  
notebook                  4.2.3                    py35_0  
numba                     0.28.1              np111py35_0  
numexpr                   2.6.1               np111py35_0  
numpy                     1.11.1                   py35_0  
odo                       0.5.0                    py35_1  
openpyxl                  2.3.2                    py35_0  
openssl                   1.0.2j                        0  
pandas                    0.19.1              np111py35_0  
partd                     0.3.6                    py35_0  
path.py                   8.2.1                    py35_0  
pathlib2                  2.1.0                    py35_0  
patsy                     0.4.1                    py35_0  
pep8                      1.7.0                    py35_0  
perl-threaded             5.22.0                       10    bioconda
pexpect                   4.0.1                    py35_0  
pickleshare               0.7.4                    py35_0  
pillow                    3.3.1                    py35_0  
pip                       8.1.2                    py35_0  
pkginfo                   1.3.2                    py35_0  
ply                       3.9                      py35_0  
prompt_toolkit            1.0.3                    py35_0  
psutil                    4.3.1                    py35_0  
ptyprocess                0.5.1                    py35_0  
py                        1.4.31                   py35_0  
pyasn1                    0.1.9                    py35_0  
pybedtools                0.7.8                    py35_1    bioconda
pycosat                   0.6.1                    py35_1  
pycparser                 2.14                     py35_1  
pycrypto                  2.6.1                    py35_4  
pycurl                    7.43.0                   py35_0  
pyflakes                  1.3.0                    py35_0  
pygments                  2.1.3                    py35_0  
pylint                    1.5.4                    py35_1  
pyopenssl                 16.0.0                   py35_0  
pyparsing                 2.1.4                    py35_0  
pyqt                      5.6.0                    py35_0  
pysam                     0.9.1.4                  py35_1    bioconda
pytables                  3.2.3.1             np111py35_0  
pytest                    2.9.2                    py35_0  
python                    3.5.2                         0  
python-dateutil           2.5.3                    py35_0  
python.app                1.2                      py35_4  
pytz                      2016.6.1                 py35_0  
pyyaml                    3.12                     py35_0  
pyzmq                     15.4.0                   py35_0  
qt                        5.6.0                         0  
qtawesome                 0.3.3                    py35_0  
qtconsole                 4.2.1                    py35_1  
qtpy                      1.1.2                    py35_0  
readline                  6.2                           2  
redis                     3.2.0                         0  
redis-py                  2.10.5                   py35_0  
requests                  2.11.1                   py35_0  
rope                      0.9.4                    py35_1  
ruamel_yaml               0.11.14                  py35_0  
samtools                  1.3.1                         5    bioconda
scikit-image              0.12.3              np111py35_1  
scikit-learn              0.17.1              np111py35_2  
scipy                     0.18.1              np111py35_0  
seaborn                   0.7.1                    py35_0  
setuptools                27.2.0                   py35_0  
simplegeneric             0.8.1                    py35_1  
singledispatch            3.4.0.3                  py35_0  
sip                       4.18                     py35_0  
six                       1.10.0                   py35_0  
snowballstemmer           1.2.1                    py35_0  
sockjs-tornado            1.0.3                    py35_0  
sphinx                    1.4.6                    py35_0  
spyder                    3.0.0                    py35_0  
sqlalchemy                1.0.13                   py35_0  
sqlite                    3.13.0                        0  
statsmodels               0.6.1               np111py35_1  
sympy                     1.0                      py35_0  
terminado                 0.6                      py35_0  
tk                        8.5.18                        0  
toolz                     0.8.0                    py35_0  
tornado                   4.4.1                    py35_0  
traitlets                 4.3.0                    py35_0  
unicodecsv                0.14.1                   py35_0  
urllib3                   1.12                     py35_0    bioconda
wcwidth                   0.1.7                    py35_0  
werkzeug                  0.11.11                  py35_0  
wheel                     0.29.0                   py35_0  
widgetsnbextension        1.2.6                    py35_0  
wrapt                     1.10.6                   py35_0  
xlrd                      1.0.0                    py35_0  
xlsxwriter                0.9.3                    py35_0  
xlwings                   0.10.0                   py35_0  
xlwt                      1.1.2                    py35_0  
xz                        5.2.2                         0  
yaml                      0.1.6                         0  

```

### R environment
```
> sessionInfo()
R version 3.3.2 (2016-10-31)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: macOS Sierra 10.12.1

locale:
[1] it_IT.UTF-8/it_IT.UTF-8/it_IT.UTF-8/C/it_IT.UTF-8/it_IT.UTF-8

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] plyr_1.8.4           Gviz_1.18.1          GenomicRanges_1.26.2 GenomeInfoDb_1.10.2  IRanges_2.8.1       
[6] S4Vectors_0.12.1     BiocGenerics_0.20.0 

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.8                   biovizBase_1.22.0             lattice_0.20-34               Rsamtools_1.26.1             
 [5] Biostrings_2.42.1             assertthat_0.1                digest_0.6.11                 mime_0.5                     
 [9] R6_2.2.0                      backports_1.0.4               acepack_1.4.1                 RSQLite_1.1-1                
[13] BiocInstaller_1.24.0          httr_1.2.1                    ggplot2_2.2.1                 zlibbioc_1.20.0              
[17] GenomicFeatures_1.26.2        lazyeval_0.2.0                data.table_1.10.0             rpart_4.1-10                 
[21] Matrix_1.2-7.1                checkmate_1.8.2               splines_3.3.2                 BiocParallel_1.8.1           
[25] AnnotationHub_2.6.4           stringr_1.1.0                 foreign_0.8-67                RCurl_1.95-4.8               
[29] biomaRt_2.30.0                munsell_0.4.3                 shiny_0.14.2                  httpuv_1.3.3                 
[33] rtracklayer_1.34.1            base64enc_0.1-3               htmltools_0.3.5               nnet_7.3-12                  
[37] SummarizedExperiment_1.4.0    tibble_1.2                    gridExtra_2.2.1               htmlTable_1.8                
[41] interactiveDisplayBase_1.12.0 Hmisc_4.0-2                   matrixStats_0.51.0            XML_3.98-1.5                 
[45] GenomicAlignments_1.10.0      bitops_1.0-6                  xtable_1.8-2                  gtable_0.2.0                 
[49] DBI_0.5-1                     magrittr_1.5                  scales_0.4.1                  stringi_1.1.2                
[53] XVector_0.14.0                latticeExtra_0.6-28           Formula_1.2-1                 RColorBrewer_1.1-2           
[57] ensembldb_1.6.2               tools_3.3.2                   dichromat_2.0-0               BSgenome_1.42.0              
[61] Biobase_2.34.0                yaml_2.1.14                   survival_2.40-1               AnnotationDbi_1.36.0         
[65] colorspace_1.3-2              cluster_2.0.5                 memoise_1.0.0                 VariantAnnotation_1.20.2     
[69] knitr_1.15.1   
```