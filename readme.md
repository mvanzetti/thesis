
# Thesis

## Candidate lists pipeline
Provide a picture of the full pipeline

### ENCODE REST API Downloader tool
Describe the tool built to download and reclassify the 437 samples containing putative enhancers from ENCODE database.
Because of a bug in the ENCODE download tool, an implementation of a data scraping tool using the ENCODE REST APIs was necessary.


### ENCODE Resource Processor
Describe the tool built to process downloaded ENCODE files and the possibile obtainable outputs

### DNase and H3K27ac ENCODE candidate enhancers
The full 437 samples list was then filtered following [this approach](http://zlab-annotations.umassmed.edu/enhancers/methods), so the list was restricted to putative enhancers of 47 human cell types showing both DNase and H3K27ac signals. The H3K27ac histone mark is associated with enhancer activity.

### dbSUPER Downloader tool
Describe the tool built to download the full list of super-enhancers in two fashions:

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
Describe the tool built to process downloaded FANTOM5 files and the possibile obtainable outputs

## Finding Overlaps
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

	In this way is possibile to understand clearly the criteria applied to build the database and choose to produce different databases changing the criteria or the sources.
- Extract the non-overlapping enhancers in an R dataframe and add the additional columns with encode=True, fantom=False
- Unify the the overlapping and non-overlapping in a new R dataframe

This pipeline could be executed for all the 47 samples simply unifing the results for each sample iteration.

The defined function is the following

```
buildEncodeOverlappingFantomBySampleDataframe <- 
  function(encode.df, fantom.df, encode_sample_type, encode_sample_term_name, overlap_type, overlap_maxgap){
    
    # subset by sample type, term name
    query.subset.df <- encode.df
    query.subset.df <- subset(query.subset.df, biosample_type == encode_sample_type)
    query.subset.df <- subset(query.subset.df, biosample_term_name == encode_sample_term_name)
    
    # build ranges
    query.subset.ranges <- makeGRangesFromDataFrame(query.subset.df, keep.extra.columns = TRUE)
    subject.ranges <- makeGRangesFromDataFrame(fantom.df, keep.extra.columns = TRUE)
    
    query.overlapping.subject <- 
      subsetByOverlaps(query.subset.ranges, subject.ranges, type = overlap_type, maxgap = overlap_maxgap)
      
    overlap.df <- as.data.frame(query.overlapping.subject)
    overlap.df['encode'] = TRUE
    overlap.df['fantom'] = TRUE
    overlap.df['overlap_type'] = overlap_type
    overlap.df['overlap_maxgap'] = overlap_maxgap
    
    overlap.df.subsetcols <- 
      subset(overlap.df, select = c('candidate_id', 'encode', 'fantom', 'overlap_type', 'overlap_maxgap'))
    
    merged.df <- merge(overlap.df.subsetcols, query.subset.df, by='candidate_id', all = TRUE)
    
    merged.df['encode'][is.na(merged.df['encode'])] <- TRUE
    merged.df['fantom'][is.na(merged.df['fantom'])] <- FALSE
    
    cat(
      sprintf(
        "Found %d candidates from ENCODE %s %s overlapping candidates from FANTOM with overlap type '%s' and maxgap %d",
        length(query.overlapping.subject), encode_sample_type, encode_sample_term_name, overlap_type, overlap_maxgap
        )
      )
    
    return(merged.df)
  }

```

Two examples of usage (on primary cell/T-cell and tissue/placenta) are:

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

### Overlaps DB pipeline: ENCODE in dbSUPER


