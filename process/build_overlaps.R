source("https://bioconductor.org/biocLite.R")

library(GenomicRanges)
library(rtracklayer)
library(IRanges)

# test range obj from csv -> using adrenal glands DNase only experiment
staging_dir <- "/Users/manuel/development/thesis/staging/ENCODE/merged/"
file1 <- "ENCFF984FSL.csv"
df1 <- read.csv(paste0(staging_dir, enhancers_filename), sep="\t")
ranges1 <- makeGRangesFromDataFrame(df1, keep.extra.columns = TRUE)

# test superposition between two different ranges -> using adrenal glands DNase + H3K27ac experiment
file2 <- "ENCFF952DET.csv"
df2 <- read.csv(paste0(staging_dir, file2), sep="\t")
ranges2 <- makeGRangesFromDataFrame(df2, keep.extra.columns = TRUE)

ovlp <- findOverlaps(subject = ranges1, query = ranges2, select = "all")
countOvlp <- countOverlaps(ranges1, ranges2, type=c("any", "start", "end", "within", "equal"))
ov <- min(length(unique(queryHits(ovlp))), length(unique(subjectHits(ovlp))))

library(VennDiagram)
grid.newpage()
draw.pairwise.venn(
  area1=length(ranges1),
  area2=length(ranges2),
  cross.area=ov,
  category=c("DNase only", "DNase + H3K27ac"),
  fill=c("steelblue", "blue3"),
  cat.cex=0.7)

ranges1.overlapping.ranges2 <- subsetByOverlaps(ranges1, ranges2)
ranges2.overlapping.ranges1 <- subsetByOverlaps(ranges2, ranges1)

cat(sprintf("1 overlapping 2 has %d results, while 2 overlapping 1 has %d results",
            length(ranges1.overlapping.ranges2), length(ranges2.overlapping.ranges1)))


common.regions = Reduce(subsetByOverlaps, list(ranges1, ranges2))

length(subsetByOverlaps(ranges1, ranges2))
length(subsetByOverlaps(ranges1, ranges2, type = 'any', maxgap = 0L))

#Assuming two candidates to be the same while
length(subsetByOverlaps(ranges1, ranges2, type = 'equal', maxgap = 10L))


cat(sprintf("regions of 1 overlapping regions from 2 are %d", length(common.regions)))

# 2. test build new dataframe with overlapping references with columns: id, overlapped_by_id,
#   where overlapped_by_id is the id of a range that overlaps the range with a certain id 
#   i.e. the id is the subject, the overlapped_by_id is the query range

from <- c(5, 2, 3, 3, 3, 2)
to <- c(11, 15, 5, 4, 5, 11)
id <- letters[1:6]

Hits(from, to, 7, 15, id)
Hits(from, to, 7, 15, id, sort.by.query=TRUE)

ovlp


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

# investigate sample types 
unique(big.df[, c('description','biosample_term_name')])
unique(big.df['biosample_type'])
unique(big.df['biosample_term_name'])

big.df.tissue <-  subset(big.df, biosample_type == "tissue")
unique(big.df.tissue['biosample_term_name'])
unique(big.df.tissue['organ_slims'])

big.df.primarycell <-  subset(big.df, biosample_type == "primary cell")
unique(big.df.primarycell['biosample_term_name'])
unique(big.df.primarycell['organ_slims'])

encode.pc.ranges <- makeGRangesFromDataFrame(big.df.primarycell, keep.extra.columns = TRUE)

# ENCODE3 primary cells enhancers overlapping FANTOM5 permissive enhancers
encode.pc.overlapping.fantom <- subsetByOverlaps(encode.pc.ranges, permissive.ranges)
cat(sprintf("regions of ENCODE primary cells overlapping regions from FANTOM are %d", 
            length(encode.pc.overlapping.fantom)))

# filter by single ENCODE sample and find overlaps: primary cell type, T-cell sample
overlap_type <-  'any'
overlap_maxgap <- 0L

big.df.pc.Tcell <-  subset(big.df.primarycell, biosample_term_name == "T-cell")
encode.pc.Tcell.ranges <- makeGRangesFromDataFrame(big.df.pc.Tcell, keep.extra.columns = TRUE)
encode.pc.Tcell.overlapping.fantom <- 
  subsetByOverlaps(encode.pc.Tcell.ranges, permissive.ranges, type = overlap_type, maxgap = overlap_maxgap)

cat(sprintf("regions of ENCODE T-cell primary cells overlapping regions from FANTOM are %d", 
            length(encode.pc.Tcell.overlapping.fantom)))

# fantom.overlapping.encode.pc.Tcell <- subsetByOverlaps(permissive.ranges, encode.pc.Tcell.ranges)
# cat(sprintf("regions of FANTOM overlapping regions from ENCODE T-cell primary cells are %d", 
#             length(fantom.overlapping.encode.pc.Tcell)))

# Build an R dataframe of single sample ENCODE enhancers overlapping FANTOM5 permissive enhancers
# select a single sample; considering now the primary cell type, T-cell sample and find overlaps
# consider the encode.pc.Tcell.overlapping.fantom and convert it to a data.frame
overlap.df <- as.data.frame(encode.pc.Tcell.overlapping.fantom)
overlap.df['encode'] = TRUE
overlap.df['fantom'] = TRUE
overlap.df['overlap_type'] = overlap_type
overlap.df['overlap_maxgap'] = overlap_maxgap

overlap.df.subsetcols <- 
  subset(overlap.df, select = c('candidate_id', 'encode', 'fantom', 'overlap_type', 'overlap_maxgap'))

merged.df <- merge(overlap.df.subsetcols, big.df.pc.Tcell, by='candidate_id', all = TRUE)

merged.df['encode'][is.na(merged.df['encode'])] <- TRUE
merged.df['fantom'][is.na(merged.df['fantom'])] <- FALSE


# group the operations in a function
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

# test the function
test.primaryCell.Tcell.df <- 
  buildEncodeOverlappingFantomBySampleDataframe(big.df, permissive.df, 'primary cell', 'T-cell', 'any', 0L)

test.tissue.placenta.df <- 
  buildEncodeOverlappingFantomBySampleDataframe(big.df, permissive.df, 'tissue', 'placenta', 'any', 0L)
