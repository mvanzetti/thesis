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


# try to open a big filtered file with DNase+H3K27ac candidate enhancers --> OK!!!
# test range obj from csv -> using adrenal glands DNase only experiment
filtered_dir <- "/Users/manuel/development/thesis/staging/ENCODE/filtered/"
filtered_file <- "filtered_20170106.csv"
big.df <- read.csv(paste0(filtered_dir, filtered_file), sep="\t")
big.ranges <- makeGRangesFromDataFrame(big.df, keep.extra.columns = TRUE)

