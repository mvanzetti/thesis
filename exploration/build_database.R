source("https://bioconductor.org/biocLite.R")

library(GenomicRanges)
library(rtracklayer)
library(IRanges)

refs_file <- read.csv("/Users/manuel/development/thesis/download/ENCODE/enhancer-like-bed_refs.csv", sep="\t")

# test: read first row
first <- refs_file[1,]
bed_filepath <- toString(row[["bed_filepath"]])
input <- import.bed(file.path(bed_filepath))


# test: read a single bed.gz file
data_dir <- "/Users/manuel/development/thesis/download/ENCODE/ENCSR002WAP/files/ENCFF257TWP"
data_dir2 <-	"/Users/manuel/development/thesis/download/ENCODE/ENCSR783FDB/files/ENCFF539CZK"
# encode_extra_cols <- c(signalValue = "numeric", pValue = "numeric",
#                           qValue = "numeric", peak = "integer")

input_test <- import.bed(file.path(data_dir, 'ENCFF257TWP.bed.gz'))
input_test2 <- import.bed(file.path(data_dir2, 'ENCFF539CZK.bed.gz'))

input_merged <- c(input_test, input_test2)

length(input_merged)
length(input_test) + length(input_test2)

input_seqnames <- seqnames(input_test)
input_ranges <- ranges(input_test)
input_strand <- strand(input_test)
input_seqinfo <- seqinfo(input_test)
input_mcols <- mcols(input_test)
input_seqlengths <- seqlengths(input_test)

df <- data.frame(seqnames=seqnames(input_test),
                 starts=start(input_test)-1,
                 ends=end(input_test),
                 names=elementMetadata(input_test)$name,
                 scores=elementMetadata(input_test)$score,
                 strands=strand(input_test),
                 source='ENCODE3',
                 assembly=refs_file[refs_file$accession=="ENCSR002WAP",]$assembly)

write.table(df, file="foo.bed", quote=F, sep="\t", row.names=F, col.names=F)

# test range obj from csv -> using adrenal glands DNase only experiment
staging_dir <- "/Users/manuel/development/thesis/staging/ENCODE/merged/"
enhancers_filename <- "ENCFF984FSL.csv"
enhancers_info <- read.csv(paste0(staging_dir, enhancers_filename), sep="\t")

enhancers_ranges <- makeGRangesFromDataFrame(enhancers_info, keep.extra.columns = TRUE)

# test superposition between two different ranges -> using adrenal glands DNase + H3K27ac experiment
enhancers_filename_2 <- "ENCFF952DET.csv"
enhancers_info_2 <- read.csv(paste0(staging_dir, enhancers_filename_2), sep="\t")

enhancers_ranges_2 <- makeGRangesFromDataFrame(enhancers_info_2, keep.extra.columns = TRUE)

ovlp <- findOverlaps(enhancers_ranges, enhancers_ranges_2, select = "all")
countOvlp <- countOverlaps(enhancers_ranges, enhancers_ranges_2, type=c("any", "start", "end", "within", "equal"))
ov <- min(length(unique(queryHits(ovlp))), length(unique(subjectHits(ovlp))))

library(VennDiagram)
grid.newpage()
draw.pairwise.venn(
  area1=length(enhancers_ranges),
  area2=length(enhancers_ranges_2),
  cross.area=ov,
  category=c("DNase only", "DNase + H3K27ac"),
  fill=c("steelblue", "blue3"),
  cat.cex=0.7)

# 1. test build new dataframe with new column with unique id like "Encyclopedia.Experiment.number" (in pyhton?)
# 2. test build new dataframe with overlapping references with columns: id, overlapped_by_id,
#   where overlapped_by_id is the id of a range that overlaps the range with a certain id 
#   i.e. the id is the subject, the overlapped_by_id is the query range
