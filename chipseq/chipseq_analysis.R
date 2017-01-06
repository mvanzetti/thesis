source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")
biocLite("BSgenome.Mmusculus.UCSC.mm9")
biocLite("biomaRt")
biocLite("Gviz")
biocLite("chipseq")
biocLite("VennDiagram")

setwd("/Users/manuel/Documents/Tesi/Chip-seq_tutorial/dev/")
dataDir <- "../data_process/"

#Quality Assessment
fls = list.files(dataDir, ".fastq$", full = TRUE)
names(fls) = sub(".fastq", "", basename(fls))

library(ShortRead)
qas = lapply(seq_along(fls), function(i, fls)
  qa(readFastq(fls[i]), names(fls)[i]), fls)
qa = do.call(rbind, qas)
rpt = report(qa, dest = 'QA_report.html')

#Chromosome length
library(BSgenome.Mmusculus.UCSC.mm9)
genome = BSgenome.Mmusculus.UCSC.mm9
si = seqinfo(genome)
si = si[paste0('chr', c(1:19, 'X', 'Y'))]

library(biomaRt)
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
               dataset = "mmusculus_gene_ensembl",
               host = "may2012.archive.ensembl.org")


#to fix Error in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]) :
#there is no package called ‘colorspace’: install.packages("colorspace", dependencies = TRUE)
fm = Gviz:::.getBMFeatureMap()
fm["symbol"] = "external_gene_id"

#we get a snapshot of the results for chromosome 6 starting at position 122530000 and ending at
#position 122900000. This region amongst others encodes a highly ES cell specific Nanog gene
#(see https://en.wikipedia.org/wiki/Homeobox_protein_NANOG). We first
#isolate gene models for this interval. The result bm is saved
library(Gviz)
bm = BiomartGeneRegionTrack(
  chromosome = 'chr6',
  genome = "mm9",
  start = 122530000,
  end = 122900000,
  biomart = mart,
  filter = list("with_ox_refseq_mrna" = TRUE),
  size = 4,
  name = "RefSeq",
  utr5 = "red3",
  utr3 = "red3",
  protein_coding = "black",
  col.line = NULL,
  cex = 7,
  collapseTranscripts = "longest",
  featureMap = fm
)

#Promoter isolation: 
#the following code is necessary to isolate gene models from the biomart databse
#The object egs contains the annotations of the most 
#external 5 and 3 prime UTRs for each gene model
listAttributes(mart)[1:3, ]
ds = useDataset('mmusculus_gene_ensembl', mart = mart)
chroms = 6
egs = getBM(
  attributes = c(
    'ensembl_gene_id',
    'external_gene_id',
    'chromosome_name',
    'start_position',
    'end_position',
    'strand'
  ),
  filters = 'chromosome_name',
  values = chroms,
  mart = ds
)

# Reading filtered chip-seq reads
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
input = import.bed(file.path(dataDir, 'ES_input_filtered_ucsc_chr6.bed'))
rep1 = import.bed(file.path(dataDir, 'H3K27ac_rep1_filtered_ucsc_chr6.bed'))
rep2 = import.bed(file.path(dataDir, 'H3K27ac_rep2_filtered_ucsc_chr6.bed'))

#Preparation of the ChIP-seq and control samples: read extension
library(chipseq)

prepareChiIPseq = function(reads) {
  frag.len = median(estimate.mean.fraglen(reads))
  cat(paste0('Median fragment size for this library is  ', round(frag.len)))
  reads.extended = resize(reads, width = frag.len)
  
  return(trim(reads.extended))
}

input = prepareChiIPseq(input)
rep1 = prepareChiIPseq(rep1)
rep2 = prepareChiIPseq(rep2)

#Binning chipseq and control data_process
#Tile the genome into non overlapping bins of size 200 bp
#using chromosome size informations of assembly mm9 

binsize = 200
bins = tileGenome(si['chr6'], tilewidth = binsize, cut.last.tile.in.chrom = TRUE)

#count reads in bins
BinChIPseq = function(reads, bins){
  mcols(bins)$score = countOverlaps(bins, reads)
  return(bins)
}

input.200bins = BinChIPseq(input, bins)
rep1.200bins = BinChIPseq(rep1, bins)
rep2.200bins = BinChIPseq(rep2, bins)

#plot coverage for 1000 bins starting from bin 200000
plot( 200000:201000, rep1.200bins$score[200000:201000],
      xlab="chr6", ylab="counts per bin" )

#the following are useful to export in bedgraph format
export(input.200bins,
       con='input_chr6.bedGraph',
       format = "bedGraph")
export(rep1.200bins,
       con='H3K27ac_rep1_chr6.bedGraph',
       format = "bedGraph")
export(rep2.200bins,
       con='H3K27ac_rep2_chr6.bedGraph',
       format = "bedGraph")

#Data visualization with Gviz
# We start with loading the gene models for chromosome 6 starting at position
# 122,530,000 and ending at position 122,900,000. 
# We focus on this region as it harbors the Nanog gene, which
# is stongly expressed in ES cells.

AT = GenomeAxisTrack()

plotTracks(c( bm, AT),
           from=122530000, to=122900000,
           transcriptAnnotation="symbol", window="auto",
           cex.title=1, fontsize=10 )

#add data_process track
input.track = DataTrack(input.200bins,
                        strand="*", genome="mm9", col.histogram='gray',
                        fill.histogram='black', name="Input", col.axis="black",
                        cex.axis=0.4, ylim=c(0,150))
rep1.track = DataTrack(rep1.200bins,
                       strand="*", genome="mm9", col.histogram='steelblue',
                       fill.histogram='black', name="Rep. 1", col.axis='steelblue',
                       cex.axis=0.4, ylim=c(0,150))
rep2.track = DataTrack(rep2.200bins,
                       strand="*", genome="mm9", col.histogram='steelblue',
                       fill.histogram='black', name="Rep. 2", col.axis='steelblue',
                       cex.axis=0.4, ylim=c(0,150))

# We observe a uniform coverage in the
# case of the input track and pronounced peaks of enrichment H3K27ac in promoter and intergenic regions.
# Importantly, H3K27ac enriched regions are easily identified.

plotTracks(c(input.track, rep1.track, rep2.track, bm, AT),
           from=122530000, to=122900000,
           transcriptAnnotation="symbol", window="auto",
           type="histogram", cex.title=0.7, fontsize=10 )

# We observe a uniform coverage in the
# case of the input track and pronounced peaks of enrichment H3K27ac in promoter and intergenic regions.
# Importantly, H3K27ac enriched regions are easily identified.

#ChipSeq Peaks

#Import peaks found with MACS
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                         qValue = "numeric")

peaks.rep1 <- import.bed(file.path(dataDir,'Rep1_peaks_ucsc_chr6.broadPeak'),
                        extraCols = extraCols_broadPeak)

peaks.rep2 <- import.bed(file.path(dataDir,'Rep2_peaks_ucsc_chr6.broadPeak'),
                         extraCols = extraCols_broadPeak)

#build tracks
peaks1.track = AnnotationTrack(peaks.rep1,
                               genome="mm9", name='Peaks Rep. 1',
                               chromosome='chr6',
                               shape='box',fill='blue3',size=2)
peaks2.track = AnnotationTrack(peaks.rep2,
                               genome="mm9", name='Peaks Rep. 2',
                               chromosome='chr6',
                               shape='box',fill='blue3',size=2)

#Visualize peaks on Nanog locus
plotTracks(c(input.track, rep1.track, peaks1.track,
             rep2.track, peaks2.track, bm, AT),
           from=122630000, to=122700000,
           transcriptAnnotation="symbol", window="auto",
           type="histogram", cex.title=0.7, fontsize=10 )

# We can see that MACS has succesfully identified regions showing high H3K27ac signal. We see that
# both biological replicates agree well, however, in some cases peaks are called only in one sample

# find the overlap between the peak sets of the two replicates
ovlp <- findOverlaps(peaks.rep1, peaks.rep2, select = "all")
countOvlp <- countOverlaps(peaks.rep1, peaks.rep2, type=c("any", "start", "end", "within", "equal"))
ov <- min(length(unique(queryHits(ovlp))), length(unique(subjectHits(ovlp))))

library(VennDiagram)
grid.newpage()
draw.pairwise.venn(
  area1=length(peaks.rep1),
  area2=length(peaks.rep2),
  cross.area=ov,
  category=c("rep1", "rep2"),
  fill=c("steelblue", "blue3"),
  cat.cex=0.7)
 
# We will focus only on peaks identified in both replicates (hereafter refered to as enriched areas). The
# enriched areas are colored in green.

enriched.regions = Reduce(subsetByOverlaps, list(peaks.rep1, peaks.rep2))

enr.reg.track = AnnotationTrack(enriched.regions,
                                genome="mm9", name='Enriched regions',
                                chromosome='chr6',
                                shape='box',fill='green3',size=2)

plotTracks(c(input.track, rep1.track, peaks1.track,
             rep2.track, peaks2.track, enr.reg.track,
             bm, AT),
           from=122630000, to=122700000,
           transcriptAnnotation="symbol", window="auto",
           type="histogram", cex.title=0.5, fontsize=10 )

# identify promoters overlapping H3K27ac peaks: identify TSS of mouse genome
egs$TSS = ifelse( egs$strand == "1", egs$start_position, egs$end_position )

# We consider regions of 200 bp around the TSS as promoters.
promoter_regions =
  GRanges(seqnames = Rle( paste0('chr', egs$chromosome_name) ),
          ranges = IRanges( start = egs$TSS - 200,
                            end = egs$TSS + 200 ),
          strand = Rle( rep("*", nrow(egs)) ),
          gene = egs$external_gene_id)

# discover how many of out the promoters overlap with a H3K27ac enriched regions.
ovlp2 = findOverlaps( enriched.regions, promoter_regions )
enriched_in_promoter.len <- length(unique(subjectHits(ovlp2)))
promoter.len <- length(promoter_regions)
cat(sprintf("%d of %d promoters are overlapped by an enriched region (%f %%)",
             enriched_in_promoter.len, promoter.len, 100*enriched_in_promoter.len/promoter.len))

ovlp2b = findOverlaps( promoter_regions, enriched.regions )
promoter_in_enriched.len <-  length( unique(subjectHits(ovlp2b)))
enriched.len <- length(enriched.regions)
cat(sprintf("%d of %d enriched regions overlap a promoter (%f %%).",
            promoter_in_enriched.len, enriched.len, 100*promoter_in_enriched.len/enriched.len))

# Is this a significant enrichment? To see, we first calculate how much chromosome 6 is part of a promotor
# region. The following command reduces the promotor list to non-overlapping intervals and sums up their
# widths
promotor_total_length = sum(width(reduce(promoter_regions)))

promotor_fraction_of_chromosome_6 = promotor_total_length / seqlengths(si)["chr6"]

cat(sprintf("%f %% of chromosome 6 length is covered by promotor regions.", 
            promotor_fraction_of_chromosome_6*100))

cat(
  sprintf("%.2f%% of promoters are overlapped by H3K27ac-enriched regions", 
          100*promoter_in_enriched.len/enriched.len),
  sprintf("even though promoters make up %.2f%% of the chromosome",
          promotor_fraction_of_chromosome_6
          ))
  
# binomial test tells if the enrichment is strong:
binom.test(promoter_in_enriched.len, enriched.len, promotor_fraction_of_chromosome_6)

#promoters ovelapped by a H3K27ac peak
pos.TSS <-  egs[ unique(queryHits(findOverlaps( promoter_regions, enriched.regions ))),]
pos.TSS[1:3,]

#Analysis of the distribution of H3K27ac around a subset of gene promoters

# we extend the analysed region to +/-1000 bp around the TSS. We divide each of these 2000 bp regions into
# 20 bins of 100 bp length each and order the bins with increasing position for genes on the '+' strand and
# decreasing for genes on the '-' strand.

#for each TSS, create a sequence of 20 bins each long 100 bp around the TSS 
# (the verse depends on the strand)
tiles = sapply( 1:nrow(pos.TSS), function(i)
  if( pos.TSS$strand[i] == "1" )
    pos.TSS$TSS[i] + seq( -1000, 900, length.out=20 )
  else
    pos.TSS$TSS[i] + seq( 900, -1000, length.out=20 ) )

# head(pos.TSS)
# pos.TSS$TSS[1] + seq(-1000, 900, length.out = 20)
# seq(-1000,900, length.out = 20)

# Now tile the promoter regions with consecutive 100bp tiles. For each region, we order the tiles
# according to the gene orientation. We create 20 tiles per promoter region.
tiles = GRanges(tilename = paste( rep( pos.TSS$ensembl_gene_id, each=20), 1:20, sep="_" ),
                seqnames = Rle( rep(paste0('chr', pos.TSS$chromosome_name), each=20) ),
                ranges = IRanges(start = as.vector(tiles),
                                 width = 100),
                strand = Rle(rep("*", length(as.vector(tiles)))),
                seqinfo=si)

#count how many reads are mapping to each tile (overlaps of tiles in rep1, rep2)
H3K27ac.p = countOverlaps(tiles, rep1) + countOverlaps(tiles, rep2)

# In the H3K27ac.p.matrix each row is a H3K27ac-enriched promoter. Each column
# corresponds to a consecutive 100bp tile of 2000 bp region around the TSS overlapping a H3K27ac peak.
# Since we have divided each promoter region in 20 tiles, we obtain a matrix with 20 columns and 784 rows
# (the number of promoters overlapping H3K27ac peak).
H3K27ac.p.matrix = matrix( H3K27ac.p, nrow=nrow(pos.TSS),
                           ncol=20, byrow=TRUE )
 
# plot the result 
# - as a heatmap 
# - as a plot of average values per each tile for all the included promoters


dev.new(width=10, height=20)

colors = colorRampPalette(c('white','red','gray','black'))(100)
layout(mat=matrix(c(1,2,0,3), 2, 2),
       widths=c(2,2,2),
       heights=c(0.5,5,0.5,5), TRUE)
par(mar=c(2,2,1.5,1))
image(seq(0, max(H3K27ac.p.matrix), length.out=100), 1,
      matrix(seq(0, max(H3K27ac.p.matrix), length.out=100),100,1),
      col = colors,
      xlab='Distance from TSS', ylab='',
      main='Number of reads', yaxt='n',
      lwd=3, axes=TRUE)

box(col='black', lwd=2)
image(x=seq(-1000, 1000, length.out=20),
      y=1:nrow(H3K27ac.p.matrix),
      z=t(H3K27ac.p.matrix[order(rowSums(H3K27ac.p.matrix)),]),
      col=colors,
      xlab='Distance from TSS (bp)',
      ylab='Promoters', lwd=2)
box(col='black', lwd=2)
abline(v=0, lwd=1, col='gray')

plot(x=seq(-1000, 1000, length.out=20),
     y=colMeans(H3K27ac.p.matrix),
     ty='b', pch=19,
     col='red4',lwd=2,
     ylab='Mean tag count',
     xlab='Distance from TSS (bp)')
abline(h=seq(1,100,by=5),
       v=seq(-1000, 1000, length.out=20),
       lwd=0.25, col='gray')
box(col='black', lwd=2)

#We observe a strong enrichment of H3K27ac modification right after the TSS and a weaker peak of
#H3K27ac at the region immediately upstream of the TSS.
