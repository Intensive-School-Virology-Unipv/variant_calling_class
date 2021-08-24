library(tidyverse)
library(Gviz)
library(VariantAnnotation)
library(GenomicFeatures)
setwd("/home/student/variantlab/calls")


##########################
## read variants into R
##########################

variants <- read_tsv("/home/student/variantlab/calls/sample1_variants.tsv")
bamfile <- "/home/student/variantlab/alignment/sample1_sorted.bam"


#######################################
## create vars in genomics coordinates
#######################################

vars <- GRanges(
  variants$REGION,
  IRanges(
    start = variants$POS,
    width = 1
  )
)

### add allele fraction for the alternative allele
vars$AF <- variants$ALT_FREQ

genomeGR <- readRDS("/home/student/DATA/variant_calling_data/refs/covid19_genome_model.RData")
genome <- "NC_045512.2"

chr <- as.character(unique(seqnames(genomeGR)))
gen <- genome(genomeGR)
gtrack <- GenomeAxisTrack()

options(ucscChromosomeNames=FALSE)
atrack <- AnnotationTrack(genomeGR, name = genome)

genetrack <- GeneRegionTrack(genomeGR,
                             genome = gen,
                             chromosome = chr,
                             name = genome,
                             showId = TRUE,
                             geneSymbol = TRUE,
                             symbol = genomeGR$gene)

variants <- AnnotationTrack(vars, chromosome = chr, genome = gen,
                            name = "Variants")

af <- DataTrack(vars, genome = gen,
                chromosome = chr,
                data = "AF",
                type = c("h", "p"),
                group = as.character(vars$Feature_Type),
                legend = TRUE,
                name = "Allele Fraction"
)

cov <- DataTrack(range = bamfile, genome = gen, type = "l",
                 name = "Coverage", window = -1,
                 chromosome = chr)


### given the size of the VM we need to limit the coordinates of the drawing

plotTracks(list(gtrack, genetrack, variants, af, cov), from = 20000, to = 29000)


