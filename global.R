library(shiny)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(ggvis)
library(ggbio)
#devtools::install_github("timelyportfolio/d3vennR")
library(d3vennR)

source("olRanges.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Prep ggbio hg19 circos files
data("CRC", package = "biovizBase")
seqlevels(hg19sub) <- paste0("chr", seqlevels(hg19sub))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ExomeDepth data all
#Exome Depth Output Annotated with gene names and cytoband
ed1anno <- read.delim("data/sample1anno_exomeDepth.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ExomeCNV data all
#ExomeCNV Output Annotated with gene names and cytoband

ec1anno <- read.delim("data/sample1anno_exomeCNV.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE)

ec1anno <- read.delim("data/sample1anno_exomeCNV.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Depth of Coverage files

doc_samp1 <-read.delim("data/sample1.sample_interval_summary",header= TRUE, sep="\t", stringsAsFactors = FALSE)

doc1 <- doc_samp1 %>%
  separate(Target, into=c("chrom", "start", "end"), sep=":|-") %>%
  mutate(start=as.numeric(start), end=as.numeric(end), 
         rang = end-start+1, coverage = average_coverage) %>%
  select(chrom, start, end, rang, coverage)


#------------------------------------------------------------------------------------
#ExomeDepthCNV results

ed1 <-  ed1anno %>%
  select(chromosome, start, end, width=end-start+1, type, nexons, BF, reads.ratio, gene=exons.hg19, cytoband)


#ExomeCNV results
ec1 <- ec1anno %>%
  mutate(type=ifelse(ratio>1,"duplication", ifelse(ratio<=1,"deletion", "")))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Functions

