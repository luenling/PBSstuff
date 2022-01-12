

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
library(ChIPQC)
browseVignettes("ChIPQC")
library(DiffBind)
install.packages("tidyverse")
library(tidyverse)
install.packages("stringr")
library(stringr)
library(readxl)
library(tidyr)
library(ggplot2)
library(reshape2)


samples=read.csv("~/sample_sheet.csv", header=TRUE)
samples$Peaks=NA
samples$PeakCaller="raw"
setwd("data")
getwd()

experiment = ChIPQC(samples, annotation = "hg19")

ChIPQCreport(experiment)

gsub("trimmo", "trimmo_markdup", samples$bamReads)

samples$bamReads=gsub("trimmo", "trimmo_markdup", samples$bamReads)
samples$bamControl=gsub("trimmo", "trimmo_markdup", samples$bamControl)

setwd("~/Sambamba/")
experiment = ChIPQC(samples, annotation = "hg19")
ChIPQCreport(experiment)

