
rm(list=ls(all=TRUE)); gc()
################################################
# 16S rRNA metagenomes
################################################

#--- choose your flavor: ---
cur <- "~/metagenome/HSN"
#cur <- ""
setwd(cur)

source("functions/FirstFunctions.R")
source("functions/LoadersFunctions.R")
source("functions/WorkingField.R")
source("functions/SecondaryFunctions.R")
source("lib/loaders.R")

#loading libraries
library(stringr)
library(ecodist)
library(stringr)
library(ecodist)
library(data.table)
library(ape)
library(ggplot2)
library(scales)
library(MASS)
library(stringr)
library(matrixStats)
library(gridExtra)
library(grid)
library(phyloseq) #source("https://bioconductor.org/biocLite.R"), biocLite("phyloseq")
library(reshape)
library(gplots)
require(MASS)

