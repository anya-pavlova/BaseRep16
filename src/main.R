
rm(list=ls(all=TRUE)); gc()
################################################
# Project Parkinson_diseas - 16S rRNA metagenomes
################################################

# choose your flavor:
cur <- "~/metagenome/BaseRep16"
#cur <- ""
setwd(cur)

# adding functions files 
source("lib/config_project.R")
source("lib/functions.R")
source("lib/functions_special.R")
source("lib/functions_special2.R")
source("lib/loaders.R")
source("lib/visualization.R")
source("src/diff_abund.R")

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


# TODO: read register of all samples and prepare some meta-data variables and files
#source("src/.R")

#--- write TOP features in .txt file
#TRA - table of the relative abundance
WriteTable <- function (TRA, outdir, type)
{
  filename<-paste(outdir,"/",type, '.txt', sep="")
  write.table(TRA, filename, quote=F, sep='\t')
}

WriteTable (family_case, OutdirCase, "family_case") 
WriteTable (genus_case, OutdirCase, "genus_case") 
WriteTable (species_case, OutdirCase, "species_case") 
WriteTable (otu_case, OutdirCase, "otu_case") 

#--- end write TOP features in .txt file

# select samples, drop some and prepare the tags
#source("src/preprocess.R")

# compare with control
#source("src/diff_abund.R")

## visualization
source("src/viz.R")




