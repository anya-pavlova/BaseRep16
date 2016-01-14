
rm(list=ls(all=TRUE)); gc()
################################################
# Project Parkinson_and_Sklerosis - 16S rRNA metagenomes
################################################

#--- choose your flavor: ---
cur <- "~/metagenome/ParkinsonAndSklerosis"
#cur <- ""
setwd(cur)

source("functions/ConfigProject.R")
source("functions/ReadQiimeFunctions.R")
source("functions/SecondaryFunctions.R")
source("functions/LoadersFunctions.R")

#---adding functions files 
#source("lib/functions.R")
#source("lib/functions_special.R")
#source("lib/functions_special2.R")
#source("lib/loaders.R")
#source("src/diff_abund.R")
#source("functions/LoadersFunctions.R")
#source("functions/SelectionFunctions.R")
#source("functions/Visualization.R")
#source("functions/WorkingField.R")

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
