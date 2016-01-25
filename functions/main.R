rm(list=ls(all=TRUE)); gc()
################################################
# 16S rRNA metagenomes
################################################

#--- choose your flavor: ---
current_directory <- "~/metagenome/BaseRep16"
#current_directory <- ""
setwd(current_directory)

source("functions/FirstFunctions.R")
source("functions/LoadersFunctions.R")
source("functions/SecondaryFunctions.R")

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

########################################################
## working with alpha, beta -diversity                ##
########################################################
flog.info("start")

path.to.config <- "Pathway/Pathway.ini"
pathway <- ReadIni(path.to.config) 
totalTable <- Load(pathway)

tableTOP <- list(TopFeatures(totalTable, type = "Type.1", subset = "case", perc = 85),
                 TopFeatures(totalTable, type = "Type.1", subset = "control", perc = 85))

WriteTable(FamilyCase, OutdirCase, "FamilyCase") 
WriteTable(GenusCase, OutdirCase, "GenusCase") 
WriteTable(SpeciesCase, OutdirCase, "SpeciesCase") 
WriteTable(OtuCase, OutdirCase, "OtuCase")

message("start working with alpha diversity")
alphaDivCase <- (totalTable$AlphaDiv[
  which(rownames(totalTable$AlphaDiv) %in% totalTable$Meta[which(totalTable$Meta[,"Type.1"] %in% "case"),"samples_name"]),])
]

WriteTable(alphaDivCase, pathway$Case$OutdirCase, "AlphaDivCaseTbl"))
AlphaDivMean <- lapply(AlphaDivCaseL[match('AlphaDiversity', names(AlphaDivCaseL))], mean)
AlphaDivSd <- lapply(AlphaDivCaseL[match('AlphaDiversity', names(AlphaDivCaseL))], sd)
AlphaDivMeanSd <- c(AlphaDivMean, AlphaDivSd)
names(AlphaDivMeanSd) <- c("mean", "sd")

WriteTable (AlphaDivMeanSd, OutdirCase, "AlphaDivMeanAndSd") 

cairo_pdf((paste(OutdirCase, '/Graphs/alpha_boxplot.pdf', sep = "/")), width = 10,  height = 10)
boxplot(alphaDivCase$AlphaDiversity)
dev.off()

#control
alphaDivCase <- (totalTable$AlphaDiv[
  which(rownames(totalTable$AlphaDiv) %in% totalTable$Meta[which(totalTable$meta[,"Type.1"] %in% "control"),"samples_name"]),])
]

message("end start working with alpha diversity")

message("start sequencing statistics")
seqStatTable <- read.table(pathway$statistics$SeqStatTblFile, comment.char = "", quote ="", as.is = T, header = T)

seqStatList <- as.list(seqStatTable)
seqStatMean <- lapply(seqStatList[-match('samp_name', names(seqStatList))], mean)
seqStatSd <- lapply(seqStatList[-match('samp_name', names(seqStatList))], sd)
names(seqStatMean) <- c("mean_init", "mean_AF", "mean_AR", "mean_Mapp")
names(seqStatSd) <- c("sd_init", "sd_AF", "sd_AR", "sd_Mapp")
seqStatMeanSd <- c(seqStatMean, seqStatSd)
seqStatMeanSd.df <- t(as.data.frame(seqStatMeanSd))

WriteTable(seqStatMeanSd.df, pathway$Case$OutdirCase, "SeqStatMeanAndSd")
message("end sequencing statistics")

#--- make MDS ---

#Family
ddHSN<-bcdist(totalTable$family)
MDS(ddHSN, pathway$Case$GraphsDir, 
    totalTable$Meta[which(rownames(totalTable$Meta) %in% totalTable$Meta[which(totalTable$Meta[,"Type.1"] %in% "case"),"samples_name"]),],
  "Type.1", "Family")

#Genus
ddHSN<-bcdist(UniteMatrices(totalTable$Genus, GenusCtrl))
MDS(ddHSN, pathway$Case$GraphsDir, pathway$Case$MetaTable, "Type.1", "Genus")

#Species
ddHSN<-bcdist(UniteMatrices(SpeciesCase, SpeciesCtrl))
MDS(ddHSN, pathway$Case$GraphsDir , totalTable$MetaTable, "Type.1", "Species")

#--- end ---
