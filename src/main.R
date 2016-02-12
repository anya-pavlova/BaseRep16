rm(list=ls(all=TRUE)); gc()
################################################
# 16S rRNA metagenomes (Base Rep)
################################################

#--- choose your flavor: ---
current_directory <- "~/metagenome/BaseRep16"
#current_directory <- ""
setwd(current_directory)

source("lib/functions.R")
source("lib/visualization.R")

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
library(futile.logger)


flog.info("start")

path.to.config <- "src/pathway.ini"
pathway <- ReadIni(path.to.config) 
totalTable <- Load(pathway)

family.case <- totalTable$family[which(rownames(totalTable$family) 
                                           %in% totalTable$Meta[which(totalTable$meta[,"Type.1"] 
                                                                      %in% "case"),"samples_name"]),]

grups <- c("case", "control")
tableTOP <- lapply(grups, function(i) TopFeatures(totalTable, type = "Type.1", subset = i, perc = 85))
names(tableTOP) <- grups


#ToDO: запись таблиц топовых значений для неограниченного числа групп по family, genus, species, otu

flog.info("writing")
WriteTable(tableTOP[["case"]][["familyTOP"]], pathway$Case$OutdirCase, "familyTOPcase") 
WriteTable(tableTOP[["case"]][["genuseTOP"]], pathway$Case$OutdirCase, "genuseTOPcase") 
WriteTable(tableTOP[["case"]][["speciesTOP"]], pathway$Case$OutdirCase,"speciesTOPcase") 
WriteTable(tableTOP[["case"]][["otuTOP"]], pathway$Case$OutdirCase, "otuTOPcase") 
# tableTop

message("start working with alpha diversity")
alphaDivCase <- (totalTable$AlphaDiv[which(rownames(totalTable$AlphaDiv) 
                                           %in% totalTable$Meta[which(totalTable$Meta[,"Type.1"] 
                                                                      %in% "case"),"samples_name"]),])]

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
  which(rownames(totalTable$AlphaDiv) 
        %in% totalTable$Meta[which(totalTable$meta[,"Type.1"] 
                                   %in% "control"),"samples_name"]),])
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
    totalTable$Meta[which(rownames(totalTable$Meta) 
                          %in% totalTable$Meta[which(totalTable$Meta[,"Type.1"]
                                                     %in% "case"),"samples_name"]),],
  "Type.1", "Family")

#Genus
ddHSN<-bcdist(UniteMatrices(totalTable$Genus, GenusCtrl))
MDS(ddHSN, pathway$Case$GraphsDir, pathway$Case$MetaTable, "Type.1", "Genus")

#Species
ddHSN<-bcdist(UniteMatrices(SpeciesCase, SpeciesCtrl))
MDS(ddHSN, pathway$Case$GraphsDir , totalTable$MetaTable, "Type.1", "Species")

