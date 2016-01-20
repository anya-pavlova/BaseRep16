########################################################
## working with alpha, beta -diversity                ##
########################################################

#--- TOP features ---
TotalTOPtable <- lapply(TotalTable, chooseTOPfeature, perc = 85)

FamilyTOP <- chooseTOPfeature(TotalTable$Family[which(rownames(TotalTable$Family)%in%),2] , perc=85)
GenusTOP <- chooseTOPfeature(TotalTable$Genus , perc=85)
SpeciesTOP <- chooseTOPfeature(TotalTable$Species , perc=85)
OtuTOP <- chooseTOPfeature(TotalTable$Otu , perc=85)

TotalTable$Meta


WriteTable (FamilyCase, OutdirCase, "FamilyCase") 
WriteTable (GenusCase, OutdirCase, "GenusCase") 
WriteTable (SpeciesCase, OutdirCase, "SpeciesCase") 
WriteTable (OtuCase, OutdirCase, "OtuCase") 
#--- end write TOP features in .txt file ---

#--- working with alpha diversity ---
AlphaDivCase <- LoadAlphaDiv(alpha_div_case)
WriteTable (AlphaDivCase, OutdirCase, "AlphaDivTbld") 

AlphaDivCaseL <- as.list(AlphaDivCase)
AlphaDivMean <- lapply(AlphaDivCaseL[match('AlphaDiversity', names(AlphaDivCaseL))], mean)
AlphaDivSd <- lapply(AlphaDivCaseL[match('AlphaDiversity', names(AlphaDivCaseL))], sd)
AlphaDivMeanSd <- c(AlphaDivMean, AlphaDivSd)
names(AlphaDivMeanSd) <- c("mean", "sd")

WriteTable (AlphaDivMeanSd, OutdirCase, "AlphaDivMeanAndSd") 

cairo_pdf((paste(OutdirCase, '/Graphs/alpha_boxplot.pdf', sep = "/")), width = 10,  height = 10)
boxplot(AlphaDivCase$AlphaDiversity)
dev.off()
#---  working with alpha diversity ---


#--- sequencing statistics--- 
SeqStatTbl <- read.table(SeqStatTblFile, comment.char = "", quote ="", as.is = T, header = T)
 
SeqStatList <- as.list(SeqStatTbl)
SeqStatMean <- lapply(SeqStatList[-match('samp_name', names(SeqStatList))], mean)
SeqStatSd <- lapply(SeqStatList[-match('samp_name', names(SeqStatList))], sd)
names(SeqStatMean) <- c("mean_init", "mean_AF", "mean_AR", "mean_Mapp")
names(SeqStatSd) <- c("sd_init", "sd_AF", "sd_AR", "sd_Mapp")
SeqStatMeanSd <- c(SeqStatMean, SeqStatSd)
SeqStatMeanSd.df <- t(as.data.frame(SeqStatMeanSd))

WriteTable (SeqStatMeanSd.df, OutdirCase, "SeqStatMeanAndSd") 
#--- end sequencing statistics---




#--- make total matrix case and control (family, genus, species) ---
#FamilyCaseControl <- UniteMatrices(FamilyCase, FamilyCtrl)
#GenusCaseControl <- UniteMatrices(GenusCase, GenusCtrl)
#SpeciesCaseControl <- UniteMatrices(SpeciesCase, SpeciesCtrl)

#--- make MDS ---

#Family
ddHSN<-bcdist(UniteMatrices(FamilyCase, FamilyCtrl))
MDS(ddHSN, GraphsDir, MetaTable, "Type.1", "Family")

#Genus
ddHSN<-bcdist(UniteMatrices(GenusCase, GenusCtrl))
MDS(ddHSN, GraphsDir, MetaTable, "Type.1", "Genus")

#Species
ddHSN<-bcdist(UniteMatrices(SpeciesCase, SpeciesCtrl))
MDS(ddHSN, GraphsDir, MetaTable, "Type.1", "Species")

#--- end ---

