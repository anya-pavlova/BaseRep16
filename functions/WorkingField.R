########################################################
## working with alpha, beta -diversity                ##
########################################################

#--- TOP features ---
#FamGenSpeList <- list (TotalTable$Family, TotalTable$Genus, TotalTable$Species)

TopFeatures <- function(FeatList)
{
  TopFeatures <- lapply(FamGenSpeList, chooseTOPfeature, perc = 85)
  
}

FamGenSpeListTOP <- TopFeatures(FamGenSpeList)
Family <- as.data.frame(FamGenSpeListTOP)

TopFeatures <- function (TotalTable, type, subset, perc)
{
  FamilyTOP <- chooseTOPfeature(TotalTable$Family[which(rownames(TotalTable$Family)%in%TotalTable$Meta[which(TotalTable$Meta[,type]%in%subset),"samples_name"]),], perc=perc)
  GenusTOP <- chooseTOPfeature(TotalTable$Genus[which(rownames(TotalTable$Genus)%in%TotalTable$Meta[which(TotalTable$Meta[,type]%in%subset),"samples_name"]),], perc=perc)
  SpeciesTOP<-chooseTOPfeature(TotalTable$Species[which(rownames(TotalTable$Species)%in%TotalTable$Meta[which(TotalTable$Meta[,type]%in%subset),"samples_name"]),], perc=perc)
  OtuTOP <- chooseTOPfeature(TotalTable$Otu[which(rownames(TotalTable$Otu)%in%TotalTable$Meta[which(TotalTable$Meta[,type]%in%subset),"samples_name"]),], perc=perc)


  list(FamilyTOP = FamilyTOP, GenusTOP = GenusTOP, SpeciesTOP = SpeciesTOP, OtuTOP = OtuTOP)
}

TableTOP<-TopFeatures(TotalTable, type= "Type.1", subset="case", perc=85)


FamilyTOP <- chooseTOPfeature(TotalTable$Family[which(rownames(TotalTable$Family)%in%),2] , perc=85)
GenusTOP <- chooseTOPfeature(TotalTable$Genus , perc=85)
SpeciesTOP <- chooseTOPfeature(TotalTable$Species , perc=85)
OtuTOP <- chooseTOPfeature(TotalTable$Otu , perc=85)

WriteTable (FamilyCase, OutdirCase, "FamilyCase") 
WriteTable (GenusCase, OutdirCase, "GenusCase") 
WriteTable (SpeciesCase, OutdirCase, "SpeciesCase") 
WriteTable (OtuCase, OutdirCase, "OtuCase") 
#--- end write TOP features in .txt file ---

#--- working with alpha diversity ---
                                  
AlphaDivCase <- (TotalTable$AlphaDiv[which(rownames(TotalTable$AlphaDiv)%in%TotalTable$Meta[which(TotalTable$Meta[,"Type.1"]%in%"case"),"samples_name"]),])]
               
WriteTable( AlphaDivCase, PathWay$Case$OutdirCase, "AlphaDivCaseTbl"))

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
ddHSN<-bcdist(TotalTable$Family)
MDS(ddHSN, PathWay$Case$GraphsDir, PathWay$Case$MetaTable, "Type.1", "Family")

#Genus
ddHSN<-bcdist(UniteMatrices(TotalTable$Genus, GenusCtrl))
MDS(ddHSN, PathWay$Case$GraphsDir, PathWay$Case$MetaTable, "Type.1", "Genus")

#Species
ddHSN<-bcdist(UniteMatrices(SpeciesCase, SpeciesCtrl))
MDS(ddHSN, PathWay$Case$GraphsDir , TotalTable$MetaTable, "Type.1", "Species")

#--- end ---

