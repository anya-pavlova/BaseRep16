########################################################
## working with alpha, beta -diversity                ##
########################################################


#--- loading case nad control ---
TotalTable <- Load(FamCaseOtuTbl, GenCaseOtuTbl, SpeCaseOtuTbl, OtuCaseTbl, FamCtrlOtuTbl, GenCtrlOtuTbl, SpeCtrlOtuTbl, OtuCtrlTbl)

FamilyCase <- TotalTable[[1]]
GenusCase <- TotalTable[[2]]
SpeciesCase <- TotalTable[[3]]
OtuCase <- TotalTable[[4]]

FamilyCtrl <- TotalTable[[6]]
GenusCtrl <- TotalTable[[7]]
SpeciesCtrl <- TotalTable[[8]]
OtuCtrl <- TotalTable[[9]]
#--- end of loading case nad control ---

#--- TOP features ---
TotalTOPtable <- lapply(TotalTable, chooseTOPfeature, perc = 85)

FamilyTOPcase <- TotalTOPtable[[1]]
GenusTOPcase <- TotalTOPtable[[2]]
SpeciesTOPcase <- TotalTOPtable[[3]]
OtuTOPcase <- TotalTOPtable[[4]]

FamilyTOPctrl <- TotalTOPtable[[6]]
GenusTOPctrl <- TotalTOPtable[[7]]
SpeciesTOPctrl <- TotalTOPtable[[8]]
OtuTOPctrl <- TotalTOPtable[[9]]


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
#FamilyCaseControl <- UniteMatrices(family,)
#GenusCaseControl <- UniteMatrices(genus,)
#SpeciesCaseControl <- UniteMatrices(species,)
#OtuCaseControl <- UniteMatrices(otu,)
