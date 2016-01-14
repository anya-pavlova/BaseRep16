########################################################
## working with alpha, beta -diversity                ##
########################################################


#--- loading case nad control ---
TotalTable <- Load(FamCaseOtuTbl, GenCaseOtuTbl, SpeCaseOtuTbl, OtuCaseTbl, FamCtrlOtuTbl, GenCtrlOtuTbl, SpeCtrlOtuTbl, OtuCtrlTbl)

FamilyCase <- TotalTable[[1]]
GenusCase <- TotalTable[[2]]
SpeciesCase <- TotalTable[[3]]
OtuCase <- TotalTable[[4]]

#FamilyCtrl <- TotalTable[[5]]
#GenusCtrl <- TotalTable[[6]]
#SpeciesCtrl <- TotalTable[[7]]
#OtuCtrl <- TotalTable[[8]]
#--- end of loading case nad control ---

#--- write TOP features in .txt file ---
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
WriteTable (mapply(c, AlphaDivMean, AlphaDivSd, SIMPLIFY=FALSE), OutdirCase, "AlphaDivMeanAndSd") 

cairo_pdf((paste(OutdirCase, '/Graphs/alpha_boxplot.pdf', sep = "/")), width = 10,  height = 10)
boxplot(AlphaDivCase$AlphaDiversity)
dev.off()
#---  working with alpha diversity ---


#--- sequencing statistics--- 
SeqStatTbl <- read.table(SeqStatTblFile, comment.char = "", quote ="", as.is = T, header = T)
 
SeqStatList <- as.list(SeqStatTbl)
SeqStatMean <- lapply(SeqStatList[-match('samp_name', names(SeqStatList))], mean)
SeqStatSd <- lapply(SeqStatList[-match('samp_name', names(SeqStatList))], sd)

WriteTable (mapply(c, SeqStatMean, SeqStatSd, SIMPLIFY=FALSE), OutdirCase, "SeqStatMeanAndSd") 
#--- end sequencing statistics---


#--- make total matrix case and control (family, genus, species) ---
#FamilyCaseControl <- UniteMatrices(family,)
#GenusCaseControl <- UniteMatrices(genus,)
#SpeciesCaseControl <- UniteMatrices(species,)
#OtuCaseControl <- UniteMatrices(otu,)
