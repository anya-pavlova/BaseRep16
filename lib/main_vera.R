rm(list=ls(all=TRUE)); gc()
################################################
# NeededSampleSize
################################################

#--- choose your flavor: ---
current.directory <- "/home/vera/Documents/metagenom/SampleSize"
setwd(current.directory)

source("lib/functions_power.R")

flog.info("start")

path.to.config <- "~/Documents/metagenom/SampleSize/src/pathway.ini"
pathway <- ReadIni(path.to.config) 

totalTable <- Load(pathway)

genus.case <- (totalTable$genus[which(rownames(totalTable$genus) 
                                                          %in% totalTable$meta[which(totalTable$meta[,"Type.1"] 
                                                                                     %in% "case"),"samples_name"]),])
genus.ctrl <- (totalTable$genus[which(rownames(totalTable$genus) 
                                        %in% totalTable$meta[which(totalTable$meta[,"Type.1"] 
                                                                   %in% "control"),"samples_name"]),])
#loading data
Case <- genus.case
Cntrl <- genus.ctrl

#собираем малопредставленные виды в один столбец
filtCase <- Data.filter(Case*6000, "data", 3, 10)
filtCntrl <- Data.filter(Cntrl*6000, "data", 3, 10)

data <- merge.features(filtCase,filtCntrl, 0)



group1 <- data[[1]]
group2 <- data[[2]]

