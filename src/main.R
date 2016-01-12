
rm(list=ls(all=TRUE)); gc()
################################################
# Project Parkinson_diseas - 16S rRNA metagenomes
################################################

# choose your flavor:
cur <- "~/metagenome/parkinson_diseas/BR_park/"
#cur <- ""
setwd(cur)

source("lib/functions.R")
source("lib/functions_special.R")
source("lib/functions_special2.R")
source("lib/loaders.R")







# TODO: read register of all samples and prepare some meta-data variables and files
#source("src/.R")





write.table(fam, '/home/anna/metagenome/parkinson_diseas/BR_park/out/fam.txt',quote=F,sep='\t')
write.table(g, '/home/anna/metagenome/parkinson_diseas/BR_park/out/g.txt',quote=F,sep='\t')
write.table(sp, '/home/anna/metagenome/parkinson_diseas/BR_park/out/sp.txt',quote=F,sep='\t')
write.table(otu, '/home/anna/metagenome/parkinson_diseas/BR_park/out/otu.txt',quote=F,sep='\t')








# select samples, drop some and prepare the tags
#source("src/preprocess.R")

# compare with control
#source("src/diff_abund.R")

## visualization
source("src/viz.R")





