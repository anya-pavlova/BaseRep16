########################################################
## MDS, HeatMaps, BoxPlot                             ##
########################################################




#meta <- cbind (1:nrow(ff), c(rep("case",length(tags_par)), rep("control",length(tags_ctrl))))
meta <- cbind (1:nrow(ff), c(rep("case",length(tags_par))))
nrow(meta)


makeMDS(dd, "/home/anna/metagenome/parkinson_diseas/BR_park/graphs", meta, 2, "case")
meta[,2]
nrow(meta)
dis<-dd
outdir="/home/anna/metagenome/parkinson_diseas/BR_park/graphs"
#length(as.vector(myMDS$points[,1]))

#nrow(sp_ctrl)
#nrow(sp_all)

## just case



meta <- cbind (1:nrow(ff_par), c(rep("case",length(tags_par))))
makeMDS(dd_par,"/home/anna/metagenome/parkinson_diseas/BR_park/graphs", meta, 2, "case")
