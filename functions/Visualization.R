########################################################
## MDS, HeatMaps, BoxPlot                             ##
########################################################

#-----make MDS for case (2D)-----
#-input: distance matrix, output directory, meta matrix,  
#-output: MDS plot

MDS<-function(dis, outdir,meta, colFact,plot_name,m_width = 5, m_height=5)
{   
  myMDS <- isoMDS(dis, k=2)
  outdirMDS<-paste(outdir,'/MDS_',plot_name,'.pdf',sep='')
  outdirMDS<-file.path(outdirMDS)
  dfMDS <- data.frame(X = as.vector(myMDS$points[,1]), Y = as.vector(myMDS$points[,2]),cols = meta[,colFact], lab = rownames(myMDS$points))
  plotMDS <- ggplot(dfMDS, aes(x=X, y=Y, color=factor(cols), label=lab)) + geom_point(size=0.3) + geom_text(size=2) + theme_bw()  	
  ggsave(outdirMDS, plotMDS, width = m_width, height=m_height)  
}
#-----end MDS-----




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
