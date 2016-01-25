########################################################
## MDS, HeatMaps, BoxPlot                             ##
########################################################
#####################################################
# produce heatmaps - for selected pairs - from the previously selected
#####################################################

# family 

ff_par_top <- FamilyTOPcase[, which(colMaxs(FamilyTOPcase) > 3)]
ff_par_top <- ff_par_top[, order(colMaxs(ff_par_top), decreasing = F)]
dd <- bcdist(FamilyCase)
hr <- hclust(dd, method="ward.D")
hc <- hclust(distSpear(t(ff_par_top)), method="ward.D")
# final pre-print change names of heatmap rows and cols
colnames(ff_par_top) <- paste("o_", unlist(data.frame(strsplit(colnames(ff_par_top), "o_"))[2,]), sep="")

# draw
cairo_pdf(paste("Graphs/heat_fam.pdf", sep=""), width=23, height=18)

#cairo_pdf(paste("graphs/heat_fam.pdf", sep=""), width=15, height=10)
heatmap.2(data.matrix(ff_par_top), Rowv=as.dendrogram(hr), 
          dendrogram="row",
          col=cols.gentleman(500), 
          cexRow=1.2, cexCol=1.2, trace='none', scale="none", margin=c(25,20), )
dev.off()

# genus

ff_par_top <- GenusTOPcase[, which(colMaxs(GenusTOPcase) > 3)]
ff_par_top <- ff_par_top[, order(colMaxs(ff_par_top), decreasing = F)]
dd <- bcdist(ff_par_top)
#dd <- bcdist(org[rownames(ff_top),]) # !
hr <- hclust(dd, method="ward.D")
hc <- hclust(distSpear(t(ff_par_top)), method="ward.D")
# final pre-print change names of heatmap rows and cols
colnames(ff_par_top) <- paste("o_", unlist(data.frame(strsplit(colnames(ff_par_top), "o_"))[2,]), sep="")
# draw
cairo_pdf(paste("Graphs/heat_genus.pdf", sep=""), width=23, height=18)
#cairo_pdf(paste("graphs/heat_g.pdf", sep=""), width=15, height=10)
heatmap.2(data.matrix(ff_par_top), Rowv=as.dendrogram(hr), 
          #Colv=as.dendrogram(hc), dendrogram="both", #Colv=as.dendrogram(hc), 
          dendrogram="row",
          #col=gray(seq(0.1, 0.99, 0.001)), #cols.yuppie(500), 
          col=cols.gentleman(500), 
          #breaks=c(0, 0.1, seq(0.2, 0.9, (0.9-0.2)/89)),
          cexRow=1.2, cexCol=1.2, trace='none', scale="none", margin=c(35,30), )
#cexRow=1.2, cexCol=1.2, trace='none', scale="none", margin=c(32,10), )
dev.off()

# species

ff_par_top <- SpeciesTOPcase[, which(colMaxs(SpeciesTOPcase) > 3)]
ff_par_top <- ff_par_top[, order(colMaxs(ff_par_top), decreasing = F)]
#dd <- bcdist(org[rownames(ff_par_top),]) # !
hr <- hclust(dd, method="ward.D")
hc <- hclust(distSpear(t(ff_par_top)), method="ward.D")
# final pre-print change names of heatmap rows and cols
colnames(ff_par_top) <- paste("o_", unlist(data.frame(strsplit(colnames(ff_par_top), "o_"))[2,]), sep="")
# draw
cairo_pdf(paste("Graphs/heat_species.pdf", sep=""), width=23, height=18)
#cairo_pdf(paste("graphs/heat_sp.pdf", sep=""), width=15, height=10)
heatmap.2(data.matrix(ff_par_top), Rowv=as.dendrogram(hr), 
          #Colv=as.dendrogram(hc), dendrogram="both", #Colv=as.dendrogram(hc), 
          dendrogram="row",
          #col=gray(seq(0.1, 0.99, 0.001)), #cols.yuppie(500), 
          col=cols.gentleman(500), 
          #breaks=c(0, 0.1, seq(0.2, 0.9, (0.9-0.2)/89)),
          cexRow=1.2, cexCol=1.2, trace='none', scale="none", margin=c(40,35), )
#cexRow=1.2, cexCol=1.2, trace='none', scale="none", margin=c(32,10), )
dev.off()




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
