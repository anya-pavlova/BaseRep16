###########################################################
### Loading case and control, alpha ad beta diversity  ###
###########################################################

#otu - operational taxonomic unit 

#--- loading case and control alpha diversity ---      
AlphaDivCase <- LoadAlphaDiv(alpha_div_case)
#--- end loading alpha diversity ---


#--- loading case beta diversity ---
FamCaseM <- FamilyCase
ddHSN <- bcdist(FamCaseM)
FamCaseMMax <- FamCaseM [, which(colMaxs(FamCaseM) > 0)]
#--- end loading beta diversity ---

#--- preparation metadate ---
#case
#meta <- cbind (1:nrow(ff), c(rep("case",length(tags_par))))
#nrow(meta)
#case+control
#meta <- cbind (1:nrow(ff), c(rep("case",length(tags_par)), rep("control",length(tags_ctrl))))
#--- end of preparation metadate ---


############## CHOOSE TAX LEVEL##############################
#sp, g, fam

TagsCaseFam <- sort(rownnames(FamilyCase))

dd <- bcdist(ff)
dd
#dd.df <- as.data.frame(dd)

mdd <- data.matrix(dd)
hr <- hclust(dd, method="ward.D") #The "ward" method has been renamed to "ward.D"; note new "ward.D2"
groups <- cutree(hr, k=2)


########################################################################################
# TEST CODE 
########################################################################################

#meta <- cbind (1:nrow(ff), c(rep("case",length(tags_par)), rep("control",length(tags_ctrl))))
meta <- cbind (1:nrow(FamCaseM), c(rep("case",length(FamCaseM))))

nrow(meta)

#--- make MDS ---
MDS(ddHSN, "/home/anna/metagenome/HSN/Graphs", meta, 2, "case")
meta[,2]
nrow(meta)
dis<-ddHSN
outdir="/home/anna/metagenome/HSN/Graphs"

meta <- cbind (1:nrow(ff_par), c(rep("case",length(tags_par))))
makeMDS(dd_par,"/home/anna/metagenome/parkinson_diseas/BR_park/graphs", meta, 2, "case")

rownames(meta)<-rownames(ff)
colSums(ff)

meta <- cbind (1:nrow(ff_par), c(rep("case",length(tags_par))))
rownames(meta)<-rownames(ff_par)

#--- end ---



#colnames(ff_top)<-str_extract(colnames(ff_top), "(?<=o__).+$")

#colnames(ff_par_top)<-str_extract(colnames(ff_par_top), "(?<=o__).+$")

colnames(ff_top) <- paste("o_", unlist(data.frame(strsplit(colnames(ff_top), "o_"))[2,]), sep="")

ff_par_top<-ff_par_top[,order(colSums(ff_par_top), decreasing = T)]


#colnames(ff_par_top) <- paste("f_", unlist(data.frame(strsplit(colnames(ff_par_top), "f_"))[2,]), sep="")

drop_unecessary_tax_names <- function(data, level="genera") {
  if (level == "genera") {
    pattern <- "g__.+"
  } else if (level == "phylum") {
    pattern <- "p__.+"
  }
  new.names <- str_extract(colnames(data), pattern)
  good.ids <- which(!is.na(new.names))
  data <- data[,good.ids]
  colnames(data) <- new.names[good.ids]
  return(data)
}

meta <- cbind (1:nrow(ff_par), c(rep("case",length(tags_par))))
rownames(meta)<-rownames(ff_par)
ff_top <- chooseTOPfeature(ff,85)
ff_par_top <- chooseTOPfeature(ff_par,85)

#colnames(ff_top)<-str_extract(colnames(ff_top), "(?<=o__).+$")
colnames(ff_par_top)<-str_extract(colnames(ff_par_top), "(?<=o__).+$")
ff_par_top<-ff_par_top[,order(colSums(ff_par_top), decreasing = T)]

colnames(ff_par_top) <- paste("g_", unlist(data.frame(strsplit(colnames(ff_par_top), "g_"))[2,]), sep="")

###### mean abundances

mean.abund.df<-as.data.frame(colMeans(ff_par))
colnames(mean.abund.df)<-"mean.percentage.abundances"
mean.abund.df$standard.deviation<-colSds(ff_par)
mean.abund.df<-mean.abund.df[order(mean.abund.df$mean.percentage.abundances, decreasing=T),]
write.table(mean.abund.df, '/home/anna/metagenome/parkinson_diseas/BR_park/out/mean_abundances.txt',quote=F,sep='\t')

##########################################################################################
# end of testing
##########################################################################################




