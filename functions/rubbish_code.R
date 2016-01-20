
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
