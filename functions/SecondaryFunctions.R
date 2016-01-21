########################################################
## Secondary functions (UniteMatrix, )               ##
########################################################

#-----make MDS (2D)-----
#-input: distance matrix, output directory, meta matrix,  
#-output: MDS plot
MDS<-function(dis, outdir, meta, colFact, plot_name,m_width = 5, m_height=5)
{   
  myMDS <- isoMDS(dis, k=2)
  outdirMDS<-paste(outdir,'/MDS_',plot_name,'.pdf',sep='')
  outdirMDS<-file.path(outdirMDS)
  dfMDS <- data.frame(X = as.vector(myMDS$points[,1]), Y = as.vector(myMDS$points[,2]),cols = meta[,colFact], lab = rownames(myMDS$points))
  plotMDS <- ggplot(dfMDS, aes(x=X, y=Y, color=factor(cols), label=lab)) + geom_point(size=0.3) + geom_text(size=2) + theme_bw()    
  ggsave(outdirMDS, plotMDS, width = m_width, height=m_height)  
}
#-----end MDS-----


#ToDo: должна проверять совпадения в именах, поиск общих элементов в матрицах, стерать элементы в 1й которые есть во 2й и их добавлять
#-----WriteTable: write TOP features in .txt file-----
#TRA - table of the relative abundance
WriteTable <- function (SomeStructure, outdir, type)
{
  filename<-paste(outdir,"/",type, '.txt', sep="")
  write.table(SomeStructure, filename, quote=F, sep='\t')
}
#-----end WriteTable: write something in .txt file-----

#-----chooseTOPfeature-----
#-choose top features with total % of abundance across all samples
#-input: feature vectors (table), percantage
#-output: table with chosen features
chooseTOPfeature<-function(dat,perc)
{
  summof<-sum(colSums(dat))
  sorted<-sort(colSums(dat), decreasing = TRUE)
  count <-0
  num<-0
  name_list<-c()  
  for (i in sorted) #in sorted by total abundance features, take those making %
  {
    count<-count+((i/summof)*100)
    num <-num+1
    if (count>perc) 
    {
      break
    }
    else 
    {
      name_list<-append(name_list, names(sorted[num]))
    }  	
  }	
  y <- dat[,name_list]
  return(y)
}
#-----end chooseTOPfeature-----

###make  distance in Spearman correlation  - for internal use
####input: feature vectors (table)
####output: class dist object
distSpear<-function (x) 
{
  result <- 1-cor(t(x), method='spearman', use='pair')
  as.dist(result)
}
#####end

####good colours for heatplot
cols.gentleman <- function(ncol=500) {
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(10, 'RdBu'))(ncol)
  return(rev(hmcol))
}
#####end

