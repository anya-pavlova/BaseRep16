########################################################
## Secondary functions (UniteMatrix, )               ##
########################################################

#-----make MDS (2D)-----
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


#-----UniteMatrices-----
#-input: two matrix  
#-output: united matrix
UniteMatrices <- function(t1, t2)
{
  uc <- sort(union(colnames(t1), colnames(t2)))
  ur <- c(rownames(t1), rownames(t2))
  t <- matrix(0, nrow = nrow(t1) + nrow(t2), ncol = length(uc))
  colnames(t) <- uc
  rownames(t) <- ur
  t[rownames(t1), colnames(t1)] <- t1
  t[rownames(t2), colnames(t2)] <- t2
  identical(t1, t[rownames(t1), colnames(t1)])
  identical(t2, t[rownames(t2), colnames(t2)])
  t
}

#ToDo: должна проверять совпадения в именах, поиск общих элементов в матрицах, стерать элементы в 1й которые есть во 2й и их добавлять

#-----WriteTable: write TOP features in .txt file-----
#TRA - table of the relative abundance
WriteTable <- function (TRA, outdir, type)
{
  filename<-paste(outdir,"/",type, '.txt', sep="")
  write.table(TRA, filename, quote=F, sep='\t')
}


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
  y<-dat[,name_list]
}
#-----end chooseTOPfeature-----
