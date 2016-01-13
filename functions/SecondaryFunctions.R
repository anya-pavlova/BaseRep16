########################################################
## Secondary functions (UniteMatrix, )               ##
########################################################


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