

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
