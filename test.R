mlist <- replicate(30, matrix(sample(1e5), nrow = 1e3), simplify = F)

for (i in 1:length(mlist)) {
  rownames(mlist[[i]]) <- paste("sample_", i, "_", 1:nrow(mlist[[i]]), sep = "")  
}

unit2matr <- function(m1, m2) {
  rbind(m1, m2)
}

system.time({
mm <- mlist[[1]]
for (i in 2:length(mlist)) {
  mm <- unit2matr(mm, mlist[[i]])
}
})

system.time({
mm <- matrix(0, nrow = length(mlist) * 1e3, ncol = 1e2)
rownames(mm) <- unlist(lapply(mlist, rownames))
for(i in 1:length(mlist)) {
  mm[rownames(mlist[[i]]), ] <- mlist[[i]]
}
})

head(mm)
