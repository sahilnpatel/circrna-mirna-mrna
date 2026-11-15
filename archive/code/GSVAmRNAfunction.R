genesets = read.delim(file.choose())
library(GSVA)
gsvaresults = gsva(as.matrix(data), genesets)
write.table(gsvaresults, sep = "\t", "gsvaresults.txt")

