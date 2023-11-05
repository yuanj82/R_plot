

data <- read.table('./data.txt', header = TRUE)
data <- as.data.frame(data)

cut_off_pvalue = 0.0000001
cut_off_logFC = 1

library(VennDiagram)

set1 <- log2FoldChange