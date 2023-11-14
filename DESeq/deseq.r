library(DESeq2)

counts = read.csv(
    'gene2.counts', 
    header = T,  
    sep = '\t', 
    row.names = "Geneid", 
    comment.char = '#', 
    check.name = F
)

# 删掉前五列
counts = counts[,-c(1:5)]

# 保留行相加大于10的数据
counts = counts[rowSums(counts)>10, ]

samples = data.frame(
    sampleID = c("C0_1", "C0_2", "C0_3", "C50_1", "C50_2", "C50_3", "C500_1", "C500_2", "C500_3"), 
    sample = c("sample1", "sample1", "sample1", "sample2", "sample2", "sample2", "sample3", "sample3", "sample3")
)

# 按照sampleID更改samples的行名
rownames(samples) = samples$sampleID

# 将因子型数据的默认排序设置为1=sample1, 2=sample2, 3=sample3
samples$sample = factor(samples$sample, levels = c('sample1', 'sample2', 'sample3'))

# 构建 DESeqDataSet 对象
dds = DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~sample)

# 计算差异倍数并获得 p 值
# parallel = TRUE 可以多线程运行，在数据量较大时建议开启
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)

# 查看样本上调还是下调
sampl1_vs_sample2 <- results(dds1, contrast = c('sample', 'sample1', 'sample2'))

sampl1_vs_sample3 <- results(dds1, contrast = c('sample', 'sample2', 'sample3'))

result1 <- data.frame(sampl1_vs_sample2, stringsAsFactors = FALSE, check.names = FALSE)
result2 <- data.frame(sampl1_vs_sample3, stringsAsFactors = FALSE, check.names = FALSE)

write.table(result1, 'sampl1_vs_sample2.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)
write.table(result2, 'sampl1_vs_sample3.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)