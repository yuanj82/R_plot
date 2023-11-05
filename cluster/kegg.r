library(clusterProfiler)

genes <- read.delim(
    'human_genes.txt',
    header = TRUE,
    stringsAsFactors = FALSE
)[[1]]

enrich.kegg <- enrichKEGG(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = 'kegg',
    ont = 'ALL',
    pAdjustMethod = 'fdr',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = FALSE
)