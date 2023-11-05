library(org.Hs.eg.db)
# library(AnnotationHub)
library(clusterProfiler)
# library(rtracklayer)

# hub <- AnnotationHub()
# query(hub, 'Homo sapiens')

# Homo sapiens，Lolium perenne

# Human <- hub[['AH5012']]

genes <- read.delim(
    'human_genes.txt',
    header = TRUE,
    stringsAsFactors = FALSE
)[[1]]

enrich.go <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = 'ALL',
    pAdjustMethod = 'fdr',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = FALSE
)

write.table(
    enrich.go,
    'entich.go.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)

dotplot(enrich.go)