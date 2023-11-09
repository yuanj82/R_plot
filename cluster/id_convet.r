library("org.Hs.eg.db") 

g_symbol <- read.delim(
    'ids.txt',
    header = TRUE,
    stringsAsFactors = FALSE
)[[1]]

g_id = mapIds(x = org.Hs.eg.db,
              keys = g_symbol,
              keytype = "SYMBOL",
              column = "ENTREZID")

write.table(
    g_id,
    'g_id.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)