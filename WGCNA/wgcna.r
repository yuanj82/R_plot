setwd('/home/hieroglyphs/work/R_plot/WGCNA')
rm(list = ls())  
Sys.setenv(LANGUAGE = "en")

library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse)
library(data.table)

counts = read.csv(
    'gene.counts', 
    header = T,  
    sep = '\t', 
    row.names = "Geneid", 
    comment.char = '#', 
    check.name = F
)

counts <- as.data.frame(counts)

# 删掉前五列
counts = counts[,-c(1:5)]

# 保留行相加大于10的数据
counts = counts[rowSums(counts)>10, ]

data <- log2(counts+1)

keep_data <- data[order(apply(data,1,mad), decreasing = T)[1:5000],]

datTraits <- data.frame(
  group = c(1, 2, 3, 1, 2, 3),
  row.names = c("C0_1", "C0_2", "C0_3", "C50_1", "C50_2", "C50_3")
)

datExpr0 <- as.data.frame(t(keep_data))


gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK

### 绘制样品的系统聚类树
if(T){
  #针对样本做聚类树
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)
  ## 若样本有性状、表型，可以添加对应颜色，查看是否聚类合理
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                  colors = rainbow(length(table(datTraits$group))), 
                                                   signed = FALSE)
  ## 绘制样品的系统聚类树及对应性状
  par(mar = c(1,4,3,1),cex=0.8)
  pdf("step1_Sample dendrogram and trait.pdf",width = 8,height = 6)
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait" )
  ## Plot a line to show the cut
  # abline(h = 23500, col = "red") #根据实际情况而定
  dev.off()
}

##若存在显著离群点；剔除掉
if(F){
  clust <- cutreeStatic(sampleTree, cutHeight = 23500, minSize = 10) # cutHeight根据实际情况而定
  table(clust)
  keepSamples <- (clust==1)
  datExpr0 <- datExpr0[keepSamples, ]
  datTraits <- datTraits[keepSamples,]
  dim(datExpr0) 
}

### 判断数据质量 : PCA进行分组查看
rm(list = ls())  
load("step1_input.Rdata")
group_list <- datTraits$group
dat.pca <- PCA(datExpr0, graph = F) 
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point","text"), #"point","text"
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE, #标签不重叠
                    col.ind = group_list, # 分组上色
                    axes.linetype=NA,  # remove axeslines
                    mean.point=F#去除分组中心点
                    ) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1) #坐标轴的纵横比
pca