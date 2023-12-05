# 数据示例：
# Term	Fold Enrichment		PValue	Count
# ribonucleoprotein complex biogenesis	0.0644 		1.43E-74	310
# mRNA processing	0.0619 		1.19E-58	298
# ribosome biogenesis	0.0428 		5.89E-56	206
# Golgi vesicle transport	0.0409 		3.24E-50	197
# RNA splicing	0.0557 		1.48E-48	268
# cytoplasmic translation	0.0264 		2.50E-46	127
# rRNA processing	0.0322 		8.50E-45	155
# establishment of protein localization to organelle	0.0532 		1.67E-44	256
# rRNA metabolic process	0.0353 		2.61E-42	170
# focal adhesion	0.0617 		2.26E-92	304
# cell-substrate junction	0.0619 		4.74E-89	305

setwd('/home/hieroglyphs/work/R_plot/DESeq')
rm(list = ls())  
Sys.setenv(LANGUAGE = "en")

library(stats) # Arch中需要先加载这些包
library(base)
library(showtext)
font_add("Times_New_Roman", "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf")  # 添加新罗马字体
showtext_auto()
# Attaching package: ‘dplyr’

# The following objects are masked from ‘package:stats’:

#     filter, lag

# The following objects are masked from ‘package:base’:

#     intersect, setdiff, setequal, union


library(dplyr)
library(ggplot2)
library(ggrepel)

GO_data <- read.table('./data.txt', header = TRUE,  sep = '\t')

GO_data <- arrange(GO_data, GO_data[,3])

GO_data$Term <- factor(GO_data$Term, levels = rev(GO_data$Term))

mytheme <- theme(
    axis.title = element_text(family = "Times_New_Roman", face = "bold", size = 25, colour = "black"),
    axis.text = element_text(family = "Times_New_Roman", face = "bold", size = 25, colour = "black"), 
    axis.line = element_line(size = 0.2, colour = "black"), 
    panel.background = element_rect(colour = "black"),
    legend.key = element_blank(),
    legend.title = element_text(family = "Times_New_Roman", face = "bold", size = 30),
    legend.text = element_text(family = "Times_New_Roman", face = "bold", size = 25),
    # legend.spacing.y = unit(0.5, "cm")  # 调整标尺与图形的距离
)

p <- ggplot(
    GO_data,
    aes(
        x = Fold.Enrichment,
        y = Term, 
        colour = 1*log10(PValue),
        size = Counts
    ))+
    geom_point()+
    scale_size(range = c(1, 4))+
    scale_colour_gradient(low = "blue", high = "red")+
    theme_bw()+
    ylab("Terms")+
    xlab("Fold.Enrichment")+
    labs(color = expression(-log[10](PValue)))+
    theme(legend.title = element_text(margin = margin(r = 30)), axis.title.x = element_text(margin = margin(t = 20)))+   # 图例与图像界面右边距
    theme(axis.text.x = element_text(family = "Times_New_Roman", face = "bold", colour = "black", angle = 0, vjust = 1))
plot <- p + mytheme
plot

ggsave(plot, filename = "bubble.png", width = 7, height = 7, dpi = 300)