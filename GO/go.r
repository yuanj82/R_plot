library(stats) # Arch中需要先加载这些包
library(base)

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
    axis.title = element_text(face = "bold", size = 14, colour = "black"),
    axis.text = element_text(face = "bold", size = 14, colour = "black"), 
    axis.line = element_line(size = 0.5, colour = "black"), 
    panel.background = element_rect(colour = "black"),
    legend.key = element_blank()
)

p <- ggplot(
    GO_data,
    aes(
        x = Fold.Enrichment,
        y = Term, 
        colour = 1*log10(PValue),
        size = Count
    ))+
    geom_point()+
    scale_size(range = c(2, 8))+
    scale_colour_gradient(low = "blue", high = "red")+
    theme_bw()+
    ylab("Terms")+
    xlab("Fold.Enrichment")+
    labs(color = expression(-log[10](PValue)))+
    theme(legend.title = element_text(margin = margin(r = 50)), axis.title.x = element_text(margin = margin(t = 20)))+
    theme(axis.text.x = element_text(face = "bold", colour = "black", angle = 0, vjust = 1))
plot <- p + mytheme
plot

ggsave(plot, filename = "GO.png", width = 15, height = 9, dpi = 300)