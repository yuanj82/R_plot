library(ggplot2)

data <- read.table('./data.txt', header = TRUE)
data <- as.data.frame(data)

# 设置pvalue和logFC的阈值
cut_off_pvalue = 0.0000001
cut_off_logFC = 1

# 根据阈值分为上调基因down和下调基因down，无差异基因为stable，并保存到change列

data$change <- ifelse (
    data$pvalue < cut_off_pvalue & abs(data$log2FoldChange) >= cut_off_logFC, 
    ifelse(data$log2FoldChange > cut_off_logFC, 'Up', 'Down'),
    'Stable'
)

p <- ggplot(
    data,
    aes(
        x = log2FoldChange,
        y = -log10(pvalue),
        colour = change,
    ))+
    geom_point(
        alpha = 0.4,
        size = 3.5,
    )+
    scale_color_manual(values = c("#546de5", "#d2dae2", "#ff4757"))+
    geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8)+
    geom_hline(yintercept = -log10(cut_off_pvalue), lty = 4, col = "black", lwd = 0.8)+
    labs(
        x = "log2(fold change)",
        y = "-log10 (p-value)"
    )+
    theme_bw()+
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank()
    )

plot <- p
plot

ggsave(plot, filename = "volcano.png", width = 10, height = 6, dpi = 300)