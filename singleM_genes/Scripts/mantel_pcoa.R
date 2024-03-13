library(ape)
library(ggplot2)
library(ggrepel)

mantel = readRDS('singleM_genes/Results/mantel_mat.rds')

mantel = lapply(mantel, function(x) {
  x[is.na(x)] = 1
  return(x)
})

names(mantel)[which(names(mantel) == "narG")] = 'narG_nxrA'
names(mantel)[which(names(mantel) == "narH")] = 'narH_nxrB'

plot_list = list()
pc1 = matrix(nrow = 51, ncol = 45) |> as.data.frame()
colnames(pc1) = names(mantel)
rownames(pc1) = c(100, 50:99)
for (i in names(mantel)) {
  ord = pcoa(as.dist(1 - mantel[[i]]), rn = c(100, 50:99))
  pc1[[i]] = ord$vectors[,'Axis.1']
  plot_df = as.data.frame(ord$vectors)
  plot_df$Cluster = as.numeric(rownames(plot_df))
  
  p = ggplot(plot_df, aes(x = Axis.1, y = Axis.2, color = Cluster)) +
    geom_point(size = 2) +
    labs(title = i) +
    theme_minimal(base_size = 15) +
    scale_color_gradient(low = 'blue', high = 'red') +
    xlab(paste0("PC1: ", (ord$values$Relative_eig[1] * 100) |> signif(digits = 2), '%')) + #add variance explained to axis
    ylab(paste0("PC2: ", (ord$values$Relative_eig[2] * 100) |> signif(digits = 2), '%')) +
    xlim(-0.4, 0.4) +
    ylim(-0.4, 0.4) +
    theme(legend.position = 'none')
  
  plot_list[[i]] = p
}

cor_values = c()
for (i in colnames(pc1)) {
  cor_value = cor(pc1[,i], as.numeric(rownames(pc1)))
  cor_values[i] = cor_value
}

table = gt::gt(data.frame(Gene = names(cor_values), Cor = cor_values))
gt::gtsave(table, 'singleM_genes/Results/pc1_vs_cluster.docx')

library(patchwork)
gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% names(plot_list), ]
plot_list = plot_list[gene_order$Gene]

combined_plot = wrap_plots(plot_list, ncol = 7)
# ggsave('singleM_genes/Figures/mantel.png', combined_plot, width = 25, height = 21)

library(ggcorrplot)
cor_plots = list()
for (i in names(mantel)) {
  df = mantel[[i]]
  cor_matrix = 2 * df - 1
  rownames(cor_matrix) = c(100, 50:99)
  colnames(cor_matrix) = rownames(cor_matrix)
  
  p = ggcorrplot(cor_matrix) +
    theme(legend.position = 'none') +
    ggtitle(i)
  
  cor_plots[[i]] = p
}

cor_plots = cor_plots[gene_order$Gene]

library(cowplot)
combined_plot2 = plot_grid(plotlist = cor_plots, ncol = 7)
ggsave('singleM_genes/Figures/mantel_cor.png', combined_plot2, width = 27, height = 21)
