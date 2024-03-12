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
for (i in names(mantel)) {
  ord = pcoa(as.dist(1 - mantel[[i]]), rn = c(100, 50:99))
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

library(patchwork)
gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% names(plot_list), ]
plot_list = plot_list[gene_order$Gene]

combined_plot = wrap_plots(plot_list, ncol = 7)
# ggsave('singleM_genes/Figures/mantel.png', combined_plot, width = 25, height = 21)
