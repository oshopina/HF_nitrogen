library(ape)
library(ggplot2)
library(ggrepel)

mantel = readRDS('singleM_genes/Results/pr_mat.rds')

mantel = lapply(mantel, function(x) {
  x[x == 1] = 0
  return(x)
})

mantel = lapply(mantel, function(x) {
  x[is.na(x)] = 1
  return(x)
})

names(mantel)[which(names(mantel) == "narG")] = 'narG_nxrA'
names(mantel)[which(names(mantel) == "narH")] = 'narH_nxrB'

library(ggcorrplot)
cor_plots = list()
for (i in names(mantel)) {
  df = mantel[[i]]
  df = 1 - df
  cor_matrix = 2 * df - 1
  rownames(cor_matrix) = c(100, 50:99)
  colnames(cor_matrix) = rownames(cor_matrix)
  
  p = ggcorrplot(cor_matrix) +
    theme(legend.position = 'none') +
    ggtitle(i)
  
  cor_plots[[i]] = p
}

gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% names(cor_plots), ]
cor_plots = cor_plots[gene_order$Gene]

library(cowplot)
combined_plot2 = plot_grid(plotlist = cor_plots, ncol = 7)
ggsave('singleM_genes/Figures/protest_cor.png', combined_plot2, width = 27, height = 21)
