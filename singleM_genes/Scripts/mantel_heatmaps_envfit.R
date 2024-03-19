library(ape)
library(ComplexHeatmap)

mantel = readRDS('singleM_genes/Results/mantel_mat.rds')

mantel = lapply(mantel, function(x) {
  x[is.na(x)] = 1
  return(x)
})

names(mantel)[which(names(mantel) == "narG")] = 'narG_nxrA'
names(mantel)[which(names(mantel) == "narH")] = 'narH_nxrB'
gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% names(mantel), ]
mantel = mantel[gene_order$Gene]

## Add number of GSVs

GSVs = read.csv2('singleM_genes/Data/GSV_number.csv', row.names = 1)
names(GSVs)[which(names(GSVs) == "narG")] = 'narG_nxrA'
names(GSVs)[which(names(GSVs) == "narH")] = 'narH_nxrB'
GSVs$numbers = c(100, 50:99)
GSVs = GSVs[order(GSVs$numbers),]
rownames(GSVs)= GSVs$numbers
GSVs = GSVs[,-ncol(GSVs)]

## Add rarefracrtion levels and PERMANOVA results

r_levels = read.csv2('singleM_genes/Data/rarefraction_levels.csv', row.names = 1)
sperm = read.csv2('singleM_genes/Results/envfit.csv', row.names = 1)
sperm$Gene[which(sperm$Gene == "narG")] = 'narG_nxrA'
sperm$Gene[which(sperm$Gene == "narH")] = 'narH_nxrB'

heatmaps = list()
for (i in names(mantel)) {
  ## Create Heatmap
  my_palette <- colorRampPalette(c('blue','white', 'red'))
  
  df = mantel[[i]] |> as.data.frame()
  df$numbers = c(100, 50:99)
  df = df[order(df$numbers),]
  df = df[, c(rownames(df), 'numbers')]
  rownames(df)= df$numbers
  df = df[,-ncol(df)]
  colnames(df) = rownames(df)
  
  temp = seq(50, 100, by = 10)
  labels = rep('', length(rownames(df)))
  labels[match(temp, rownames(df))] = temp
  
  ## Add clustering
  clust = hclust(as.dist(1 - df))
  clust = cutree(clust, k = 5)
  clust = as.factor(clust)
  names(clust) = paste('c.', names(clust))
  
  raref = r_levels[i,]
  perm = sperm[sperm$Gene == i,]
  perm$numbers = c(100, 50:99) |> as.numeric()
  perm = perm[order(perm$numbers),]
  perm$asterisk <- ifelse(perm$p_value < 0.05, "*", "")
  perm$R2 = as.numeric(perm$R2)
  perm$p_value = as.numeric(perm$p_value)
  
  ha = rowAnnotation(Group = clust,
                     col = list(Group = c(`1` = '#7fc97f', `2` = '#beaed4', `3` = '#fdc086', `4` = '#386cb0', `5` = '#ffff99')),
                     gp = gpar(col = "black"),
                     show_legend = F)
  
  ha1 = HeatmapAnnotation(`GSV number` = anno_barplot(GSVs[,i]),
                          height = unit(1, 'cm'))
  
  col_fun = circlize::colorRamp2(c(min(perm$R2), max(perm$R2)), c('yellow', 'brown'))
  ha2 = HeatmapAnnotation(R2 = perm$R2, 
                          col = list(R2 = col_fun),
                          p = anno_text(perm$asterisk),
                          annotation_name_side = 'left',
                          show_legend = F)
  
  hm = Heatmap(as.matrix(df),
          row_order = rownames(df),
          column_order = colnames(df),
          row_labels = labels,
          column_labels = labels,
          right_annotation = ha,
          top_annotation = ha1,
          bottom_annotation = ha2,
          show_heatmap_legend = F,
          column_title = paste0(i, '(', raref, ')'))
  
  heatmaps[[i]] = hm
}

library(ggplotify)
heatmaps = lapply(heatmaps, function(x){
  as.ggplot(x)
})

library(cowplot)
combined_plot2 = plot_grid(plotlist = heatmaps, ncol = 7)
## ggsave('singleM_genes/Figures/mantel_heatmap_envfit.png', combined_plot2, width = 30, height = 27)

# pdf('singleM_genes/Figures/mantel_heatmaps.pdf')
# for (i in seq_along(heatmaps)) {
#   print(heatmaps[[i]])
# }
# dev.off()
