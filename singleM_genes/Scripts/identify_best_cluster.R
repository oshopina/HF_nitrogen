library(openxlsx)

stats_df = read.xlsx('singleM_genes/Results/stats_by_cluster.xlsx')

temp = stats_df[stats_df$p_perm <= 0.01,]
temp = temp[temp$p_value_shannon <= 0.01,]

gene_tables = split(temp, temp$Gene)

best_cluster = data.frame()
for (i in names(gene_tables)) {
  df = gene_tables[[i]]
  df$sum_r2 = (2*df$R2_perm) + df$R2_shannon
  df = dplyr::top_n(df, 5, sum_r2)
  df$sum_F = (2*df$F_perm) + df$F_value_shannon
  best_cluster = rbind(best_cluster, df[df$sum_F == max(df$sum_F),])
}

gene_order = read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% best_cluster$Gene,]
rownames(best_cluster) = best_cluster$Gene
rownames(gene_order) = gene_order$Gene
best_cluster = best_cluster[rownames(gene_order),]

write.xlsx(best_cluster, 'singleM_genes/Results/best_cluster.xlsx')

length(unique(stats_df$Gene))
length(unique(stats_df$cluster))
length(unique(temp$Gene))
length(unique(temp$cluster))
 