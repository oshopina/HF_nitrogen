
stats_df = openxlsx::read.xlsx('singleM_genes/Results/stats_by_cluster.xlsx')

temp = stats_df[stats_df$p_perm <= 0.01,]
temp = temp[temp$p_value_sobs <= 0.01,]
temp = temp[temp$p_value_shannon <= 0.01,]
temp = temp[temp$p_value_chao1 <= 0.05,]

gene_tables = split(temp, temp$Gene)

best_cluster = data.frame()
for (i in names(gene_tables)) {
  df = gene_tables[[i]]
  df$sum_r2 = df$R2_perm + df$R2_sobs + df$R2_shannon + df$R2_chao1
  df = dplyr::top_n(df, 5, sum_r2)
  df$sum_F = df$F_perm + df$F_value_sobs + df$F_value_shannon + df$F_value_chao1
  best_cluster = rbind(best_cluster, df[df$sum_F == max(df$sum_F),])
}

# openxlsx::write.xlsx(best_cluster, 'singleM_genes/Results/best_cluster.xlsx')

length(unique(stats_df$Gene))
length(unique(stats_df$cluster))
length(unique(temp$Gene))
length(unique(temp$cluster))
 