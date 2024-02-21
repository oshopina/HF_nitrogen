library(gsveasyr)
library(ComplexHeatmap)

env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
load_alpha_diversity('singleM_genes/Data/', env.df = env,
                     changeSampleName = T, refColumn = env$gsveasy_sample, 
                     clustlvls = 50:100)

metrics = c('C_sobs', 'C_shannon', 'C_chao1')
alpha_clusters = mget(ls(pattern = 'GSV'))
rm(list = ls(pattern = 'GSV'))

alpha_results = list()
for (v in metrics) {
  fs = c()
  ps = c()
  alpha_r2 = data.frame()
  for (j in names(alpha_clusters)) {
    anova = auto_aov_fixed(alpha_clusters[[j]][[v]], ~ pH, alpha_clusters[[j]]$C.env)
    gene_names = unique(anova$Results$Data)
    fs = anova$Results$F_value %>% na.omit()
    ps = anova$Results$p_value %>%  na.omit()
    r2s = c()
    for (i in gene_names) {
      r2 = summary(lm(
        alpha_clusters[[j]][[v]][[i]] ~ pH,
        alpha_clusters[[j]]$C.env
      ))$r.squared
      r2s[i] = r2
    }
    alpha_r2 = rbind(alpha_r2, cbind(
      Cluster = j,
      R2 = r2s,
      F_value = fs,
      p_value = ps,
      Gene = names(r2s)
    ))
  }
  alpha_results[[v]] = alpha_r2
}

heatmaps = list()
for (i in c('R2', 'F_value')) {
  heatmap_table = pivot_wider(alpha_results$C_shannon[,c('Cluster', i, 'Gene')],
                              names_from = Cluster,
                              values_from = !!sym(i)) %>% as.data.frame()
  rownames(heatmap_table) = heatmap_table$Gene
  heatmap_table = heatmap_table[,-1]
  heatmap_table[] = lapply(heatmap_table, as.numeric)
  heatmap_table = scale(t(heatmap_table))
  
  col_fun = colorRampPalette(c('white', 'yellow', 'brown'))
  gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
  gene_order = gene_order[gene_order$Gene %in% colnames(heatmap_table),]
  
  plot = Heatmap(t(heatmap_table), cluster_rows = F, cluster_columns = F, row_order = gene_order$Gene,
                 column_order = c(names(alpha_clusters[2:51]), names(alpha_clusters[1])),
                 column_labels = c(100,50:99),
                 name = i, column_title = 'Shannon diversity', col = col_fun(100))
  heatmaps[i] = plot
}


        