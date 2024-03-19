library(gsveasyr)
library(ComplexHeatmap)

# env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
# rownames(env) = env$Hoosfield.ID
# load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
#                         changeSampleName = T, refColumn = env$gsveasy_sample, 
#                         clustlvls = c(50:100))
# clusters = mget(ls(pattern = 'GSV'))
# rm(list = ls(pattern = 'GSV'))
# genes = unique(names(clusters$GSV_rarefied_gsvtables100))
# 
# hellinger_diversity = function(otu_table) {
#   dist_matrix = dist(vegan::decostand(otu_table, method = 'hellinger'),
#                      method = 'euc') %>% as.matrix() 
#   return(dist_matrix)
# }
# 
# cluster_table = data.frame()
# r2 = c()
# p = c()
# 
# for (j in names(clusters)) {
#   print(j)
#   for (i in genes) {
#     print(i)
#     df = clusters[[j]][[i]]
#     df = df[rownames(df) %in% rownames(env), ]
#     ph_table = env[rownames(df), ]
#     
#     distance = hellinger_diversity(df)
#     PC = wcmdscale(distance)
#     
#     if (any(distance != 0)) {
#       fit = envfit(PC, ph_table$pH)
#       r2[i] = fit$vectors$r
#       p[i] = fit$vectors$pvals
#     } else {
#       r2[i] = 0
#       p[i] = 1
#     }
#     
#     }
#     cluster_table = rbind(cluster_table, cbind(
#       Cluster = j,
#       R2 = r2,
#       p_value = p,
#       Gene = names(r2)
#     ))
# }
# 
# save.image('../HF_nitrogen_paper/envfit.RData')
## write.csv2(cluster_table, 'singleM_genes/Results/envfit.csv')

load('../HF_nitrogen_paper/envfit.RData')
cluster_table$Gene[cluster_table$Gene == 'narG'] = 'narG_nxrA'
cluster_table$Gene[cluster_table$Gene == 'narH'] = 'narH_nxrB'

heatmaps = list()
for (i in colnames(cluster_table[,2:3])) {
  heatmap_table = pivot_wider(cluster_table[,c('Cluster', i, 'Gene')], names_from = Cluster, values_from = sym(i)) %>% as.data.frame()
  rownames(heatmap_table) = heatmap_table$Gene
  heatmap_table = heatmap_table[, -1]
  heatmap_table[] = lapply(heatmap_table, as.numeric)
  
    heatmap_table = scale(t(heatmap_table))
  
  col_fun = colorRampPalette(c('white', 'yellow', 'brown'))
  gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
  gene_order = gene_order[gene_order$Gene %in% unique(cluster_table$Gene),]
  
  plot = Heatmap(t(heatmap_table), cluster_rows = F, cluster_columns = F, row_order = gene_order$Gene,
                 column_order = c(names(clusters[2:51]), names(clusters[1])),
                 column_labels = c(100,50:99),
                 name = i, column_title = 'envfit', col = col_fun(100))
  heatmaps[i] = plot
}

