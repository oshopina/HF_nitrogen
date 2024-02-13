library(gsveasyr)

env = read.csv('graftM_genes/Data/mag_env_for_shotgun_samples.csv')
rownames(env) = env$Hoosfield.ID
load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
                     changeSampleName = T, refColumn = env$gsveasy_sample, 
                     clustlvls = c(50:100))
clusters = mget(ls(pattern = 'GSV'))
rm(list = ls(pattern = 'GSV'))
genes = unique(names(clusters$GSV_rarefied_gsvtables100))

cluster_table = data.frame()
for (j in names(clusters)) {
  print(paste('Started cluster:', j))
  genes_df = data.frame()
  for (i in genes) {
    print(paste('Started gene:', i))
    df = clusters[[j]][[i]]
    df = df[rownames(df) %in% rownames(env),]
    ph_table = env[rownames(df),]
    
    perm = adonis2(df ~ pH, data = ph_table, permutations = 1000)[1,]
    perm = data.frame(perm[,3:5], cluster = j, Gene = i)
    genes_df = rbind(genes_df, perm)
  }
  rownames(genes_df) = genes_df$Gene
  cluster_table = rbind(cluster_table, genes_df)
}

colnames(cluster_table) = c('R2_perm', 'F_perm', 'p_perm', 'cluster', 'Gene')

alpha_table = read.csv2('singleM_genes/Results/alpha_div_by_cluster_res.csv')
alpha_table$cluster = gsub("\\D", "", alpha_table$cluster)
cluster_table$cluster = gsub("\\D", "", cluster_table$cluster)
cluster_table$Gene[cluster_table$Gene == 'narG'] = "narG_nxrA"
cluster_table$Gene[cluster_table$Gene == 'narH'] = "narH_nxrB"
combined_df = merge(cluster_table, alpha_table)

write.csv2(cluster_table,'singleM_genes/Results/permanova_by_cluster_res.csv')
openxlsx::write.xlsx(combined_df, 'singleM_genes/Results/stats_by_cluster.xlsx')
