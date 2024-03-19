library(gsveasyr)

env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
                     changeSampleName = T, refColumn = env$gsveasy_sample, 
                     clustlvls = c(50:100))
clusters = mget(ls(pattern = 'GSV'))
rm(list = ls(pattern = 'GSV'))
genes = unique(names(clusters$GSV_rarefied_gsvtables100))

hellinger_diversity = function(otu_table) {
  dist_matrix = dist(vegan::decostand(otu_table, method = 'hellinger'),
                     method = 'euc') %>% as.matrix() 
  return(dist_matrix)
}

distance_tables = list()
for (j in names(clusters)) {
  print(paste('Started cluster:', j))
  for (i in genes) {
    print(paste('Started gene:', i))
    df = clusters[[j]][[i]]
    distance = hellinger_diversity(df)
    distance_tables[[j]][[i]] = distance
  }
}

cluster_table = data.frame()
for (j in names(distance_tables)) {
  print(paste('Started cluster:', j))
  genes_df = data.frame()
  for (i in genes) {
    print(paste('Started gene:', i))
    df = distance_tables[[j]][[i]]
    df = df[rownames(df) %in% rownames(env),]
    df = df[,colnames(df) %in% rownames(env)]
    ph_table = env[rownames(df),]
    
      perm = adonis2(df ~ pH, data = ph_table, permutations = 1000)[1,]
      perm = data.frame(perm[,3:5], cluster = j, Gene = i)
      genes_df = rbind(genes_df, perm)
  }
  rownames(genes_df) = genes_df$Gene
  cluster_table = rbind(cluster_table, genes_df)
}

colnames(cluster_table) = c('R2_perm', 'F_perm', 'p_perm', 'cluster', 'Gene')

write.csv2(cluster_table,'singleM_genes/Results/permanova_by_hellinger_res.csv')

