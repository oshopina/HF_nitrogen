library(gsveasyr)

env = read.csv('graftM_genes/Data/mag_env_for_shotgun_samples.csv')
rownames(env) = env$Hoosfield.ID
load_alpha_diversity('singleM_genes/Data/', env.df = env,
                     changeSampleName = T, refColumn = env$gsveasy_sample, 
                     clustlvls = 50:100)

metrics = c('C_sobs', 'C_shannon', 'C_chao1')
alpha_clusters = mget(ls(pattern = 'GSV'))
rm(list = ls(pattern = 'GSV'))

alpha_results = list()
for (v in metrics) {
  alpha_r2 = data.frame()
  for (j in names(alpha_clusters)) {
    anova = auto_aov_fixed(alpha_clusters[[j]][[v]], ~ pH, alpha_clusters[[j]]$C.env)$Results
    signif = anova[!is.na(anova$p_value),]
    gene_names = as.character(signif$Data)
    
    r2s = c()
    for (i in gene_names) {
      r2 = summary(lm(
        alpha_clusters[[j]][[v]][[i]] ~ pH,
        alpha_clusters[[j]]$C.env
      ))$r.squared
      r2s[i] = r2
    }
    
    temp = cbind(signif, R2 = r2s, cluster = j)
    alpha_r2 = rbind(alpha_r2, temp)
  }
  alpha_results[[v]] = alpha_r2
}

sobs_df = alpha_results$C_sobs[,c('Data', 'F_value', 'p_value', 'R2', 'cluster')]
shannon_df = alpha_results$C_shannon[,c('Data', 'F_value', 'p_value', 'R2', 'cluster')]
chao1_df = alpha_results$C_chao1[,c('Data', 'F_value', 'p_value', 'R2', 'cluster')]

colnames(sobs_df) = c('Gene', 'F_value_sobs', 'p_value_sobs', 'R2_sobs', 'cluster')
colnames(shannon_df) = c('Gene', 'F_value_shannon', 'p_value_shannon', 'R2_shannon', 'cluster')
colnames(chao1_df) = c('Gene', 'F_value_chao1', 'p_value_chao1', 'R2_chao1', 'cluster')

combined_df = merge(sobs_df, shannon_df)
combined_df = merge(combined_df, chao1_df)

write.csv2(combined_df, 'singleM_genes/Results/alpha_div_by_cluster_res.csv')

