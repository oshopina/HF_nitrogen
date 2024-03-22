library(gsveasyr)
library(changepoint)
library(ggridges)
library(ggplot2)
library(viridis)
library(changepoint.geo)

env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
env <- env[order(env$pH),]
load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
                        changeSampleName = T, refColumn = env$gsveasy_sample, 
                        clustlvls = c(50:100))

clusters = mget(ls(pattern = 'GSV'))
rm(list = ls(pattern = 'GSV'))
genes = unique(names(clusters$GSV_rarefied_gsvtables100))

change_points_df = data.frame()
for (j in names(clusters)) {
  print(j)
  for (i in genes) {
    ############################### Change point analysis #########################
    print(i)
    df = clusters[[j]][[i]]
    env_temp = env[,c('Hoosfield.ID', 'pH')]
    env_temp = env_temp[rownames(env_temp) %in% rownames(df),]
    df = df[rownames(env_temp), ]
    env_temp = env_temp[rownames(env) %in% rownames(df),]
    cpt = geomcp(df)
    
    dist_cpt = cpt.meanvar(
      distance(cpt),
      method = "PELT",
      penalty = 'CROPS',
      pen.value = c(5, 500)
    )
    
    dist_var = cpts.full(dist_cpt)
    ######################### WRITE CHOSEN NUMBER OF CHANGEPOINTS FOR ANGLE AND DISTANCE
    if(all(is.na(dist_var))) next
    cp_dist = sum(!is.na(dist_var[nrow(dist_var) - 2,]))
    
    ang_cpt = cpt.meanvar(
      angle(cpt),
      method = "PELT",
      penalty = 'CROPS',
      pen.value = c(5, 500)
    )
    
    ang_var = cpts.full(ang_cpt)
    ######################### WRITE CHOSEN NUMBER OF CHANGEPOINTS FOR ANGLE AND DISTANCE
    cp_angle = sum(!is.na(ang_var[nrow(ang_var) - 2,]))
    
    change_points = c(cpts(dist_cpt, cp_dist)[!is.na(cpts(dist_cpt, cp_dist))], cpts(ang_cpt, cp_angle)[!is.na(cpts(ang_cpt, cp_angle))])
    samples = rownames(df[change_points,])
    flattened_df = cbind(`Hoosfield.ID` = samples,
                         Cluster = j,
                         Gene = i)
    flattened_df = replicate(100, flattened_df, simplify = FALSE)
    flattened_df = do.call(rbind, flattened_df)
    flattened_df = merge(flattened_df, env_temp, all.y = T)
    flattened_df$Gene = i
    flattened_df$Cluster = j
    
    change_points_df = rbind(change_points_df, flattened_df)
  }
}

change_points_df$Gene[change_points_df$Gene == 'narG'] = 'narG_nxrA'
change_points_df$Gene[change_points_df$Gene == 'narH'] = 'narH_nxrB'

gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% unique(change_points_df$Gene),]
change_points_df$Gene = factor(change_points_df$Gene, levels = rev(gene_order$Gene))
change_points_df = change_points_df[,-c(1:2)]

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))
joyplot = ggplot(change_points_df, aes(x = pH, y = Gene, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 3) +
  scale_fill_gradientn(colours = mypal.pH(256)) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlim(3.5,8.0) +
  ggtitle('GSVs composition')


