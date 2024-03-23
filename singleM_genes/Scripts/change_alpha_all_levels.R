library(gsveasyr)
library(changepoint)
library(ggridges)
library(ggplot2)
library(viridis)

env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
env <- env[order(env$pH),]
load_alpha_diversity('singleM_genes/Data/', env.df = env,
                     changeSampleName = T, refColumn = env$gsveasy_sample,
                     clustlvls = 50:100)

clusters = mget(ls(pattern = 'GSV'))
rm(list = ls(pattern = 'GSV'))
genes = unique(names(clusters$GSV_alpha_div_cluster100$C_shannon))


change_points_df = data.frame()
for (j in names(clusters)) {
  print(j)
  for (i in genes) {
    ############################### Change point analysis #########################
    print(i)
    df = clusters[[j]][['C_shannon']][,i]

    if(all(is.na(df))) next
    names(df) = rownames(clusters[[j]][['C_shannon']])
    df = df[!is.na(df)]
    env_temp = env[,c('Hoosfield.ID', 'pH')]
    env_temp = env_temp[rownames(env_temp) %in% names(df),]
    df = df[rownames(env_temp)]


    dist_cpt = cpt.meanvar(
      df,
      method = "PELT",
      penalty = 'CROPS',
      pen.value = c(5, 500),
      minseglen = 20
    )

    dist_var = cpts.full(dist_cpt)
    ######################### WRITE CHOSEN NUMBER OF CHANGEPOINTS FOR ANGLE AND DISTANCE
    if(all(is.na(dist_var))) next
    penalties = pen.value.full(dist_cpt) |> cumsum() |> diff()
    cutoff = which(penalties < 50)
    if(length(cutoff) == 0 || nrow(dist_var) == 1) {
      cp_dist = sum(!is.na(dist_var[1,]))
    } else cp_dist = sum(!is.na(dist_var[max(cutoff) + 1,]))

    change_points = c(cpts(dist_cpt, cp_dist)[!is.na(cpts(dist_cpt, cp_dist))])
    if(length(change_points) == 0) next
    samples = names(df)[change_points]

    flattened_df = cbind(`Hoosfield.ID` = samples,
                         Cluster = j,
                         Gene = i)
    flattened_df = merge(flattened_df, env_temp, all.x = T)
    
    change_points_df = rbind(change_points_df, flattened_df)
  }
}


change_points_df$Gene[change_points_df$Gene == 'narG'] = 'narG_nxrA'
change_points_df$Gene[change_points_df$Gene == 'narH'] = 'narH_nxrB'

gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% unique(change_points_df$Gene),]
change_points_df$Gene = factor(change_points_df$Gene, levels = rev(gene_order$Gene))


mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))
joyplot = ggplot(change_points_df, aes(x = pH, y = Gene, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 3) +
  scale_fill_gradientn(colours = mypal.pH(256)) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlim(3.5,8.0) +
  ggtitle('Shannon diversity')
