library(gsveasyr)
library(changepoint)
library(ggridges)
library(ggplot2)
library(viridis)
library(changepoint.geo)

# env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
# rownames(env) = env$Hoosfield.ID
# env <- env[order(env$pH),]
# load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
#                         changeSampleName = T, refColumn = env$gsveasy_sample,
#                         clustlvls = c(50:100))
# 
# clusters = mget(ls(pattern = 'GSV'))
# rm(list = ls(pattern = 'GSV'))
# genes = unique(names(clusters$GSV_rarefied_gsvtables100))
# 
# load('../HF_nitrogen_paper/heatmaps_all_levels.Rdata')
# 
# change_points_df = data.frame()
# for (j in 1:length(clusters)) {
#   print(j)
#   for (i in genes) {
#     ############################### Change point analysis #########################
#     print(i)
#     df = clusters[[j]][[i]]
#     env_temp = env[,c('Hoosfield.ID', 'pH')]
#     env_temp = env_temp[rownames(env_temp) %in% rownames(df),]
#     df = df[rownames(env_temp), ] |> as.data.frame()
#     env_temp = env_temp[rownames(env_temp) %in% rownames(df),]
# 
#     df_mat = df |> as.matrix() |> vegan::decostand(method = 'normalize', MARGIN = 2)
# 
#     cpt = geomcp(df_mat)
# 
#     dist_cpt = cpt.meanvar(
#       distance(cpt),
#       method = "PELT",
#       penalty = 'CROPS',
#       pen.value = c(5, 500),
#       minseglen = 20
#     )
# 
#     dist_var = cpts.full(dist_cpt)
# 
#     if (all(is.na(dist_var)))
#       next
#     penalties = pen.value.full(dist_cpt)
# 
#     if (nrow(dist_var) == 1) {
#       cp_dist = sum(!is.na(dist_var[1,]))
#     } else {
#       cutoff = cpt.meanvar(
#         penalties,
#         method = "PELT",
#         penalty = 'CROPS',
#         pen.value = c(5, 500)
#       ) |> cpts(ncpt = 1)
#       cp_dist = sum(!is.na(dist_var[cutoff[1],]))
#     }
# 
#     ang_cpt = cpt.meanvar(
#       angle(cpt),
#       method = "PELT",
#       penalty = 'CROPS',
#       pen.value = c(5, 500),
#       minseglen = 20
#     )
# 
#     ang_var = cpts.full(ang_cpt)
# 
#     if (all(is.na(ang_var)))
#       next
#     penalties = pen.value.full(ang_cpt)
#     if (nrow(ang_var) == 1) {
#       cp_angle = sum(!is.na(ang_var[1,]))
#     } else {
#       cutoff = cpt.meanvar(
#         penalties,
#         method = "PELT",
#         penalty = 'CROPS',
#         pen.value = c(5, 500)
#       ) |> cpts(ncpt = 1)
#       cp_angle = sum(!is.na(ang_var[cutoff[1],]))
#     }
# 
#     change_points = c(cpts(dist_cpt, cp_dist)[!is.na(cpts(dist_cpt, cp_dist))],
#     cpts(ang_cpt, cp_angle)[!is.na(cpts(ang_cpt, cp_angle))])
# 
#     if(length(change_points) == 0) next
#     samples = rownames(df)[change_points]
# 
#     flattened_df = cbind(`Hoosfield.ID` = samples,
#                          Cluster = j,
#                          Gene = i)
#     flattened_df = merge(flattened_df, env_temp, all.x = T)
# 
#     change_points_df = rbind(change_points_df, flattened_df)
#   }
# }
# 
# change_points_df$Gene[change_points_df$Gene == 'narG'] = 'narG_nxrA'
# change_points_df$Gene[change_points_df$Gene == 'narH'] = 'narH_nxrB'
# 
# gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
# gene_order = gene_order[gene_order$Gene %in% unique(change_points_df$Gene),]
# change_points_df$Gene = factor(change_points_df$Gene, levels = rev(gene_order$Gene))
# 
# saveRDS(change_points_df, '../HF_nitrogen_paper/change_points.rds')

change_points_df = readRDS('../HF_nitrogen_paper/change_points.rds')

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))
joyplot = ggplot(change_points_df, aes(x = pH, y = Gene, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 3) +
  scale_fill_gradientn(colours = mypal.pH(256)) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlim(3.5,8.0) +
  ggtitle('GSVs composition')

# ggsave('singleM_genes/Figures/change_point.svg', joyplot, width = 5, height = 7)

binplot = ggplot(change_points_df, aes(x = pH, y = Gene, fill = after_stat(x))) +
  geom_density_ridges_gradient(stat = 'binline', bins = 30, scale = 3) +
  scale_fill_gradientn(colours = mypal.pH(256)) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlim(3.5,8.0) +
  ggtitle('GSVs composition')

# ggsave('singleM_genes/Figures/change_point_bin.svg', binplot, width = 5, height = 7)

frequency_table = table(change_points_df[, c(3,4)]) |> as.data.frame()

a = frequency_table[frequency_table$Gene == 'nirK' & frequency_table$Freq != 0,]
