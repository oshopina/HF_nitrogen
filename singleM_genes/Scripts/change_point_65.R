library(gsveasyr)
library(changepoint)
library(ggridges)
library(ggplot2)
library(viridis)

env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
                        changeSampleName = T, refColumn = env$gsveasy_sample, 
                        clustlvls = c(65, 85))

################# Choose the best clustering level for each gene ###############
genes65 = GSV_rarefied_gsvtables65
genes85 = GSV_rarefied_gsvtables85
rm(list = ls(pattern = 'GSV'))
###############################################################################
env_sorted = env[order(env$pH),]
env_sorted = env_sorted[rownames(env_sorted) %in% rownames(genes$anfH),]
env_sorted$change_point = c(1:117)
sample_list = data.frame(ID = env_sorted$Hoosfield.ID)
pH_list = env_sorted[,c('change_point', 'pH')]

## Change point analysis
result = data.frame()

for (i in names(genes65)) {
  print(i)
  df = genes65[[i]] %>% as.data.frame()
  df$ID = rownames(df)
  df = merge(df, sample_list, all.y = T)
  df[is.na(df) == T] = 0
  rownames(df) = df$ID
  df = df[sample_list$ID,]
  df = df[,-1]
  cpt = cpt.meanvar(t(df), method="PELT", minseglen = 20)
  names(cpt) = colnames(as.data.frame(df))
  points = lapply(cpt, cpts)
  flattened_vector <- unlist(points)
  flattened_vector = sort(flattened_vector)
  flattened_df = data.frame(GSV = names(flattened_vector), change_point = flattened_vector)
  flattened_df = rbind(flattened_df, flattened_df, flattened_df, flattened_df, flattened_df, flattened_df, flattened_df) ## Add scale to tipping points
  flattened_df = merge(flattened_df, pH_list, all.y = T)
  flattened_df$Gene = i
  if (length(flattened_vector) > 0) {
    result = rbind(result, flattened_df) 
  }
}

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))

gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% unique(result$Gene),]
result$Gene = factor(result$Gene, levels = rev(gene_order$Gene))

joyplot = ggplot(result, aes(x = pH, y = Gene, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_gradientn(colours = mypal.pH(256)) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlim(3.5,8.0)

########################### Shannon diversity bargraph #########################

load_alpha_diversity('singleM_genes/Data/', env.df = env,
                     changeSampleName = T, refColumn = env$gsveasy_sample, 
                     clustlvls = c(65, 85))


######################### Combine best clustering levels #######################
shannon65 = data.frame(GSV_alpha_div_cluster65$C_shannon,
                     ID = GSV_alpha_div_cluster65$C.env$Hoosfield.ID)
shannon85 = data.frame(GSV_alpha_div_cluster85$C_shannon,
                       ID = GSV_alpha_div_cluster85$C.env$Hoosfield.ID)

rm(list = ls(pattern = 'GSV'))

shannon65 = shannon65[shannon65$ID %in% env$Hoosfield.ID,]
rownames(shannon65) = shannon65$ID
shannon65 = shannon65[,-ncol(shannon65)]

shannon85 = shannon85[shannon85$ID %in% env$Hoosfield.ID,]
rownames(shannon85) = shannon85$ID
shannon85 = shannon85[,-ncol(shannon85)]
###############################################################################

df = shannon65
df$ID = rownames(df)
df = merge(df, sample_list, all.y = T)
df[is.na(df) == T] = 0
rownames(df) = df$ID
df = df[sample_list$ID,]
df = df[,-1]
cpt = cpt.meanvar(t(df), method="PELT", minseglen = 20)
names(cpt) = colnames(as.data.frame(df))
points = lapply(cpt, cpts)
flattened_vector <- unlist(points)
flattened_vector = sort(flattened_vector)
flattened_df = data.frame(GSV = names(flattened_vector), change_point = flattened_vector)
flattened_df = rbind(flattened_df, flattened_df, flattened_df, flattened_df, flattened_df, flattened_df, flattened_df) ## Add scale to tipping points
flattened_df = merge(flattened_df, pH_list, all.y = T)
flattened_df$Alpha = 'Shannon'

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))

alpha_plot = ggplot(flattened_df) +
  geom_histogram(aes(x = change_point), fill = mypal.pH(117), binwidth = 1) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank()) +
  xlab('') +
  ylab('Shannon diversity')

library(patchwork)

alpha_plot + joyplot + plot_layout(ncol = 1, nrow = 2, heights = c(1, 9))


#### 85 ####

## Change point analysis
result = data.frame()

for (i in names(genes85)) {
  print(i)
  df = genes85[[i]] %>% as.data.frame()
  df$ID = rownames(df)
  df = merge(df, sample_list, all.y = T)
  df[is.na(df) == T] = 0
  rownames(df) = df$ID
  df = df[sample_list$ID,]
  df = df[,-1]
  cpt = cpt.meanvar(t(df), method="PELT", minseglen = 20)
  names(cpt) = colnames(as.data.frame(df))
  points = lapply(cpt, cpts)
  flattened_vector <- unlist(points)
  flattened_vector = sort(flattened_vector)
  flattened_df = data.frame(GSV = names(flattened_vector), change_point = flattened_vector)
  flattened_df = rbind(flattened_df, flattened_df, flattened_df, flattened_df, flattened_df, flattened_df, flattened_df) ## Add scale to tipping points
  flattened_df = merge(flattened_df, pH_list, all.y = T)
  flattened_df$Gene = i
  if (length(flattened_vector) > 0) {
    result = rbind(result, flattened_df) 
  }
}

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))

gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% unique(result$Gene),]
result$Gene = factor(result$Gene, levels = rev(gene_order$Gene))

joyplot = ggplot(result, aes(x = pH, y = Gene, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_gradientn(colours = mypal.pH(256)) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlim(3.5,8.0)


df = shannon85
df$ID = rownames(df)
df = merge(df, sample_list, all.y = T)
df[is.na(df) == T] = 0
rownames(df) = df$ID
df = df[sample_list$ID,]
df = df[,-1]
cpt = cpt.meanvar(t(df), method="PELT", minseglen = 20)
names(cpt) = colnames(as.data.frame(df))
points = lapply(cpt, cpts)
flattened_vector <- unlist(points)
flattened_vector = sort(flattened_vector)
flattened_df = data.frame(GSV = names(flattened_vector), change_point = flattened_vector)
flattened_df = rbind(flattened_df, flattened_df, flattened_df, flattened_df, flattened_df, flattened_df, flattened_df) ## Add scale to tipping points
flattened_df = merge(flattened_df, pH_list, all.y = T)
flattened_df$Alpha = 'Shannon'

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))

alpha_plot = ggplot(flattened_df) +
  geom_histogram(aes(x = change_point), fill = mypal.pH(117), binwidth = 1) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank()) +
  xlab('') +
  ylab('Shannon diversity')

library(patchwork)

alpha_plot + joyplot + plot_layout(ncol = 1, nrow = 2, heights = c(1, 9))
