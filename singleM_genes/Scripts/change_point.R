library(gsveasyr)
library(changepoint)
library(ggridges)
library(ggplot2)
library(viridis)

env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
best_cluster = openxlsx::read.xlsx('singleM_genes/Results/best_cluster.xlsx')
load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
                        changeSampleName = T, refColumn = env$gsveasy_sample, 
                        clustlvls = unique(best_cluster$cluster))

################# Choose the best clustering level for each gene ###############
genes = list(anfH = GSV_rarefied_gsvtables70$anfH,
                     asnB = GSV_rarefied_gsvtables59$asnB,
                     gdhA = GSV_rarefied_gsvtables67$gdhA,
                     glnA = GSV_rarefied_gsvtables51$glnA,
                     glsA = GSV_rarefied_gsvtables99$glsA,
                     gltB = GSV_rarefied_gsvtables80$gltB,
                     gltD = GSV_rarefied_gsvtables68$gltD,
                     gltS = GSV_rarefied_gsvtables81$gltS,
                     gudB = GSV_rarefied_gsvtables62$gudB,
                     hzsA = GSV_rarefied_gsvtables60$hzsA,
                     hzsB = GSV_rarefied_gsvtables62$hzsB,
                     hzsC = GSV_rarefied_gsvtables74$hzsC,
                     napA = GSV_rarefied_gsvtables53$napA,
                     narB = GSV_rarefied_gsvtables58$narB,
                     narC = GSV_rarefied_gsvtables68$narC,
                     narG_nxrA = GSV_rarefied_gsvtables67$narG,
                     narH_nxrB = GSV_rarefied_gsvtables59$narH,
                     narI = GSV_rarefied_gsvtables68$narI,
                     narJ = GSV_rarefied_gsvtables80$narJ,
                     nasA = GSV_rarefied_gsvtables51$nasA,
                     nasB = GSV_rarefied_gsvtables60$nasB,
                     nifH = GSV_rarefied_gsvtables76$nifH,
                     nirA = GSV_rarefied_gsvtables87$nirA,
                     nirB = GSV_rarefied_gsvtables72$nirB,
                     nirD = GSV_rarefied_gsvtables64$nirD,
                     nirK = GSV_rarefied_gsvtables60$nirK,
                     nod = GSV_rarefied_gsvtables74$nod,
                     norB = GSV_rarefied_gsvtables56$norB,
                     norC = GSV_rarefied_gsvtables65$norC,
                     nosZ = GSV_rarefied_gsvtables51$nosZ,
                     NR = GSV_rarefied_gsvtables53$NR,
                     nrfA = GSV_rarefied_gsvtables57$nrfA,
                     nrfB = GSV_rarefied_gsvtables52$nrfB,
                     nrfC = GSV_rarefied_gsvtables52$nrfC,
                     nrfD = GSV_rarefied_gsvtables51$nrfD,
                     nrfH = GSV_rarefied_gsvtables67$nrfH,
                     ureA = GSV_rarefied_gsvtables57$ureA,
                     ureB = GSV_rarefied_gsvtables87$ureB,
                     ureC = GSV_rarefied_gsvtables89$ureC,
                     ureD = GSV_rarefied_gsvtables60$ureD,
                     ureE = GSV_rarefied_gsvtables51$ureE,
                     ureF = GSV_rarefied_gsvtables75$ureF,
                     ureG = GSV_rarefied_gsvtables73$ureG,
                     ureJ = GSV_rarefied_gsvtables84$ureJ,
                     vnfH = GSV_rarefied_gsvtables51$vnfH)

rm(list = ls(pattern = 'GSV'))
###############################################################################
env_sorted = env[order(env$pH),]
env_sorted = env_sorted[rownames(env_sorted) %in% rownames(genes$anfH),]
env_sorted$change_point = c(1:117)
sample_list = data.frame(ID = env_sorted$Hoosfield.ID)
pH_list = env_sorted[,c('change_point', 'pH')]

## Change point analysis
result = data.frame()

for (i in names(genes)) {
  print(i)
  df = genes[[i]] %>% as.data.frame()
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
                     clustlvls = unique(best_cluster$cluster))


######################### Combine best clustering levels #######################
shannon = data.frame(anfH = GSV_alpha_div_cluster70$C_shannon$anfH,
                     asnB = GSV_alpha_div_cluster59$C_shannon$asnB,
                     gdhA = GSV_alpha_div_cluster67$C_shannon$gdhA,
                     glnA = GSV_alpha_div_cluster51$C_shannon$glnA,
                     glsA = GSV_alpha_div_cluster99$C_shannon$glsA,
                     gltB = GSV_alpha_div_cluster80$C_shannon$gltB,
                     gltD = GSV_alpha_div_cluster68$C_shannon$gltD,
                     gltS = GSV_alpha_div_cluster81$C_shannon$gltS,
                     gudB = GSV_alpha_div_cluster62$C_shannon$gudB,
                     hzsA = GSV_alpha_div_cluster60$C_shannon$hzsA,
                     hzsB = GSV_alpha_div_cluster62$C_shannon$hzsB,
                     hzsC = GSV_alpha_div_cluster74$C_shannon$hzsC,
                     napA = GSV_alpha_div_cluster53$C_shannon$napA,
                     narB = GSV_alpha_div_cluster58$C_shannon$narB,
                     narC = GSV_alpha_div_cluster68$C_shannon$narC,
                     narG_nxrA = GSV_alpha_div_cluster67$C_shannon$narG,
                     narH_nxrB = GSV_alpha_div_cluster59$C_shannon$narH,
                     narI = GSV_alpha_div_cluster68$C_shannon$narI,
                     narJ = GSV_alpha_div_cluster80$C_shannon$narJ,
                     nasA = GSV_alpha_div_cluster51$C_shannon$nasA,
                     nasB = GSV_alpha_div_cluster60$C_shannon$nasB,
                     nifH = GSV_alpha_div_cluster76$C_shannon$nifH,
                     nirA = GSV_alpha_div_cluster87$C_shannon$nirA,
                     nirB = GSV_alpha_div_cluster72$C_shannon$nirB,
                     nirD = GSV_alpha_div_cluster64$C_shannon$nirD,
                     nirK = GSV_alpha_div_cluster60$C_shannon$nirK,
                     nod = GSV_alpha_div_cluster74$C_shannon$nod,
                     norB = GSV_alpha_div_cluster56$C_shannon$norB,
                     norC = GSV_alpha_div_cluster65$C_shannon$norC,
                     nosZ = GSV_alpha_div_cluster51$C_shannon$nosZ,
                     NR = GSV_alpha_div_cluster53$C_shannon$NR,
                     nrfA = GSV_alpha_div_cluster57$C_shannon$nrfA,
                     nrfB = GSV_alpha_div_cluster52$C_shannon$nrfB,
                     nrfC = GSV_alpha_div_cluster52$C_shannon$nrfC,
                     nrfD = GSV_alpha_div_cluster51$C_shannon$nrfD,
                     nrfH = GSV_alpha_div_cluster67$C_shannon$nrfH,
                     ureA = GSV_alpha_div_cluster57$C_shannon$ureA,
                     ureB = GSV_alpha_div_cluster87$C_shannon$ureB,
                     ureC = GSV_alpha_div_cluster89$C_shannon$ureC,
                     ureD = GSV_alpha_div_cluster60$C_shannon$ureD,
                     ureE = GSV_alpha_div_cluster51$C_shannon$ureE,
                     ureF = GSV_alpha_div_cluster75$C_shannon$ureF,
                     ureG = GSV_alpha_div_cluster73$C_shannon$ureG,
                     ureJ = GSV_alpha_div_cluster84$C_shannon$ureJ,
                     vnfH = GSV_alpha_div_cluster51$C_shannon$vnfH,
                     ID = GSV_alpha_div_cluster51$C.env$Hoosfield.ID)

rm(list = ls(pattern = 'GSV'))

shannon = shannon[shannon$ID %in% env$Hoosfield.ID,]
rownames(shannon) = shannon$ID
shannon = shannon[,-ncol(shannon)]

###############################################################################

df = shannon
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
