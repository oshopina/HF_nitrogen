library(gsveasyr)
library(changepoint)
library(ggridges)
library(ggplot2)
library(viridis)

env = read.csv('graftM_genes/Data/mag_env_for_shotgun_samples.csv')
rownames(env) = env$Hoosfield.ID
best_cluster = openxlsx::read.xlsx('singleM_genes/Results/best_cluster.xlsx')
load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
                        changeSampleName = T, refColumn = env$gsveasy_sample, 
                        clustlvls = unique(best_cluster$cluster))

################# Choose the best clustering level for each gene ###############
genes = list(anfH = GSV_rarefied_gsvtables73$anfH,
             asnB = GSV_rarefied_gsvtables67$asnB,
             gdhA = GSV_rarefied_gsvtables86$gdhA,
             glnA = GSV_rarefied_gsvtables51$glnA,
             glsA = GSV_rarefied_gsvtables99$glsA,
             gltB = GSV_rarefied_gsvtables96$gltB,
             gltD = GSV_rarefied_gsvtables70$gltD,
             gltS = GSV_rarefied_gsvtables86$gltS,
             gudB = GSV_rarefied_gsvtables72$gudB,
             hzsA = GSV_rarefied_gsvtables79$hzsA,
             hzsB = GSV_rarefied_gsvtables71$hzsB,
             hzsC = GSV_rarefied_gsvtables74$hzsC,
             napA = GSV_rarefied_gsvtables66$napA,
             narB = GSV_rarefied_gsvtables81$narB,
             narC = GSV_rarefied_gsvtables68$narC,
             narG_nxrA = GSV_rarefied_gsvtables78$narG,
             narH_nxrB = GSV_rarefied_gsvtables67$narH,
             narI = GSV_rarefied_gsvtables93$narI,
             narJ = GSV_rarefied_gsvtables92$narJ,
             nasA = GSV_rarefied_gsvtables75$nasA,
             nasB = GSV_rarefied_gsvtables66$nasB,
             nifH = GSV_rarefied_gsvtables78$nifH,
             nirA = GSV_rarefied_gsvtables96$nirA,
             nirB = GSV_rarefied_gsvtables76$nirB,
             nirD = GSV_rarefied_gsvtables97$nirD,
             nirK = GSV_rarefied_gsvtables68$nirK,
             nod = GSV_rarefied_gsvtables80$nod,
             norB = GSV_rarefied_gsvtables77$norB,
             norC = GSV_rarefied_gsvtables70$norC,
             nosZ = GSV_rarefied_gsvtables88$nosZ,
             NR = GSV_rarefied_gsvtables95$NR,
             nrfA = GSV_rarefied_gsvtables57$nrfA,
             nrfB = GSV_rarefied_gsvtables95$nrfB,
             nrfC = GSV_rarefied_gsvtables64$nrfC,
             nrfD = GSV_rarefied_gsvtables80$nrfD,
             nrfH = GSV_rarefied_gsvtables86$nrfH,
             ureA = GSV_rarefied_gsvtables90$ureA,
             ureB = GSV_rarefied_gsvtables86$ureB,
             ureC = GSV_rarefied_gsvtables89$ureC,
             ureD = GSV_rarefied_gsvtables60$ureD,
             ureE = GSV_rarefied_gsvtables74$ureE,
             ureF = GSV_rarefied_gsvtables83$ureF,
             ureG = GSV_rarefied_gsvtables89$ureG,
             ureJ = GSV_rarefied_gsvtables69$ureJ,
             vnfH = GSV_rarefied_gsvtables82$vnfH)
rm(list = ls(pattern = 'GSV'))
###############################################################################
env_sorted = env[order(env$pH),]
env_sorted = env_sorted[rownames(env_sorted) %in% rownames(genes$anfH),]
env_sorted$change_point = c(1:119)
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
  cpt = cpt.meanvar(t(df), method="PELT",
                    penalty = "Manual", pen.value = 30, minseglen = 30)
  names(cpt) = colnames(as.data.frame(df))
  points = lapply(cpt, cpts)
  flattened_vector <- unlist(points)
  flattened_vector = sort(flattened_vector)
  flattened_df = data.frame(GSV = names(flattened_vector), change_point = flattened_vector)
  flattened_df = rbind(flattened_df, flattened_df, flattened_df, flattened_df, flattened_df, flattened_df, flattened_df) ## Add scale to tipping points
  flattened_df = merge(flattened_df, pH_list, all.y = T)
  flattened_df$Gene = i
  result = rbind(result, flattened_df)
}

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))

gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% names(genes),]
result$Gene = factor(result$Gene, levels = rev(gene_order$Gene))

ggplot(result, aes(x = pH, y = Gene, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_gradientn(colours = mypal.pH(256)) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlim(3.5,8.0)




