library(gsveasyr)
library(ggridges)
library(ggplot2)
library(viridis)
library(changepoint)

env = read.csv('graftM_genes/Data/mag_env_for_shotgun_samples.csv')
rownames(env) = env$Hoosfield.ID
best_cluster = openxlsx::read.xlsx('singleM_genes/Results/best_cluster.xlsx')
load_alpha_diversity('singleM_genes/Data/', env.df = env,
                     changeSampleName = T, refColumn = env$gsveasy_sample, 
                     clustlvls = unique(best_cluster$cluster))

######################### Combine best clustering levels #######################
sobs = data.frame(anfH = GSV_alpha_div_cluster73$C_sobs$anfH,
                  asnB = GSV_alpha_div_cluster67$C_sobs$asnB,
                  gdhA = GSV_alpha_div_cluster86$C_sobs$gdhA,
                  glnA = GSV_alpha_div_cluster51$C_sobs$glnA,
                  glsA = GSV_alpha_div_cluster99$C_sobs$glsA,
                  gltB = GSV_alpha_div_cluster96$C_sobs$gltB,
                  gltD = GSV_alpha_div_cluster70$C_sobs$gltD,
                  gltS = GSV_alpha_div_cluster86$C_sobs$gltS,
                  gudB = GSV_alpha_div_cluster72$C_sobs$gudB,
                  hzsA = GSV_alpha_div_cluster79$C_sobs$hzsA,
                  hzsB = GSV_alpha_div_cluster71$C_sobs$hzsB,
                  hzsC = GSV_alpha_div_cluster74$C_sobs$hzsC,
                  napA = GSV_alpha_div_cluster66$C_sobs$napA,
                  narB = GSV_alpha_div_cluster81$C_sobs$narB,
                  narC = GSV_alpha_div_cluster68$C_sobs$narC,
                  narG_nxrA = GSV_alpha_div_cluster78$C_sobs$narG,
                  narH_nxrB = GSV_alpha_div_cluster67$C_sobs$narH,
                  narI = GSV_alpha_div_cluster93$C_sobs$narI,
                  narJ = GSV_alpha_div_cluster92$C_sobs$narJ,
                  nasA = GSV_alpha_div_cluster75$C_sobs$nasA,
                  nasB = GSV_alpha_div_cluster66$C_sobs$nasB,
                  nifH = GSV_alpha_div_cluster78$C_sobs$nifH,
                  nirA = GSV_alpha_div_cluster96$C_sobs$nirA,
                  nirB = GSV_alpha_div_cluster76$C_sobs$nirB,
                  nirD = GSV_alpha_div_cluster97$C_sobs$nirD,
                  nirK = GSV_alpha_div_cluster68$C_sobs$nirK,
                  nod = GSV_alpha_div_cluster80$C_sobs$nod,
                  norB = GSV_alpha_div_cluster77$C_sobs$norB,
                  norC = GSV_alpha_div_cluster70$C_sobs$norC,
                  nosZ = GSV_alpha_div_cluster88$C_sobs$nosZ,
                  NR = GSV_alpha_div_cluster95$C_sobs$NR,
                  nrfA = GSV_alpha_div_cluster57$C_sobs$nrfA,
                  nrfB = GSV_alpha_div_cluster95$C_sobs$nrfB,
                  nrfC = GSV_alpha_div_cluster64$C_sobs$nrfC,
                  nrfD = GSV_alpha_div_cluster80$C_sobs$nrfD,
                  nrfH = GSV_alpha_div_cluster86$C_sobs$nrfH,
                  ureA = GSV_alpha_div_cluster90$C_sobs$ureA,
                  ureB = GSV_alpha_div_cluster86$C_sobs$ureB,
                  ureC = GSV_alpha_div_cluster89$C_sobs$ureC,
                  ureD = GSV_alpha_div_cluster60$C_sobs$ureD,
                  ureE = GSV_alpha_div_cluster74$C_sobs$ureE,
                  ureF = GSV_alpha_div_cluster83$C_sobs$ureF,
                  ureG = GSV_alpha_div_cluster89$C_sobs$ureG,
                  ureJ = GSV_alpha_div_cluster69$C_sobs$ureJ,
                  vnfH = GSV_alpha_div_cluster82$C_sobs$vnfH)
rownames(sobs) = GSV_alpha_div_cluster51$C.env$Hoosfield.ID

shannon = data.frame(anfH = GSV_alpha_div_cluster73$C_shannon$anfH,
                     asnB = GSV_alpha_div_cluster67$C_shannon$asnB,
                     gdhA = GSV_alpha_div_cluster86$C_shannon$gdhA,
                     glnA = GSV_alpha_div_cluster51$C_shannon$glnA,
                     glsA = GSV_alpha_div_cluster99$C_shannon$glsA,
                     gltB = GSV_alpha_div_cluster96$C_shannon$gltB,
                     gltD = GSV_alpha_div_cluster70$C_shannon$gltD,
                     gltS = GSV_alpha_div_cluster86$C_shannon$gltS,
                     gudB = GSV_alpha_div_cluster72$C_shannon$gudB,
                     hzsA = GSV_alpha_div_cluster79$C_shannon$hzsA,
                     hzsB = GSV_alpha_div_cluster71$C_shannon$hzsB,
                     hzsC = GSV_alpha_div_cluster74$C_shannon$hzsC,
                     napA = GSV_alpha_div_cluster66$C_shannon$napA,
                     narB = GSV_alpha_div_cluster81$C_shannon$narB,
                     narC = GSV_alpha_div_cluster68$C_shannon$narC,
                     narG_nxrA = GSV_alpha_div_cluster78$C_shannon$narG,
                     narH_nxrB = GSV_alpha_div_cluster67$C_shannon$narH,
                     narI = GSV_alpha_div_cluster93$C_shannon$narI,
                     narJ = GSV_alpha_div_cluster92$C_shannon$narJ,
                     nasA = GSV_alpha_div_cluster75$C_shannon$nasA,
                     nasB = GSV_alpha_div_cluster66$C_shannon$nasB,
                     nifH = GSV_alpha_div_cluster78$C_shannon$nifH,
                     nirA = GSV_alpha_div_cluster96$C_shannon$nirA,
                     nirB = GSV_alpha_div_cluster76$C_shannon$nirB,
                     nirD = GSV_alpha_div_cluster97$C_shannon$nirD,
                     nirK = GSV_alpha_div_cluster68$C_shannon$nirK,
                     nod = GSV_alpha_div_cluster80$C_shannon$nod,
                     norB = GSV_alpha_div_cluster77$C_shannon$norB,
                     norC = GSV_alpha_div_cluster70$C_shannon$norC,
                     nosZ = GSV_alpha_div_cluster88$C_shannon$nosZ,
                     NR = GSV_alpha_div_cluster95$C_shannon$NR,
                     nrfA = GSV_alpha_div_cluster57$C_shannon$nrfA,
                     nrfB = GSV_alpha_div_cluster95$C_shannon$nrfB,
                     nrfC = GSV_alpha_div_cluster64$C_shannon$nrfC,
                     nrfD = GSV_alpha_div_cluster80$C_shannon$nrfD,
                     nrfH = GSV_alpha_div_cluster86$C_shannon$nrfH,
                     ureA = GSV_alpha_div_cluster90$C_shannon$ureA,
                     ureB = GSV_alpha_div_cluster86$C_shannon$ureB,
                     ureC = GSV_alpha_div_cluster89$C_shannon$ureC,
                     ureD = GSV_alpha_div_cluster60$C_shannon$ureD,
                     ureE = GSV_alpha_div_cluster74$C_shannon$ureE,
                     ureF = GSV_alpha_div_cluster83$C_shannon$ureF,
                     ureG = GSV_alpha_div_cluster89$C_shannon$ureG,
                     ureJ = GSV_alpha_div_cluster69$C_shannon$ureJ,
                     vnfH = GSV_alpha_div_cluster82$C_shannon$vnfH)
rownames(shannon) = GSV_alpha_div_cluster51$C.env$Hoosfield.ID

chao1 = data.frame(anfH = GSV_alpha_div_cluster73$C_chao1$anfH,
                   asnB = GSV_alpha_div_cluster67$C_chao1$asnB,
                   gdhA = GSV_alpha_div_cluster86$C_chao1$gdhA,
                   glnA = GSV_alpha_div_cluster51$C_chao1$glnA,
                   glsA = GSV_alpha_div_cluster99$C_chao1$glsA,
                   gltB = GSV_alpha_div_cluster96$C_chao1$gltB,
                   gltD = GSV_alpha_div_cluster70$C_chao1$gltD,
                   gltS = GSV_alpha_div_cluster86$C_chao1$gltS,
                   gudB = GSV_alpha_div_cluster72$C_chao1$gudB,
                   hzsA = GSV_alpha_div_cluster79$C_chao1$hzsA,
                   hzsB = GSV_alpha_div_cluster71$C_chao1$hzsB,
                   hzsC = GSV_alpha_div_cluster74$C_chao1$hzsC,
                   napA = GSV_alpha_div_cluster66$C_chao1$napA,
                   narB = GSV_alpha_div_cluster81$C_chao1$narB,
                   narC = GSV_alpha_div_cluster68$C_chao1$narC,
                   narG_nxrA = GSV_alpha_div_cluster78$C_chao1$narG,
                   narH_nxrB = GSV_alpha_div_cluster67$C_chao1$narH,
                   narI = GSV_alpha_div_cluster93$C_chao1$narI,
                   narJ = GSV_alpha_div_cluster92$C_chao1$narJ,
                   nasA = GSV_alpha_div_cluster75$C_chao1$nasA,
                   nasB = GSV_alpha_div_cluster66$C_chao1$nasB,
                   nifH = GSV_alpha_div_cluster78$C_chao1$nifH,
                   nirA = GSV_alpha_div_cluster96$C_chao1$nirA,
                   nirB = GSV_alpha_div_cluster76$C_chao1$nirB,
                   nirD = GSV_alpha_div_cluster97$C_chao1$nirD,
                   nirK = GSV_alpha_div_cluster68$C_chao1$nirK,
                   nod = GSV_alpha_div_cluster80$C_chao1$nod,
                   norB = GSV_alpha_div_cluster77$C_chao1$norB,
                   norC = GSV_alpha_div_cluster70$C_chao1$norC,
                   nosZ = GSV_alpha_div_cluster88$C_chao1$nosZ,
                   NR = GSV_alpha_div_cluster95$C_chao1$NR,
                   nrfA = GSV_alpha_div_cluster57$C_chao1$nrfA,
                   nrfB = GSV_alpha_div_cluster95$C_chao1$nrfB,
                   nrfC = GSV_alpha_div_cluster64$C_chao1$nrfC,
                   nrfD = GSV_alpha_div_cluster80$C_chao1$nrfD,
                   nrfH = GSV_alpha_div_cluster86$C_chao1$nrfH,
                   ureA = GSV_alpha_div_cluster90$C_chao1$ureA,
                   ureB = GSV_alpha_div_cluster86$C_chao1$ureB,
                   ureC = GSV_alpha_div_cluster89$C_chao1$ureC,
                   ureD = GSV_alpha_div_cluster60$C_chao1$ureD,
                   ureE = GSV_alpha_div_cluster74$C_chao1$ureE,
                   ureF = GSV_alpha_div_cluster83$C_chao1$ureF,
                   ureG = GSV_alpha_div_cluster89$C_chao1$ureG,
                   ureJ = GSV_alpha_div_cluster69$C_chao1$ureJ,
                   vnfH = GSV_alpha_div_cluster82$C_chao1$vnfH)
rownames(chao1) = GSV_alpha_div_cluster51$C.env$Hoosfield.ID

alpha = list(sobs, shannon, chao1)
names(alpha) = c('Sobs', 'Shannon', 'Chao1')
rm(list = ls(pattern = 'GSV'), sobs, shannon, chao1)

###############################################################################

env_sorted = env[order(env$pH),]
env_sorted = env_sorted[rownames(env_sorted) %in% rownames(alpha$Sobs),]
env_sorted$change_point = c(1:120)
sample_list = data.frame(ID = env_sorted$Hoosfield.ID)
pH_list = env_sorted[,c('change_point', 'pH')]

## Change point analysis
result = data.frame()

for (i in names(alpha)) {
  print(i)
  df = alpha[[i]] %>% as.data.frame()
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
  result = rbind(result, flattened_df)
}

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))


ggplot(result, aes(x = pH, y = Gene, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_gradientn(colours = mypal.pH(256)) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlim(3.5,8.0) +
  ylab('Alpha-diversity')
