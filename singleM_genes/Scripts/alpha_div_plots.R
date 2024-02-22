library(gsveasyr)
library(ggridges)
library(ggplot2)
library(ggpmisc)

env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
best_cluster = openxlsx::read.xlsx('singleM_genes/Results/best_cluster.xlsx')
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

env_sorted = env[order(env$pH), ]
df = shannon
df = df[env_sorted$Hoosfield.ID, ]

anova = auto_aov_fixed(df, ~ pH, env_sorted)$Results

df$pH = env_sorted$pH
df_longer = tidyr::pivot_longer(df,
                                cols = c(1:(ncol(df) - 1)),
                                names_to = 'Gene',
                                values_to = 'Diversity') %>%
  drop_na()

gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% colnames(df), ]
df_longer$Gene = factor(df_longer$Gene, levels = gene_order$Gene)

ggplot(df_longer, aes(x = pH, y = Diversity)) +
  geom_point() +
  geom_smooth() +
  theme_minimal(base_size = 15) +
  facet_wrap(vars(Gene), scales = 'free_y') +
  ylab('Shannon diversity')



