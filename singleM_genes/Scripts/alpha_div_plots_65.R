library(gsveasyr)
library(ggridges)
library(ggplot2)
library(ggpmisc)

env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
load_alpha_diversity('singleM_genes/Data/', env.df = env,
                        changeSampleName = T, refColumn = env$gsveasy_sample, 
                        clustlvls = c(65, 85))


######################### Combine best clustering levels #######################
shannon65 = GSV_alpha_div_cluster65$C_shannon
shannon85 = GSV_alpha_div_cluster85$C_shannon

rm(list = ls(pattern = 'GSV'))
###############################################################################

#### 65 ####
env_sorted = env[order(env$pH), ]
df = shannon65
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

#### 85 ####

df = shannon85
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

