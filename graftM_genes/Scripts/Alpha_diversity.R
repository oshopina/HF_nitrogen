library(gsveasyr)
library(ComplexHeatmap)
library(changepoint)
library(gt)
library(ggplot2)

## Read data
env = read.csv('graftM_genes/Data/mag_env_for_shotgun_samples.csv')
GSV_gene_count_list = load_graftm_gene_count(
  n_folder = 'graftM_genes/Data/',
  min_num_reads = 1000000,
  changeSampleName = TRUE,
  env.df = env,
  refColumn = env$gsveasy_sample
)

## Sort data
env_sorted = GSV_gene_count_list$env[order(GSV_gene_count_list$env$pH),]
df_sorted = GSV_gene_count_list$N_per_million[rownames(env_sorted),]
rownames(env_sorted) = env_sorted$Hoosfield.ID
rownames(df_sorted) = env_sorted$Hoosfield.ID

## Alpha diversity metrics
shannon <- diversity(df_sorted, index = "shannon")
sobs <- rowSums(decostand(df_sorted, method = 'pa'))
chao <- estimateR(round(df_sorted))[2,]
alpha <- data.frame(shannon, sobs, chao, env_sorted$pH, samples = names(shannon))
colnames(alpha) <- c('shannon', 'sobs', 'chao', 'pH', 'samples')

## Plot alpha diversity

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))

for (i in colnames(alpha[,1:3])) {
  plot1 = ggplot(alpha, aes_string(y = i, x = 'pH', color = 'pH')) +
    geom_point(size = 2.5) +
    scale_color_gradientn(colours = mypal.pH(256)) +
    theme_classic() +
    theme(legend.position = 'none')
  
  ggsave(paste0('graftM_genes/Figures/', i, '.png'), plot1, height = 4, width = 5)
}

## ANOVA for each gene
anova = auto_aov_fixed(alpha[,1:3], ~ pH, env_sorted)$Results
anova = subset(anova, str_detect(Parameter, 'pH'))[, c('Data', 'F_value', 'p_value', 'Signif')]
rownames(anova) = anova$Data
anova$F_value = round(anova$F_value, digits = 2)
anova$p_value = round(anova$p_value, digits = 4)

## Save ANOVA results to a Word document
# gtsave(gt::gt(anova), 'graftM_genes/Results/anova_alpha.docx')
