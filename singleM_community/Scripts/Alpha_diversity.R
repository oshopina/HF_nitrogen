library(gsveasyr)
library(gt)
library(ggplot2)

## Read data
env <- read.csv('graftM_genes/Data/mag_env_for_shotgun_samples.csv')
rownames(env) <- env$Hoosfield.ID
env <- env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap
env = env[env$Hoosfield.ID != 'H076',]

## Load singleM community data
load_singlem_cmm('singleM_community/Data/otu_condensed_table.csv', env_df = env)
otu.bac.nr <- otu.bac.nr[env$gsveasy_sample,] ## sort otu table by env table
rownames(otu.bac.nr) <- rownames(env)

## Alpha diversity metrics
shannon <- diversity(otu.bac.nr, index = "shannon")
sobs <- rowSums(decostand(otu.bac.nr, method = 'pa'))
chao <- estimateR(round(otu.bac.nr))[2,]
alpha <- data.frame(shannon, sobs, chao, env$pH, samples = names(shannon))
colnames(alpha) <- c('shannon', 'sobs', 'chao', 'pH', 'samples')

## Plot alpha diversity

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))

for (i in colnames(alpha[,1:3])) {
  plot1 = ggplot(alpha, aes(y = !!sym(i), x = pH, color = pH)) +
    geom_point(size = 2.5) +
    scale_color_gradientn(colours = mypal.pH(256)) +
    theme_classic() +
    theme(legend.position = 'none')
  
  ggsave(paste0('singleM_community/Figures/', i, '.png'), plot1, height = 4, width = 5)
}

## ANOVA for each gene
anova = auto_aov_fixed(alpha[,1:3], ~ pH, env)$Results
anova = subset(anova, str_detect(Parameter, 'pH'))[, c('Data', 'F_value', 'p_value', 'Signif')]
rownames(anova) = anova$Data
anova$F_value = round(anova$F_value, digits = 2)
anova$p_value = round(anova$p_value, digits = 4)

## Save ANOVA results to a Word document
# gtsave(gt::gt(anova), 'singleM_community/Results/anova_alpha.docx')
