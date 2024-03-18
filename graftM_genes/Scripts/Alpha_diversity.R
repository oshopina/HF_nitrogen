library(gsveasyr)
library(ComplexHeatmap)
library(changepoint)
library(gt)
library(ggplot2)

## Read data
env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
GSV_gene_count_list = load_graftm_gene_count(
  n_folder = 'graftM_genes/Data/',
  min_num_reads = 1000000,
  changeSampleName = TRUE,
  env.df = env,
  refColumn = env$gsveasy_sample
)

## Sort data
env_sorted = GSV_gene_count_list$env[order(GSV_gene_count_list$env$pH),] |> na.omit()
df_sorted = GSV_gene_count_list$N_per_million[rownames(env_sorted),]
rownames(env_sorted) = env_sorted$Hoosfield.ID
rownames(df_sorted) = env_sorted$Hoosfield.ID

## Alpha diversity metrics
shannon <- diversity(df_sorted, index = "shannon")
sobs <- rowSums(decostand(df_sorted, method = 'pa'))
chao <- estimateR(round(df_sorted))[2,]
alpha <- data.frame(shannon, sobs, chao, env_sorted$pH, samples = names(shannon))
colnames(alpha) <- c('shannon', 'sobs', 'chao', 'pH', 'samples')

## Change point analysis

cpt = cpt.meanvar(
  alpha$shannon,
  method = "PELT",
  penalty = 'CROPS',
  pen.value = c(5, 500)
)

pen.value.full(cpt)
var = cpts.full(cpt)
tail(var)
plot(cpt, diagnostic = T)
plot(cpt, ncpts = 2)
cp = readline('Number of changepoint for distance data: ') |> as.numeric()
change_points = cpts(cpt, cp)[!is.na(cpts(cpt, cp))]

## ANOVA
anova = auto_aov_fixed(alpha[,1:3], ~ pH, env_sorted)$Results
anova = subset(anova, str_detect(Parameter, 'pH'))[, c('Data', 'F_value', 'p_value', 'Signif')]
rownames(anova) = anova$Data
anova$F_value = round(anova$F_value, digits = 2)
anova$p_value = round(anova$p_value, digits = 4)

## Save ANOVA results to a Word document
# gtsave(gt::gt(anova), 'singleM_community/Results/anova_alpha.docx')

## Plot alpha diversity
mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))

change_point_position = alpha$pH[change_points]
ggplot(alpha, aes(y = shannon, x = pH, color = pH)) +
  geom_point(size = 2.5) +
  geom_line() +
  scale_color_gradientn(colours = mypal.pH(256)) +
  theme_classic() +
  theme(legend.position = 'none') +
  geom_vline(xintercept = change_point_position, color = 'red', size = 1) +
  annotate("text", x = change_point_position - 0.07, y = 3.1, label = 'Change point', 
           angle = 90) +
  annotate("text", x = 4, y = 3.25, label = 'p-value < 0.001 (ANOVA)') +
  ylab('Shannon ')
