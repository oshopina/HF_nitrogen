library(gsveasyr)
library(ComplexHeatmap)
library(gt)

## Read data
env = read.csv('graftM_genes/Data/mag_env_for_shotgun_samples.csv')
gene_ontology = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
GSV_gene_count_list = load_graftm_gene_count(
  n_folder = 'graftM_genes/Data/',
  min_num_reads = 1000000,
  changeSampleName = TRUE,
  env.df = env,
  refColumn = env$gsveasy_sample
)
genes_abs = GSV_gene_count_list$N_per_million/1000000 * GSV_gene_count_list$env$Bact_CN_gsoil

## ANOVA for each gene
anova = auto_aov_fixed(genes_abs, ~ pH, GSV_gene_count_list$env)$Results
anova = subset(anova, str_detect(Parameter, 'pH'))[, c('Data', 'F_value', 'p_value', 'Signif')]
rownames(anova) = anova$Data
anova = anova[gene_ontology$Gene,]
anova$F_value = round(anova$F_value, digits = 2)
anova$p_value = round(anova$p_value, digits = 4)

## Save ANOVA results to a Word document
#gtsave(gt(anova), 'graftM_genes/Figures/anova_qPCR.docx')

## Sort environment data by pH
env_sorted = GSV_gene_count_list$env[order(GSV_gene_count_list$env$pH),]

## Set up color palette and scale gene count data
my_palette = colorRampPalette(c('white', 'black'))
df_scaled = t(scale(sqrt(genes_abs)))
df_scaled = df_scaled[,rownames(env_sorted)]
colnames(df_scaled) = env_sorted$Hoosfield.ID
df_scaled = df_scaled[gene_ontology$Gene,]

## Define color palette for pH levels and names for pH levels
col_fun = circlize::colorRamp2(
  c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
  c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a", "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2")
)
names_for_pH = c(3.7, rep("", 18), 4, rep("", 15), 4.5, rep("", 15), 5, rep("", 13), 5.5, rep("", 7), 6, rep("", 10), 
                 6.5, rep("", 12), 7, rep("", 13), 7.5, rep("", 6), 8.0)

# Create HeatmapAnnotation
ha = HeatmapAnnotation(
  pH_labels = anno_text(names_for_pH, rot = 0),
  empty = anno_empty(border = FALSE, height = unit(3, "mm")),
  pH = env_sorted$pH,
  col = list(pH = col_fun),
  show_legend = FALSE,
  gp = gpar(col = "black")
)

# Create heatmap
Heatmap(
  df_scaled,
  row_order = rownames(df_scaled),
  column_order = env_sorted$Hoosfield.ID,
  col = my_palette(100),
  show_column_names = F,
  row_split = factor(gene_ontology$Function, levels = unique(gene_ontology$Function)),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10),
  top_annotation = ha,
  show_heatmap_legend = FALSE
)
