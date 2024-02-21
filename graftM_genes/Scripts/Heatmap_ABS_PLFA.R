library(gsveasyr)
library(ComplexHeatmap)
library(gt)
library(openxlsx)

## Read data
env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
gene_ontology = read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
plfa = read.xlsx('graftM_genes/Data/PLFA_indicators.xlsx')
rownames(plfa) = plfa$X1
GSV_gene_count_list = load_graftm_gene_count(
  n_folder = 'graftM_genes/Data/',
  min_num_reads = 1000000,
  changeSampleName = TRUE,
  env.df = env,
  refColumn = env$gsveasy_sample
)
plfa = plfa[GSV_gene_count_list$env$Hoosfield.ID,]
plfa = plfa[!is.na(rowSums(plfa[,-1])),]
genes_abs = GSV_gene_count_list$N_per_million
genes_abs = genes_abs[env$Hoosfield.ID,]
genes_abs = genes_abs[rownames(plfa),]
env = env[rownames(genes_abs),]

## According to literature there is 2.0*10-5 mg of bacterial PLFA per cell so calculating the copy number
plfa_CN = (plfa$biomass - plfa$fungi)/0.00002
genes_abs = genes_abs/1000000 * plfa_CN

############################### Change point analysis ##########################

cpt <- cpt.meanvar(t(genes_abs), method = "PELT", minseglen = 20)
names(cpt) <- colnames(genes_abs)
points <- lapply(cpt, cpts)
points <- Filter(function(x) length(x) > 0, points)
flattened_vector <- unlist(points) ## Flatten the list into a vector
frequency_table <- table(flattened_vector) ## Compute the frequency of each number
frequency_genes_abs <- data.frame(Number = names(frequency_table), Frequency = as.numeric(frequency_table)) ## Convert frequency_table to a data frame
## Fill samples for which change points were not detected with 0s
missing_data <- data.frame(Number = 1:104) 
complete_data <- merge(missing_data, frequency_genes_abs, all.x = TRUE)
complete_data$Frequency[is.na(complete_data$Frequency)] <- 0
rownames(complete_data) <- rownames(genes_abs)

## Sort environment data by pH
env_sorted = env[order(env$pH),]

## Set up color palette and scale gene count data
my_palette = colorRampPalette(c('white', 'black'))
genes_abs_scaled = t(scale(sqrt(genes_abs)))
genes_abs_scaled = genes_abs_scaled[,rownames(env_sorted)]
genes_abs_scaled = genes_abs_scaled[gene_ontology$Gene,]

## Define color palette for pH levels and names for pH levels
col_fun = circlize::colorRamp2(
  c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
  c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a", "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2")
)
names_for_pH = c(3.7, rep("", 12), 4, rep("", 16), 4.5, rep("", 13), 5, rep("", 8), 5.5, rep("", 6), 6, rep("", 8), 
                 6.5, rep("", 12), 7, rep("", 13), 7.5, rep("", 6), 8.0)

# Create HeatmapAnnotation
ha = HeatmapAnnotation(
  pH = env_sorted$pH,
  empty = anno_empty(border = FALSE, height = unit(2, "mm")),
  pH_labels = anno_text(names_for_pH, rot = 0),
  col = list(pH = col_fun),
  show_legend = FALSE,
  gp = gpar(col = "black")
)

ha2 <- HeatmapAnnotation(
  `change point frequency` = anno_barplot(complete_data$Frequency, bar_width = 1, axis_param = list(side = 'right')),
  height = unit(2, "cm"),
  annotation_name_rot = 90,
  annotation_label = 'Change\n point\n frequency',
  annotation_name_offset = unit(0.3, 'cm'),
  annotation_name_side = 'left',
  annotation_name_gp = gpar(fontsize = 9)
)

# Create heatmap
Heatmap(
  genes_abs_scaled,
  row_order = rownames(genes_abs_scaled),
  column_order = env_sorted$Hoosfield.ID,
  col = my_palette(100),
  show_column_names = FALSE,
  row_split = factor(gene_ontology$Function, levels = unique(gene_ontology$Function)),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10),
  bottom_annotation = ha,
  top_annotation = ha2,
  show_heatmap_legend = FALSE
)

## ANOVA for each gene
anova = auto_aov_fixed(genes_abs, ~ pH, env)$Results
anova = subset(anova, str_detect(Parameter, 'pH'))[, c('Data', 'F_value', 'p_value', 'Signif')]
rownames(anova) = anova$Data
anova = anova[gene_ontology$Gene,]
anova$F_value = round(anova$F_value, digits = 2)
anova$p_value = round(anova$p_value, digits = 4)

## Save ANOVA results to a Word document
# gtsave(gt::gt(anova), 'graftM_genes/Results/anova_PLFA.docx')

# write.xlsx(genes_abs, 'graftM_genes/Results/PLFA_ABS_counts.xlsx', rowNames = T)
