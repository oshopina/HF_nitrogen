library(gsveasyr)
library(ComplexHeatmap)
library(changepoint)

env = read.csv('graftM_genes/Data/mag_env_for_shotgun_samples.csv')
rownames(env) = env$Hoosfield.ID
env = env[order(env$pH),]
load_singlem_cmm('singleM_community/Data/otu_condensed_table.csv', env_df = env)
otu.bac.nr = otu.bac.nr[env$gsveasy_sample,]
rownames(otu.bac.nr) = rownames(env)

anova = auto_aov_fixed(otu.bac.nr, ~ pH, env_df = env)$Results
sign_otus = anova[anova$Signif != " ",]
df = otu.bac.nr[,sign_otus$Data]
df = df[, apply(df, 2, max) >= 150] ## 150 = 1.5%

## Change point analysis

cpt = cpt.mean(t(df), method="PELT",
                  penalty = "Manual", pen.value = 100, minseglen = 10)
names(cpt) = colnames(df)

points = lapply(cpt, cpts)
points = Filter(function(x) length(x) > 0, points)

# Flatten the list into a vector
flattened_vector <- unlist(points)

# Compute the frequency of each number
frequency_table <- table(flattened_vector)

# Convert frequency_table to a data frame
frequency_df <- data.frame(Number = names(frequency_table), Frequency = as.numeric(frequency_table))
missing_data = data.frame(Number = 1:120)
complete_data = merge(missing_data, frequency_df, all.x = TRUE)
complete_data$Frequency[is.na(complete_data$Frequency)] = 0
rownames(complete_data) = rownames(df)


tax = taxonomy.all[colnames(df),]
tax = tax[order(tax$Phylum, tax$Class, tax$Order, tax$Family, tax$Genus),]
tax[is.na(tax) == T] = 'Unc.'

tax$phylum_label = tax$Phylum
tax$phylum_label[duplicated(tax$phylum_label)] <- ""
tax$class_label = paste0('c_', tax$Class)
tax$class_label[duplicated(tax$class_label)] <- ""
tax$order_label = paste0('o_', tax$Order)
tax$order_label[duplicated(tax$order_label)] <- ""
tax$f_g_label = paste0('f_', tax$Family, 'g_', tax$Genus)
tax$f_g_label[duplicated(tax$f_g_label)] <- ""
tax$otu_label = paste0('[', substring(tax$OTU,2), ']')

df = df[,tax$OTU]

my_palette = colorRampPalette(c('white', 'black'))
df_scaled = t(scale(sqrt(df)))

## Define color palette for pH levels and names for pH levels
col_fun = circlize::colorRamp2(
  c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
  c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a", "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2")
)
names_for_pH = c(3.7, rep("", 18), 4, rep("", 15), 4.5, rep("", 15), 5, rep("", 13), 5.5, rep("", 7), 6, rep("", 10), 
                 6.5, rep("", 12), 7, rep("", 13), 7.5, rep("", 7), 8.0)


# Create HeatmapAnnotation
ha = HeatmapAnnotation(
  pH_labels = anno_text(names_for_pH, rot = 0),
  empty = anno_empty(border = FALSE, height = unit(3, "mm")),
  pH = env$pH,
  col = list(pH = col_fun),
  show_legend = FALSE,
  gp = gpar(col = "black")
)

ha_c = rowAnnotation(
  phylum = anno_text(tax$phylum_label, gp = gpar(fontface = 'bold')),
  class = anno_text(tax$class_label),
  order = anno_text(tax$order_label)
)

ha_f = rowAnnotation(
  otu = anno_text(tax$otu_label),
  family = anno_text(tax$f_g_label)
)

ha2 = HeatmapAnnotation(`change point frequency` = anno_barplot(complete_data$Frequency, bar_width = 1, 
                                                                axis_param = list(side = 'right')),
                        height = unit(2, "cm"), annotation_name_rot = 90, annotation_label = 'Change\n point\n frequency',
                        annotation_name_offset = unit(0.3, 'cm'), annotation_name_side = 'left',
                        annotation_name_gp = gpar(fontsize = 9))

# Create heatmap
Heatmap(
  df_scaled,
  row_order = rownames(df_scaled),
  column_order = env$Hoosfield.ID,
  col = my_palette(100),
  show_column_names = FALSE,
  show_row_names = F,
  row_title = NULL,
  row_split = tax$Phylum,
  top_annotation = ha2,
  bottom_annotation = ha,
  left_annotation = ha_c,
  right_annotation = ha_f,
  show_heatmap_legend = FALSE
)
