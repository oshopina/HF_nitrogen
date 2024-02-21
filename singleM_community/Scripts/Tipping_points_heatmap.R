library(gsveasyr)
library(ComplexHeatmap)
library(changepoint)
library(gt)

## Read data
env <- read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) <- env$Hoosfield.ID
env <- env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap
env = env[env$Hoosfield.ID != 'H076',]

## Load singleM community data
load_singlem_cmm('singleM_community/Data/otu_condensed_table.csv', env_df = env)
otu.bac.nr <- otu.bac.nr[env$gsveasy_sample,] ## sort otu table by env table
rownames(otu.bac.nr) <- rownames(env)

perm = vegan::adonis2(otu.bac.nr ~ pH, env)
# gtsave(gt(perm), 'singleM_community/Results/PERMANOVA.docx')

## Calculate distances
hellinger_diversity = function(otu_table) {
  dist_matrix = dist(vegan::decostand(otu_table, method = 'hellinger'),
                     method = 'euc') %>% as.matrix() 
  return(dist_matrix)
}
beta_diversity = hellinger_diversity(otu.bac.nr)

## Change point analysis
cpt = cpt.meanvar(beta_diversity, method="PELT", minseglen = 20)
names(cpt) = rownames(otu.bac.nr)

points = lapply(cpt, cpts)
points = Filter(function(x) length(x) > 0, points)

# Flatten the list into a vector
flattened_vector <- unlist(points)

# Compute the frequency of each number
frequency_table <- table(flattened_vector)

# Convert frequency_table to a data frame
frequency_df <- data.frame(Number = names(frequency_table), Frequency = as.numeric(frequency_table))
missing_data = data.frame(Number = 1:117)
complete_data = merge(missing_data, frequency_df, all.x = TRUE)
complete_data$Frequency[is.na(complete_data$Frequency)] = 0
rownames(complete_data) = rownames(beta_diversity)


## Build heatmap
my_palette = colorRampPalette(c('#100d12', '#980011', '#fe8d2f', '#fff1a2'))
col_fun = circlize::colorRamp2(c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
                               c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a",
                                 "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2"))
names_for_pH = c(3.7, rep("", 16), 4, rep("", 15), 4.5, rep("", 15), 5, rep("", 13), 5.5, rep("", 7), 6, rep("", 10), 
                 6.5, rep("", 12), 7, rep("", 13), 7.5, rep("", 6), 8.0)

## Create HeatmapAnnotation
ha = HeatmapAnnotation(
  pH = env$pH,
  empty = anno_empty(border = FALSE, height = unit(3, "mm")),
  pH_labels = anno_text(names_for_pH, rot = 0),
  col = list(pH = col_fun),
  show_legend = FALSE,
  gp = gpar(col = "black")
)

ha1 = rowAnnotation(pH_labels = anno_text(names_for_pH, rot = 0), pH = env$pH,
                    col = list(pH = col_fun), show_legend = F, 
                    show_annotation_name = F, gp = gpar(col = "black"))

ha2 = HeatmapAnnotation(`change point frequency` = anno_barplot(complete_data$Frequency, bar_width = 1, 
                                                                axis_param = list(side = 'right')),
                        height = unit(2, "cm"), annotation_name_rot = 90, annotation_label = 'Change\n point\n frequency',
                        annotation_name_offset = unit(0.3, 'cm'), annotation_name_side = 'left',
                        annotation_name_gp = gpar(fontsize = 9))

Heatmap(
  beta_diversity,
  row_order = rownames(beta_diversity),
  column_order = env$Hoosfield.ID,
  col = my_palette(100),
  show_column_names = F,
  show_row_names = F,
  bottom_annotation = ha,
  left_annotation = ha1,
  top_annotation = ha2,
  show_heatmap_legend = FALSE
)
