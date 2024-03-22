library(gsveasyr)
library(changepoint)
library(changepoint.geo)
library(vegan)
library(circlize)
library(ComplexHeatmap)

env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
env <- env[order(env$pH),]
load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
                        changeSampleName = T, refColumn = env$gsveasy_sample, 
                        clustlvls = c(80))
genes = unique(names(GSV_rarefied_gsvtables80))

gc()

for (ngene in 1:length(genes)) {
j = genes[ngene]
print(j)
df = GSV_rarefied_gsvtables80[[j]] |> as.data.frame()
env_temp = env[rownames(env) %in% rownames(df),]
df = df[env_temp$Hoosfield.ID,]

############################### Change point analysis ##########################
df_mat = df |> as.matrix()

cpt = geomcp(df_mat)

dist_cpt = cpt.meanvar(
  distance(cpt),
  method = "PELT",
  penalty = 'CROPS',
  pen.value = c(5, 500),
  minseglen = 20
)

dist_var = cpts.full(dist_cpt)
######################### WRITE CHOSEN NUMBER OF CHANGEPOINTS FOR ANGLE AND DISTANCE
if(all(is.na(dist_var))) next
cp_dist = sum(!is.na(dist_var[nrow(dist_var) - 2,]))

ang_cpt = cpt.meanvar(
  angle(cpt),
  method = "PELT",
  penalty = 'CROPS',
  pen.value = c(5, 500),
  minseglen = 20
)

ang_var = cpts.full(ang_cpt)
######################### WRITE CHOSEN NUMBER OF CHANGEPOINTS FOR ANGLE AND DISTANCE
cp_angle = sum(!is.na(ang_var[nrow(ang_var) - 2,]))

max_values = apply(df, 2, max)
max_values = order(max_values, decreasing = T)
df_plot = df[, max_values[1:50]]

## ANOVA for each gene
anova = auto_aov_fixed(df_plot, ~ pH, env_temp)$Results
anova = subset(anova, str_detect(Parameter, 'pH'))[, c('Data', 'F_value', 'p_value', 'Signif')]
rownames(anova) = anova$Data

gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
anova$F_value = round(anova$F_value, digits = 2)
anova$p_value = round(anova$p_value, digits = 4)
anova$Signif <- gsub("\\*+", "*", anova$Signif)

################################# Heatmap #####################################
## Set up color palette and scale gene count data
my_palette = colorRampPalette(c(
  'white',
  '#eeeeee',
  '#aaaaaa',
  '#444444',
  '#3a3a3a',
  '#2d2d2d',
  'black'
))
df_scaled = t(scale(sqrt(df_plot)))
df_scaled = df_scaled[, rownames(env_temp)]

medians = apply(df_plot, 2, function(x) {
  round(mean(x), digits = 1)
})
median_col = colorRamp2(c(0, max(medians)), c('white', 'aquamarine4'))

## Define color palette for pH levels and names for pH levels
col_fun = colorRamp2(
  c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
  c(
    "#9e0142",
    "#d53e4f",
    "#f46d43",
    "#fdae61",
    "#fee08a",
    "#e6f598",
    "#aadda4",
    "#66a2a5",
    "#3288ad",
    "#5e4fa2"
  )
)
# Get unique pH values
unique_pH <- sort(unique(env_temp$pH))
names_for_pH <- rep("", length(env_temp$pH))
approx_positions <- c(3.7, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.5)
# Iterate over each approximate pH value
for (pH_value in approx_positions) {
  # Find the index in unique_pH that is closest to the approximate pH value
  closest_index <- which.min(abs(unique(env_temp$pH) - pH_value))
  
  # Find the index in env_temp$pH that corresponds to the first occurrence of the unique pH value
  first_occurrence_index <- which(env_temp$pH == unique(env_temp$pH)[closest_index])[1]
  
  # Set the corresponding position in names_for_pH to the approximate pH value
  names_for_pH[first_occurrence_index] <- pH_value
}

# Create HeatmapAnnotation
ha = HeatmapAnnotation(
  pH = env_temp$pH,
  empty = anno_empty(border = FALSE, height = unit(2, "mm")),
  pH_labels = anno_text(names_for_pH, rot = 0),
  col = list(pH = col_fun),
  show_legend = FALSE,
  gp = gpar(col = "black")
)

ha2 = HeatmapAnnotation(
  change_point_dist = anno_lines(data.set(dist_cpt),
                                 axis = F),
  change_point_ang = anno_lines(data.set(ang_cpt),
                                axis = F),
  height = unit(3, "cm"),
  show_annotation_name = F
)

ha_c <- rowAnnotation(
  Median = medians,
  col = list(Median = median_col),
  show_legend = F,
  show_annotation_name = F
)

ha_f <- rowAnnotation(anova = anno_text(anova$Signif))

lgd = Legend(
  title = 'Mean\nfrequency',
  col_fun = median_col,
  at = c(0, max(medians)),
  direction = 'vertical',
  title_gp = gpar(fontsize = 9, fontface = 'bold')
)

grid.newpage()

# Create heatmap
heatmap_name <- paste("heatmap_", j, sep = "")
draw(Heatmap(
  df_scaled,
  column_order = env_temp$Hoosfield.ID,
  col = my_palette(100),
  show_column_names = FALSE,
  bottom_annotation = ha,
  top_annotation = ha2,
  left_annotation = ha_c,
  right_annotation = ha_f,
  show_heatmap_legend = FALSE,
  column_title = j
))

draw(lgd, x = unit(0.04, "npc"), y = unit(0.95, "npc"))

change_points_dist = cpts(dist_cpt, cp_dist)[!is.na(cpts(dist_cpt, cp_dist))]

#creating vertical line for each change point in distance
for (i in change_points_dist) {
  decorate_annotation("change_point_dist", {
    grid.lines(
      x = unit(c(i, i), 'native'),
      y = unit(c(min(
        data.set(dist_cpt)
      ), max(
        data.set(dist_cpt)
      )), 'native'),
      gp = gpar(col = "red", lwd = 3)
    )
  })
}

#Adding title
decorate_annotation("change_point_dist", {
  grid.text(
    "Mean\nchange\npoint",
    x = unit(1, "npc"),
    y = unit(0.5, "npc"),
    rot = -90,
    vjust = -0.1,
    gp = gpar(fontsize = 9)
  )
})

change_points_ang = cpts(ang_cpt, cp_angle)[!is.na(cpts(ang_cpt, cp_angle))]

#creating vertical line for each change point in angle
for (i in change_points_ang) {
  decorate_annotation("change_point_ang", {
    grid.lines(
      x = unit(c(i, i), 'native'),
      y = unit(c(min(
        data.set(ang_cpt)
      ), max(
        data.set(ang_cpt)
      )), 'native'),
      gp = gpar(col = "red", lwd = 3)
    )
  })
}

#Adding title
decorate_annotation("change_point_ang", {
  grid.text(
    "Variance\nchange\npoint",
    x = unit(1, "npc"),
    y = unit(0.5, "npc"),
    rot = -90,
    vjust = -0.1,
    gp = gpar(fontsize = 9)
  )
})
 
p = recordPlot()
plot.new()
png(paste0('singleM_genes/Figures/gene_plots/', j, '.png'), height = 900, width = 700)
print(p)
dev.off()
}



