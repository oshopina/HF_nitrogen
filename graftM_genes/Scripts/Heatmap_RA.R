library(gsveasyr)
library(ComplexHeatmap)
library(gt)
library(changepoint)

## Read data
env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
gene_ontology = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
env = env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap
env = env[env$Hoosfield.ID != 'H076',]

GSV_gene_count_list = load_graftm_gene_count(
  n_folder = 'graftM_genes/Data/',
  min_num_reads = 1000000,
  changeSampleName = TRUE,
  env.df = env,
  refColumn = env$gsveasy_sample
)

df = GSV_gene_count_list$N_per_million
df = df[env$Hoosfield.ID,]

############################### Change point analysis ##########################
df_mat = df |> as.matrix()

cpt = geomcp(df_mat)
plot(cpt)

dist_cpt = cpt.meanvar(
  distance(cpt),
  method = "PELT",
  penalty = 'CROPS',
  pen.value = c(5, 500)
)

pen.value.full(dist_cpt)
dist_var = cpts.full(dist_cpt)
tail(dist_var)
plot(dist_cpt, diagnostic = T)
plot(dist_cpt, ncpts = 1)

ang_cpt = cpt.meanvar(
  angle(cpt),
  method = "PELT",
  penalty = 'CROPS',
  pen.value = c(5, 500)
)

pen.value.full(ang_cpt)
ang_var = cpts.full(ang_cpt)
tail(ang_var)
plot(ang_cpt, diagnostic = T)
plot(ang_cpt, ncpts = 2)


## Set up color palette and scale gene count data
my_palette = colorRampPalette(c('white', 'black'))
df_scaled = t(scale(sqrt(df)))
df_scaled = df_scaled[,rownames(env)]
df_scaled = df_scaled[gene_ontology$Gene,]

## Define color palette for pH levels and names for pH levels
col_fun = circlize::colorRamp2(
  c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
  c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a", "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2")
)
names_for_pH = c(3.7, rep("", 17), 4, rep("", 15), 4.5, rep("", 15), 5, rep("", 13), 5.5, rep("", 7), 6, rep("", 10), 
                 6.5, rep("", 12), 7, rep("", 13), 7.5, rep("", 5), 8.0)

# Create HeatmapAnnotation
ha = HeatmapAnnotation(
  pH = env$pH,
  empty = anno_empty(border = FALSE, height = unit(2, "mm")),
  pH_labels = anno_text(names_for_pH, rot = 0),
  col = list(pH = col_fun),
  show_legend = FALSE,
  gp = gpar(col = "black")
)

ha2 = HeatmapAnnotation(
  change_point_dist = anno_lines(
    data.set(dist_cpt),
    axis_param = list(side = 'right', at = c(9, 12))
  ),
  change_point_ang = anno_lines(
    data.set(ang_cpt),
    axis_param = list(side = 'right', at = c(0.08, 0.18))
  ),
  height = unit(3, "cm"), 
  show_annotation_name = F
)

# Create heatmap
Heatmap(
  df_scaled,
  row_order = rownames(df_scaled),
  column_order = env$Hoosfield.ID,
  col = my_palette(100),
  show_column_names = FALSE,
  row_split = factor(gene_ontology$Function, levels = unique(gene_ontology$Function)),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10),
  bottom_annotation = ha,
  top_annotation = ha2,
  show_heatmap_legend = FALSE
)

change_points_dist = cpts(dist_cpt, 1)[!is.na(cpts(dist_cpt, 1))]

#creating vertical line for each change point in distance
for(i in change_points_dist) {
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
    x = unit(0, "npc"),
    y = unit(0.5, "npc"),
    rot = 90,
    vjust = -0.1, 
    gp = gpar(fontsize = 9)
  )
})

change_points_ang = cpts(ang_cpt, 2)[!is.na(cpts(ang_cpt, 2))]

#creating vertical line for each change point in angle
for(i in change_points_ang) {
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
    x = unit(0, "npc"),
    y = unit(0.5, "npc"),
    rot = 90,
    vjust = -0.1,
    gp = gpar(fontsize = 9)
  )
})


## ANOVA for each gene
anova = auto_aov_fixed(df, ~ pH, env)$Results
anova = subset(anova, str_detect(Parameter, 'pH'))[, c('Data', 'F_value', 'p_value', 'Signif')]
rownames(anova) = anova$Data
anova = anova[gene_ontology$Gene,]
anova$F_value = round(anova$F_value, digits = 2)
anova$p_value = round(anova$p_value, digits = 4)

## Save ANOVA results to a Word document
#gtsave(gt::gt(anova), 'graftM_genes/Results/anova_ra.docx')
