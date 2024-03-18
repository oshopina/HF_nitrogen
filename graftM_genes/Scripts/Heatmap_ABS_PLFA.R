library(gsveasyr)
library(ComplexHeatmap)
library(gt)
library(openxlsx)

## Read data
env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
env = env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap
env = env[env$Hoosfield.ID != 'H076',]

gene_ontology = read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')

plfa = read.xlsx('graftM_genes/Data/PLFA_indicators.xlsx')
rownames(plfa) = plfa$X1
plfa = plfa[plfa$X1 != 'H097',]
plfa = plfa[env$Hoosfield.ID,]
plfa = plfa[!is.na(rowSums(plfa[,-1])),]


GSV_gene_count_list = load_graftm_gene_count(
  n_folder = 'graftM_genes/Data/',
  min_num_reads = 1000000,
  changeSampleName = TRUE,
  env.df = env,
  refColumn = env$gsveasy_sample
)

genes_abs = GSV_gene_count_list$N_per_million
genes_abs = genes_abs[env$Hoosfield.ID,]
genes_abs = genes_abs[rownames(plfa),]

env = env[rownames(genes_abs),]

## According to literature there is 2.0*10-5 mg of bacterial PLFA per cell so calculating the copy number
plfa_CN = (plfa$biomass - plfa$fungi)/0.00002
genes_abs = genes_abs/1000000 * plfa_CN

############################### Change point analysis ##########################

df_mat = sqrt(log2(genes_abs + 1)) |> as.matrix()
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
plot(dist_cpt, ncpts = 2)
######################### WRITE CHOSEN NUMBER OF CHANGEPOINTS FOR ANGLE AND DISTANCE
cp_dist = readline('Number of changepoint for distance data: ') |> as.numeric()

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
######################### WRITE CHOSEN NUMBER OF CHANGEPOINTS FOR ANGLE AND DISTANCE
cp_angle = readline('Number of changepoint for angle data: ') |> as.numeric()

## ANOVA for each gene
anova = auto_aov_fixed(genes_abs, ~ pH, env)$Results
anova = subset(anova, str_detect(Parameter, 'pH'))[, c('Data', 'F_value', 'p_value', 'Signif')]
rownames(anova) = anova$Data
anova = anova[gene_ontology$Gene,]
anova$F_value = round(anova$F_value, digits = 2)
anova$p_value = round(anova$p_value, digits = 4)
anova$Signif <- gsub("\\*+", "*", anova$Signif)

## Save ANOVA results to a Word document
#gtsave(gt::gt(anova), 'graftM_genes/Results/anova_ra.docx')

################################# Heatmap #####################################
## Set up color palette and scale gene count data
my_palette = colorRampPalette(c('white', 'black'))
genes_abs_scaled = t(scale(sqrt(log2(genes_abs + 1))))
genes_abs_scaled = genes_abs_scaled[,rownames(env)]
genes_abs_scaled = genes_abs_scaled[gene_ontology$Gene,]

medians = apply(genes_abs[,gene_ontology$Gene], 2, function(x){
  median(x)
})
median_col = circlize::colorRamp2(c(0, 80), c('white', 'aquamarine4'))

## Define color palette for pH levels and names for pH levels
col_fun <- circlize::colorRamp2(
  c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
  c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a", "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2")
)
names_for_pH = c(3.7, rep("", 11), 4, rep("", 15), 4.5, rep("", 13), 5, rep("", 10), 5.5, rep("", 6), 6, rep("", 8), 
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
    axis = F
  ),
  change_point_ang = anno_lines(
    data.set(ang_cpt),
    axis = F
  ),
  height = unit(3, "cm"), 
  show_annotation_name = F
)

ha_c <- rowAnnotation(
  Median = medians,
  col = list(Median = median_col),
  show_legend = F, 
  show_annotation_name = F
)

ha_f <- rowAnnotation(
  anova = anno_text(anova$Signif)
)

lgd = Legend(title = 'Median frequency', col_fun = median_col, at = c(0, 30, 80),
             direction = 'horizontal')

# Create heatmap
Heatmap(
  genes_abs_scaled,
  row_order = rownames(genes_abs_scaled),
  column_order = env$Hoosfield.ID,
  col = my_palette(100),
  show_column_names = FALSE,
  row_split = factor(gene_ontology$Function, levels = unique(gene_ontology$Function)),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10),
  bottom_annotation = ha,
  top_annotation = ha2,
  left_annotation = ha_c,
  right_annotation = ha_f,
  show_heatmap_legend = FALSE
)

draw(lgd, x = unit(0.10, "npc"), y = unit(0.95, "npc"))
change_points_dist = cpts(dist_cpt, cp_dist)[!is.na(cpts(dist_cpt, cp_dist))]

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

change_points_ang = cpts(ang_cpt, cp_angle)[!is.na(cpts(ang_cpt, cp_angle))]

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



