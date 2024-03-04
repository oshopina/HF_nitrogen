library(gsveasyr)
library(ComplexHeatmap)
library(changepoint)
library(changepoint.geo)
library(gt)

################################ Read data #####################################
env <- read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) <- env$Hoosfield.ID
env <- env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap
env = env[env$Hoosfield.ID != 'H076',]

## Load singleM community data
load_singlem_cmm('singleM_community/Data/otu_condensed_table.csv', env_df = env)
otu.bac.nr <- otu.bac.nr[env$gsveasy_sample,] ## sort otu table by env table
rownames(otu.bac.nr) <- rownames(env)

## PERMANOVA
perm = vegan::adonis2(otu.bac.nr ~ pH, env)
# gtsave(gt(perm), 'singleM_community/Results/PERMANOVA.docx')

## Calculate distances
hellinger_diversity = function(otu_table) {
  dist_matrix = dist(vegan::decostand(otu_table, method = 'hellinger'),
                     method = 'euc') %>% as.matrix() 
  return(dist_matrix)
}
beta_diversity = hellinger_diversity(otu.bac.nr)

############################### Change point analysis #########################
cpt = geomcp(beta_diversity)
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


################################ Build heatmap #################################
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

change_points_dist = cpts(dist_cpt, 2)[!is.na(cpts(dist_cpt, 2))]

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
