library(gsveasyr)
library(ComplexHeatmap)
library(changepoint)
library(changepoint.geo)

################################ Prepare data ##################################

## Read environmental data
env <- read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) <- env$Hoosfield.ID
env <- env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap

## Load singleM community data
load_singlem_cmm('singleM_community/Data/otu_condensed_table.csv', env_df = env)
otu.bac.nr <- otu.bac.nr[env$gsveasy_sample,] ## sort otu table by env table
rownames(otu.bac.nr) <- rownames(env)


df <- otu.bac.nr[, apply(otu.bac.nr, 2, max) >= 150]  # Keep columns with max value >= 150 (1.5%)

############################### Change point analysis ##########################

df_mat = df |> as.matrix() |> t()
cpt = geomcp(df)
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


######################## Prepare taxonomy labels for heatmap ##################

tax <- taxonomy.all[colnames(df),]
tax <- tax[order(tax$Phylum, tax$Class, tax$Order, tax$Family, tax$Genus),] ## Order taxonomy alphabetically
tax[is.na(tax)] <- 'Unc.'

tax$phylum_label <- tax$Phylum
tax$phylum_label[duplicated(tax$phylum_label)] <- ""
tax$class_label <- paste0('c_', tax$Class)
tax$class_label[duplicated(tax$class_label)] <- ""
tax$order_label <- paste0('o_', tax$Order)
tax$order_label[duplicated(tax$order_label)] <- ""
tax$f_g_label <- paste0('f_', tax$Family, ';g_', tax$Genus)
tax$f_g_label[duplicated(tax$f_g_label)] <- ""
tax$otu_label <- paste0('[', substring(tax$OTU, 2), ']')

df <- df[, tax$OTU]

############################## Heatmap ########################################

## Scale data for heatmap
df_scaled <- t(scale(sqrt(df)))

## Create bottom heatmap annotation (pH ribbon)
col_fun <- circlize::colorRamp2(
  c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
  c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a", "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2")
)
names_for_pH <- c(3.7, rep("", 16), 4, rep("", 15), 4.5, rep("", 15), 5, rep("", 13), 5.5, rep("", 7), 6, rep("", 10), 
                  6.5, rep("", 12), 7, rep("", 13), 7.5, rep("", 6), 8.0)

ha <- HeatmapAnnotation(
  pH = env$pH,
  empty = anno_empty(border = FALSE, height = unit(3, "mm")),
  pH_labels = anno_text(names_for_pH, rot = 0),
  col = list(pH = col_fun),
  show_legend = FALSE,
  gp = gpar(col = "black")
)

## Create taxonomy annotations

ha_c <- rowAnnotation(
  phylum = anno_text(tax$phylum_label, gp = gpar(fontface = 'bold')),
  class = anno_text(tax$class_label),
  order = anno_text(tax$order_label)
)

ha_f <- rowAnnotation(
  otu = anno_text(tax$otu_label),
  family = anno_text(tax$f_g_label)
)

## Create additional annotation for change point frequency
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

## Create Heatmap
my_palette <- colorRampPalette(c('white', 'black'))

Heatmap(
  df_scaled,
  row_order = rownames(df_scaled),
  column_order = env$Hoosfield.ID,
  col = my_palette(100),
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_title = NULL,
  row_split = tax$Phylum,
  top_annotation = ha2,
  bottom_annotation = ha,
  left_annotation = ha_c,
  right_annotation = ha_f,
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


## Perform ANOVA analysis
anova <- auto_aov_fixed(df, ~ pH, env_df = env)$Results
anova = subset(anova, str_detect(Parameter, 'pH'))[, c('Data', 'F_value', 'p_value', 'Signif')]
rownames(anova) = anova$Data
anova = anova[colnames(df),]
anova$F_value = round(anova$F_value, digits = 2)
anova$p_value = round(anova$p_value, digits = 4)

## Save ANOVA results to a Word document
# gtsave(gt::gt(anova), 'singleM_community/Results/anova_ra.docx')
