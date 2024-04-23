library(gsveasyr)
library(ggplot2)

env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
load_alpha_diversity('singleM_genes/Data/', env.df = env,
                        changeSampleName = T, refColumn = env$gsveasy_sample, 
                        clustlvls = seq(50, 100, by = 5))
load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
                        changeSampleName = T, refColumn = env$gsveasy_sample, 
                        clustlvls = 50)

## Get rarefraction levels 
r_levels = c()
for (i in names(GSV_rarefied_gsvtables50)) {
  level = rowSums(GSV_rarefied_gsvtables50[[i]])[1]
  r_levels[i] = level
}
names(r_levels)[which(names(r_levels) == 'narG')] = 'narG_nxrA'
names(r_levels)[which(names(r_levels) == 'narH')] = 'narH_nxrB'

# write.csv2(r_levels, 'singleM_genes/Data/rarefraction_levels.csv')
rm(GSV_rarefied_gsvtables50)
################# Create separate tables for each gene #########################
all_tables = mget(ls(pattern = 'GSV'))
genes = colnames(all_tables$GSV_alpha_div_cluster100$C_shannon)
rm(list = ls(pattern = 'GSV'))

alpha_by_gene = list()
for (j in genes) {
  df = data.frame(matrix(nrow = 120, ncol = 0))
  for (i in names(all_tables)) {
    alpha = all_tables[[i]][['C_shannon']][[j]]
    df = cbind(df, alpha)
  }
  colnames(df) = c(100, seq(50, 95, by = 5))
  rownames(df) = rownames(all_tables$GSV_alpha_div_cluster100$C_shannon)
  alpha_by_gene[[j]] = df
}

################# Create a plot for each gene ##################################

gene_plots = list()
for (j in names(alpha_by_gene)) {
  print(j)
  df = alpha_by_gene[[j]]
  raref = r_levels[j]
  if (all(is.na(df)) == F) {
    df = df[rownames(env), ]
    df$pH = env$pH
    df = pivot_longer(df, cols = 1:11, names_to = 'Cluster')
    df$Cluster = factor(df$Cluster, levels = seq(50, 100, by = 5))
    
    p = ggplot(df, aes(x = pH, y = value, col = Cluster)) +
      geom_point() +
      geom_smooth(aes(fill = Cluster), method = loess, col = 'black')
    
    ggplot_data <- ggplot_build(p)
    predicted_values <- ggplot_data$data[[2]]
    
    # Initialize an empty correlation matrix
    cor_matrix <- matrix(NA, nrow = 11, ncol = 11)
    
    # Loop through all pairwise combinations of groups
    for (group1 in 1:11) {
      for (group2 in 1:11) {
        # Subset the data for the two groups
        values1 <-
          predicted_values[predicted_values$group == group1,]
        values2 <-
          predicted_values[predicted_values$group == group2,]
        
        # Calculate the correlation coefficient and store it in the matrix
        cor_matrix[group1, group2] <-
          cor(values1[['y']], values2[['y']])
        if(is.na(cor_matrix[group1, group2])) {
          cor_matrix[group1, group2] = 0
        }
        cor_matrix[group2, group1] <-
          cor_matrix[group1, group2]  # for symmetry
      }
    }
    
    p2 <- ggplot(df, aes(x = pH, y = value, col = Cluster)) +
      geom_point(show.legend = F) +
      geom_smooth(
        aes(fill = Cluster),
        method = loess,
        col = 'black',
        show.legend = FALSE
      ) +
      annotate(
        geom = 'text',
        x = 6,
        y = min(df$value, na.rm = TRUE) + 0.05,
        label = paste0(
          'ρ = ',
          round(mean(cor_matrix), digits = 2)),
        size = 5,
        fontface = 'bold'
      ) +
      theme_classic(base_size = 17) +
      ggtitle(paste0(j, ' (r.level = ',raref, ')')) +
      ylab('Shannon')
    
    gene_plots[[j]] = p2
  }
}

gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% names(gene_plots), ]
gene_plots = gene_plots[gene_order$Gene]

library(cowplot)
combined_plot2 = plot_grid(plotlist = gene_plots, ncol = 9)
# ggsave('singleM_genes/Figures/5_levels1.svg', device = 'svg', combined_plot2, width = 30, height = 17)

dummy_plot = ggplot(df, aes(x = pH, y = value, col = Cluster)) +
  geom_point(size = 5) +
  theme_classic(base_size = 20) +
  theme(legend.position="bottom") +
  guides(color = guide_legend(nrow = 1)) +
  labs(color = "Clustering level")

# ggsave('singleM_genes/Figures/5_legend.svg', device = 'svg', dummy_plot)

# pdf('singleM_genes/Figures/5_levels.pdf')
# for (i in seq_along(gene_plots)) {
#   print(gene_plots[[i]])
# }
# dev.off()
            