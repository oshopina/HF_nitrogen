library(gsveasyr)
library(ggplot2)

env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
load_alpha_diversity('singleM_genes/Data/', env.df = env,
                        changeSampleName = T, refColumn = env$gsveasy_sample, 
                        clustlvls = 50:100)
load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
                        changeSampleName = T, refColumn = env$gsveasy_sample, 
                        clustlvls = 50)

## Get rarefraction levels 
r_levels = c()
for (i in names(GSV_rarefied_gsvtables50)) {
  level = rowSums(GSV_rarefied_gsvtables50[[i]])[1]
  r_levels[i] = level
}
################# Create separate tables for each gene #########################
all_tables = mget(ls(pattern = 'GSV'))
rm(list = ls(pattern = 'GSV'))
genes = colnames(all_tables$GSV_alpha_div_cluster100$C_shannon)

alpha_by_gene = list()
for (j in genes) {
  df = data.frame(matrix(nrow = 120, ncol = 0))
  for (i in names(all_tables)) {
    alpha = all_tables[[i]][['C_shannon']][[j]]
    df = cbind(df, alpha)
  }
  colnames(df) = c(100, 50:99)
  rownames(df) = rownames(all_tables$GSV_alpha_div_cluster100$C_shannon)
  alpha_by_gene[[j]] = df
}

################# Create a plot for each gene ##################################

gene_plots = list()
for (j in genes) {
  print(j)
  df = alpha_by_gene[[j]]
  raref = r_levels[j]
  if (all(is.na(df)) == F) {
    df = df[rownames(env), ]
    df$pH = env$pH
    df = pivot_longer(df, cols = 1:51, names_to = 'Cluster')
    df$Cluster = factor(df$Cluster, levels = 50:100)
    
    p = ggplot(df, aes(x = pH, y = value, col = Cluster)) +
      geom_point() +
      geom_smooth(aes(fill = Cluster), method = loess, col = 'black')
    
    ggplot_data <- ggplot_build(p)
    predicted_values <- ggplot_data$data[[2]]
    
    # Initialize an empty correlation matrix
    cor_matrix <- matrix(NA, nrow = 51, ncol = 51)
    
    # Loop through all pairwise combinations of groups
    for (group1 in 1:51) {
      for (group2 in 1:51) {
        # Subset the data for the two groups
        values1 <-
          predicted_values[predicted_values$group == group1,]
        values2 <-
          predicted_values[predicted_values$group == group2,]
        
        # Calculate the correlation coefficient and store it in the matrix
        cor_matrix[group1, group2] <-
          cor(values1[['y']], values2[['y']])
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
        x = 7,
        y = min(df$value, na.rm = TRUE),
        label = paste0(
          'Mean correlation = ',
          round(mean(cor_matrix), digits = 2),
          '+-',
          round(sd(cor_matrix), digits = 2)
        )
      ) +
      theme_classic(base_size = 15) +
      ggtitle(paste0(j, ' (',raref, ')')) +
      ylab('Shannon')
    
    gene_plots[[j]] = p2
  }
}

library(patchwork)
gene_order = openxlsx::read.xlsx('graftM_genes/Data/nitrogen_genes.xlsx')
gene_order = gene_order[gene_order$Gene %in% names(gene_plots), ]
gene_plots = gene_plots[gene_order$Gene]

combined_plot = wrap_plots(gene_plots, ncol = 7)
# ggsave('singleM_genes/Figures/all_levels.png', combined_plot, width = 21, height = 21)

# pdf('singleM_genes/Figures/all_levls.pdf')
# for (i in seq_along(gene_plots)) {
#   print(gene_plots[[i]])
# }
# dev.off()
            