library(gsveasyr)
library(vegan)

# env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
# rownames(env) = env$Hoosfield.ID
# load_rarefied_gsvtables('singleM_genes/Data/', env.df = env,
#                         changeSampleName = T, refColumn = env$gsveasy_sample, 
#                         clustlvls = 50:100)
# 
# 
# hellinger_diversity = function(otu_table) {
#   dist_matrix = dist(vegan::decostand(otu_table, method = 'hellinger'),
#                      method = 'euc') %>% as.matrix() 
#   return(dist_matrix)
# }
# rarefied_tables = mget(ls(pattern = 'GSV'))
# rm(list = ls(pattern = 'GSV'))
# genes = names(rarefied_tables$GSV_rarefied_gsvtables100)

load('../HF_nitrogen_paper/mantel.RData')

matrix_list = list()
for (v in genes) {
  print(v)
  cor_matrix <- matrix(NA, nrow = 51, ncol = 51)
  # Loop through the groups
  for (i in 1:50) {
    group1 <- names(rarefied_tables)[i]
    # Subset the data for the first group
    dist1 <- hellinger_diversity(rarefied_tables[[group1]][[v]])
    # Start the inner loop from i + 1 to avoid redundant comparisons
    for (j in seq(i + 1, 51)) {
      group2 <- names(rarefied_tables)[j]
      # Subset the data for the second group
      dist2 <- hellinger_diversity(rarefied_tables[[group2]][[v]])
      # Print the current comparison
      cat("Comparing:", group1, "and", group2, "\n")
      # Calculate the Mantel correlation coefficient and store it in the matrix
      mantel_r <- mantel(dist1, dist2)
      # Store the result in the upper triangular part of the correlation matrix
      cor_matrix[i, j] <- mantel_r$statistic
      cor_matrix[j, i] <- mantel_r$statistic
    }
    colnames(cor_matrix) <- names(rarefied_tables)
    rownames(cor_matrix) <- names(rarefied_tables)
  }
  matrix_list[[v]] = cor_matrix
}

pr_list <- list()
for (v in genes) {
  print(v)
  cor_matrix <- matrix(NA, nrow = 51, ncol = 51)
  
  # Loop through the groups
  for (i in 1:50) {
    group1 <- names(rarefied_tables)[i]
    # Subset the data for the first group
    dist1 <- hellinger_diversity(rarefied_tables[[group1]][[v]])
    
    # Start the inner loop from i + 1 to avoid redundant comparisons
    for (j in seq(i + 1, 51)) {
      group2 <- names(rarefied_tables)[j]
      # Subset the data for the second group
      dist2 <- hellinger_diversity(rarefied_tables[[group2]][[v]])
      
      # Print the current comparison
      cat("Comparing:", group1, "and", group2, "\n")
      
      tryCatch({
        # Calculate the Mantel correlation coefficient and store it in the matrix
        mantel_r <- protest(dist1, dist2)
        # Store the result in the upper triangular part of the correlation matrix
        cor_matrix[i, j] <- mantel_r$ss
        cor_matrix[j, i] <- mantel_r$ss
      }, error = function(e) {
        # Handle the error and set the matrix elements to 1
        cat("Error occurred:", conditionMessage(e), "\n")
        cor_matrix[i, j] <<- 1
        cor_matrix[j, i] <<- 1
      })
    }
  }
  colnames(cor_matrix) <- names(rarefied_tables)
  rownames(cor_matrix) <- names(rarefied_tables)
  pr_list[[v]] <- cor_matrix
}

