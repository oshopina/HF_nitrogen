library(gsveasyr)
library(kohonen)

################################ Prepare data ##################################

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
df.m = data.matrix(df)

set.seed(4)
grid.size = 10
som.grid = somgrid(xdim = grid.size, ydim = grid.size, topo = 'hexagonal', toroidal = T)
som.model = som(df.m, grid = som.grid)

col_fun <- circlize::colorRamp2(
  c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
  c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a", "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2")
)

ph.colors = env |> select(Hoosfield.ID, pH) |> mutate(col = col_fun(pH))

umap = plot(som.model, type = "dist.neighbours", palette.name = terrain.colors)

colors1 = circlize::colorRamp2(seq(min(umap), max(umap), (max(umap) - min(umap)) / 10), terrain.colors(11))

plot(
  som.model,
  type = "mapping",
  pchs = 20,
  col = ph.colors$col,
  cex = 3,
  bgcol = colors1(umap)) 


som.model.kmeans <- kmeans(sqrt(som.model$codes[[1]]), 3)$cluster

add.cluster.boundaries(som.model, som.model.kmeans)



