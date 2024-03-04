library(gsveasyr)
library(Umatrix)
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

df.m = df |> as.matrix()
df.m[is.na(df.m)] = 0

col_fun <- circlize::colorRamp2(
  c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
  c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a", "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2")
)
ph.colors = env |> select(Hoosfield.ID, pH) |> mutate(col = col_fun(pH))

set.seed(1)
res = esomTrain(df.m,
                Lines = 100,
                Columns = 100,
                Epochs = 30)

png('C:/Users/Ольга/Dropbox/Olga Shopina/HF_nitrogen/graftM_genes/Figures/esom1.png')
plotMatrix(
  res$Umatrix,
  res$BestMatches,
  Cls = 1:nrow(df),
  ClsColors = ph.colors$col,
  BmSize = 5,
  YellowCircle = T,
  TransparentContours = T
)
dev.off()

class_manual = iClassification(res$Umatrix, res$BestMatches)

png('C:/Users/Ольга/Dropbox/Olga Shopina/HF_nitrogen/graftM_genes/Figures/esom2.png')
plotMatrix(
  res$Umatrix,
  res$BestMatches,
  Cls = class_manual$Cls,
  BmSize = 5,
  YellowCircle = T, 
  DrawLegend = T,
  TransparentContours = T,
  TransparentOcean = T
)
dev.off()


class.1 = env[class_manual$Cls == 3,]
summary(class.1)

showMatrix3D(
  res$Umatrix)
