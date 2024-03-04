library(gsveasyr)
library(Umatrix)
################################ Prepare data ##################################

## Read environmental data
env <- read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) <- env$Hoosfield.ID
env <- env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap

## Load singleM community data
load_singlem_cmm('singleM_community/Data/otu_condensed_table.csv', env_df = env)
otu.bac.nr <- otu.bac.nr[env$gsveasy_sample,] ## sort otu table by env table
rownames(otu.bac.nr) <- rownames(env)

df = otu.bac.nr
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

plotMatrix(
  res$Umatrix,
  res$BestMatches,
  Cls = 1:nrow(df),
  ClsColors = ph.colors$col,
  BmSize = 5,
  YellowCircle = T,
  TransparentContours = T
)

class_manual = iClassification(res$Umatrix, res$BestMatches)

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

class.1 = env[class_manual$Cls == 4,]
summary(class.1)

showMatrix3D(
  res$Umatrix)
