library(openxlsx)
library(ggplot2)
library(ggpmisc)
library(dplyr)

plfa = read.xlsx('graftM_genes/Results/PLFA_ABS_counts.xlsx', rowNames = T)
plfa_indices = read.xlsx('graftM_genes/Data/PLFA_indicators.xlsx')
rownames(plfa_indices) = plfa_indices$X1
qPCR_all = read.xlsx('graftM_genes/Results/qPCR_ABS_counts.xlsx', rowNames = T)
qPCR = read.csv('graftM_genes/Data/N and 16s and ITS.csv', row.names = 1)
env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID
GSV_gene_count_list = load_graftm_gene_count(
  n_folder = 'graftM_genes/Data/',
  min_num_reads = 1000000,
  changeSampleName = TRUE,
  env.df = env,
  refColumn = env$gsveasy_sample
)
ra = GSV_gene_count_list$N_per_million[env$Hoosfield.ID,]

## Compare qPCR and qPCR estimate

env_sorted = env[order(env$pH),]
qPCR_sorted = qPCR[rownames(env_sorted),] %>% na.omit()
plfa_ind_sorted = plfa_indices[rownames(env_sorted),] %>% na.omit()

plot(qPCR_sorted$pH, qPCR_sorted$X16s.copynumber_g.dry.soil)

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))

temp = env_sorted[rownames(env_sorted) %in% rownames(qPCR_sorted),]
ggplot(qPCR_sorted, aes(y = X16s.copynumber_g.dry.soil, x = temp$pH, color = temp$pH)) +
  geom_point(size = 2.5) +
  scale_color_gradientn(colours = mypal.pH(256)) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('pH') +
  ylab('qPCR')

ggplot(env_sorted, aes(y = Bact_CN_gsoil, x = as.factor(env_sorted$pH), color = pH)) +
  geom_point(size = 2.5) +
  scale_color_gradientn(colours = mypal.pH(256)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_blank()) +
  xlab('pH') +
  ylab('qPCR estimate')

plfa_ind_sorted = plfa_ind_sorted[rownames(plfa_ind_sorted) %in% rownames(qPCR),]
temp = env_sorted[rownames(env_sorted) %in% rownames(plfa_ind_sorted),]
ggplot(plfa_ind_sorted, aes(y = biomass, x = temp$pH, color = temp$pH)) +
  geom_point(size = 2.5) +
  scale_color_gradientn(colours = mypal.pH(256)) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('pH') +
  ylab('PLFA biomass')

## Compare individual genes abundances

qPCR = qPCR[rownames(qPCR) %in% rownames(plfa),]
plfa = plfa[rownames(qPCR),]
qPCR_all = qPCR_all[rownames(qPCR),]
ra = ra[rownames(qPCR),]

qPCR = qPCR[,c(7,8,11,13,14)]
colnames(qPCR) = c('narH_nxrB', 'amoA', 'napA', 'narG_nxrA', 'norB')

## Loop for qPCR estimates
for (i in colnames(qPCR)) {
  
  df = data.frame(qPCR = qPCR[i], estimate = qPCR_all[i])
  colnames(df) = c('qPCR', 'qPCR estimate')
  
  plot1 = ggplot(df, aes(x = qPCR, y = `qPCR estimate`)) +
    geom_point() +
    stat_poly_line() +
    stat_poly_eq() +
    ggtitle(i) +
    theme_classic()
  
  ggsave(paste0('graftM_genes/Figures/qPCR_vs_estimate/', i, '_qPCR.png'), plot1,
         height = 4, width = 5)
}

## Loop for PLFA estimates
for (i in colnames(qPCR)) {
  
  df = data.frame(qPCR = log(qPCR[i]), estimate = log(plfa[i]))
  colnames(df) = c('qPCR', 'PLFA estimate')
  
  plot1 = ggplot(df, aes(x = qPCR, y = `PLFA estimate`)) +
    geom_point() +
    stat_poly_line() +
    stat_poly_eq(use_label(c('R2', 'P'))) +
    ggtitle(i) +
    theme_classic()
  
  ggsave(paste0('graftM_genes/Figures/qPCR_vs_estimate/', i, '_PLFA.png'), plot1,
         height = 4, width = 5)
}

## Loop for Relative abundance
for (i in colnames(qPCR)) {
  
  df = data.frame(qPCR = qPCR[i], estimate = ra[i])
  colnames(df) = c('qPCR', 'Relative abundance')
  
  plot1 = ggplot(df, aes(x = qPCR, y = `Relative abundance`)) +
    geom_point() +
    stat_poly_line() +
    stat_poly_eq(use_label(c('R2', 'P'))) +
    ggtitle(i) +
    theme_classic()
  
  ggsave(paste0('graftM_genes/Figures/qPCR_vs_estimate/', i, '_RA.png'), plot1,
         height = 4, width = 5)
}

###################### Compare models ##########################################
env_qpcr = env[rownames(qPCR),]
model = lm(qPCR$X16s.copynumber_g.dry.soil ~ pH + I(pH^2), env_qpcr)
model1 = lm(qPCR$X16s.copynumber_g.dry.soil ~ pH, env_qpcr)
model2 = lm(qPCR$X16s.copynumber_g.dry.soil ~ I(pH^2), env_qpcr)
summary(model)
summary(model1)
summary(model2)

model = lm(plfa_indices$bacteria ~ pH + I(pH^2), env_qpcr)
model1 = lm(plfa_indices$bacteria ~ pH, env_qpcr)
model2 = lm(plfa_indices$bacteria ~ I(pH^2), env_qpcr)
summary(model)
summary(model1)
summary(model2)


plot(qPCR$X16s.copynumber_g.dry.soil ~ pH, env_qpcr)
curve(predict(model, newdata = data.frame(pH = x)), add = T)
