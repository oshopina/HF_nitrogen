library(openxlsx)
library(ggplot2)
library(ggpmisc)
library(dplyr)

plfa = read.xlsx('graftM_genes/Results/PLFA_ABS_counts.xlsx', rowNames = T)
plfa_indices = read.xlsx('graftM_genes/Data/PLFA_indicators.xlsx')
rownames(plfa_indices) = plfa_indices$X1
qPCR_all = read.xlsx('graftM_genes/Results/qPCR_ABS_counts.xlsx', rowNames = T)
qPCR = read.csv('graftM_genes/Data/N and 16s and ITS.csv', row.names = 1)
env = read.csv('graftM_genes/Data/mag_env_for_shotgun_samples.csv')
rownames(env) = env$Hoosfield.ID

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
  
  df = data.frame(qPCR = qPCR[i], estimate = plfa[i])
  colnames(df) = c('qPCR', 'PLFA estimate')
  
  plot1 = ggplot(df, aes(x = qPCR, y = `PLFA estimate`)) +
    geom_point() +
    stat_poly_line() +
    stat_poly_eq() +
    ggtitle(i) +
    theme_classic()
  
  ggsave(paste0('graftM_genes/Figures/qPCR_vs_estimate/', i, '_PLFA.png'), plot1,
         height = 4, width = 5)
}
