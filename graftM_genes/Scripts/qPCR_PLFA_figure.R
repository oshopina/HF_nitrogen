library(openxlsx)
library(ggplot2)
library(ggpmisc)
library(dplyr)

plfa = read.xlsx('graftM_genes/Results/PLFA_ABS_counts.xlsx', rowNames = T)
qPCR = read.csv('graftM_genes/Data/N and 16s and ITS.csv', row.names = 1)
env = read.csv('graftM_genes/Data/mag_env_no_outliers.csv')
rownames(env) = env$Hoosfield.ID

env_sorted = env[order(env$pH),]
qPCR_sorted = qPCR[rownames(env_sorted),] %>% na.omit()
plfa_sorted = plfa[rownames(env_sorted),] %>% na.omit()

mypal.pH <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"))

############################### Genes by qPCR ##################################
temp = env_sorted[rownames(qPCR_sorted),]

amoa = ggplot(qPCR_sorted[rownames(temp),], aes(y = amoAbaccopiesgdrysoil, x = temp$pH, color = temp$pH)) +
  geom_point(size = 2.5) +
  scale_color_gradientn(colours = mypal.pH(256)) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('pH') +
  ylab('amoA') +
  geom_smooth()

napa = ggplot(qPCR_sorted[rownames(temp),], aes(y = NapAqtygdrysoil, x = temp$pH, color = temp$pH)) +
  geom_point(size = 2.5) +
  scale_color_gradientn(colours = mypal.pH(256)) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('pH') +
  ylab('napA') + 
  geom_smooth()

norb = ggplot(qPCR_sorted[rownames(temp),], aes(y = cNorBgdrysoil, x = temp$pH, color = temp$pH)) +
  geom_point(size = 2.5) +
  scale_color_gradientn(colours = mypal.pH(256)) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('pH') +
  ylab('norB') +
  geom_smooth()

temp1 = env_sorted[rownames(plfa_sorted),]

amoa_plfa = ggplot(plfa_sorted, aes(y = amoA, x = temp1$pH, color = temp1$pH)) +
  geom_point(size = 2.5) +
  scale_color_gradientn(colours = mypal.pH(256)) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('pH') +
  ylab('amoA') +
  geom_smooth()

napa_plfa = ggplot(plfa_sorted, aes(y = napA, x = temp1$pH, color = temp1$pH)) +
  geom_point(size = 2.5) +
  scale_color_gradientn(colours = mypal.pH(256)) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('pH') +
  ylab('napA') +
  geom_smooth()

norb_plfa = ggplot(plfa_sorted, aes(y = norB, x = temp1$pH, color = temp1$pH)) +
  geom_point(size = 2.5) +
  scale_color_gradientn(colours = mypal.pH(256)) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('pH') +
  ylab('norB') +
  geom_smooth()


###############################################################################

library(patchwork)

combined_plot = amoa + amoa_plfa + napa + napa_plfa + norb + norb_plfa + plot_layout(ncol = 2, nrow = 3)

ggsave('graftM_genes/Figures/qPCR_combined_plot.svg', combined_plot, 
       height = 10, width = 10)
