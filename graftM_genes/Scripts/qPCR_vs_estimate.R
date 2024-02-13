library(openxlsx)
library(ggplot2)
library(ggpmisc)

plfa = read.xlsx('graftM_genes/Results/PLFA_ABS_counts.xlsx', rowNames = T)
qPCR_all = read.xlsx('graftM_genes/Results/qPCR_ABS_counts.xlsx', rowNames = T)
qPCR = read.csv('graftM_genes/Data/N and 16s and ITS.csv', row.names = 1)
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
