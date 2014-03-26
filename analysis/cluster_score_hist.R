
library(ggplot2)

data <- read.csv('../src/PLM_g_1_bossi_hpa_all.csv')

   
    fig <- ggplot(data, aes(x=PPI, y=ClusterScore, group)) +
            geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
