library(ggplot2)
library(ggthemes)
library(viridis)
library(ggridges)
library(tidyverse)
source('C:/Liang/Code/Pro2/R_p2/seprate_vec.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/theme_henrik.R', echo=TRUE)

setwd("C:/Liang/PhdIntroProject2/smcplot")

tree = 2
fileunsort_name = paste0('smcdataa',tree,'.csv')
gammadata = read.csv(fileunsort_name)
gammadata_matrix = as.matrix(t(gammadata))


gamma_df = as.data.frame(as.table(gammadata_matrix))[,2:3]
colnames(gamma_df) = c('Iteration','Samples')
levels(gamma_df$Iteration) = c(1:30)

# Mountain plot of the evolution of SMC
ggplot(gamma_df, aes(x = Samples, y = Iteration)) + 
  geom_density_ridges_gradient(aes(fill = factor(as.integer(gamma_df$Iteration) %% 2)),
                               scale = 3,size =0.25,alpha = .4, color = "lightblue") +
  theme_ridges()+theme_henrik(grid='', legend.position='none')+
  scale_fill_manual(name = 'Richness',values = c('0' = '#2A7FFF', '1' = '#5599FF'))+
  labs(title = 'Evolution of SMC')+
  theme(legend.position='none')+
  scale_y_discrete(limits = rev(levels(gamma_df$Iteration)))


# Heatmap of SMC
finvec=c()
for(i in c(1:30)){
  startp = (i-1)*10000+1
  endp = i*10000
  focalvec = gamma_df$Samples[startp:endp]
  finvec1 = seprate_vec(vec = focalvec,length =2000,binwidth = 0.001)
  finvec = rbind(finvec, cbind(finvec1,i))
}

htdf = data.frame(finvec)
colnames(htdf) = c('Frequence','value','Iteration')
# htdf$Iteration = as.character(htdf$Iteration)

p <- ggplot(htdf, aes(x = value, y = ordered(Iteration, levels =rev(sort(unique(htdf$Iteration)))),
                      fill = Frequence))

p_out <- p + geom_tile() +
  scale_fill_viridis_c(option = "A") +theme_hc()+
  labs(title = "Evolution of SMC",
       x = expression(alpha), y = "Iteration", fill = "Number of samples") +
  theme(legend.position = "bottom", legend.box.just = "bottom")


p_out




