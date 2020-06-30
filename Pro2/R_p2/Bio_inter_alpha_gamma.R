library(ggthemes)
library(DDD)
library(ggplot2)


whale_data <- c(3.186932  , 3.14145   , 3.19434156, 3.23044892, 2.83016414,
                2.96988321, 3.14904428, 3.08974852, 3.03899361, 3.01836305,
                3.37965384, 3.24482357, 3.03590789, 3.06096064, 2.731966 )
mean(whale_data)
omega <- function(deltaz,alpha, whale){
  contri <- exp(-exp(-alpha*(deltaz-whale)^2))
  return(contri)
}
gamma_term <- function(deltaz,gamma){
  contri <- exp(1-gamma*(deltaz-3.05)^2)
  return(contri)
}



alpha.range = 99.31  #seq(0,1,0.01)
gamma.range = 4.84
deltaz.range = seq(2.5, 3.7,0.01)

omega_com <- NULL
for (species in c(1:length(whale_data))) {
  species_trait <- whale_data[species]
  omega.temp <-  sapply(deltaz.range, omega, alpha = alpha.range, whale = species_trait)
  omega_com <- rbind(omega_com, cbind(deltaz.range, omega.temp, species))
}
gamma.temp <- cbind(deltaz.range, gamma_term(deltaz = deltaz.range, gamma = gamma.range))


omega_df <- as.data.frame(omega_com)
names(omega_df) <- c('deltaz','value','species')

gamma_df <- as.data.frame(gamma.temp)
names(gamma_df) <- c('deltaz','value')


combined_plot <- ggplot(data=omega_df, aes(x=deltaz, y = value, group = species))+
  geom_line(aes(color = 'Competition niches'),lwd = 1)+
  ylab(bquote('Effect on the population growth ' ~ frac(N[t+1],N[t])))+
  geom_line(data = gamma_df,aes(y = value ,color ='Environmental niche'), lwd = 2)+
  theme_bw()+ xlab(bquote(Log[10] ~ mu))+
  theme(text = element_text(size=20),legend.position = c(0.8, 0.9), 
        legend.title = element_blank())+
  scale_x_continuous( breaks = seq(2.5,4,0.5))+
  scale_y_continuous( breaks = seq(0,2.5,0.5), labels = c('0', '50%', '100%', '150%', '200%', '250%'))+
  annotate("text", x = 3.1, y = 0.35, label = bquote(alpha==.(alpha.range)),size = 6)+
  annotate("text", x = 3.3, y = 1.5, label = bquote(gamma==.(gamma.range)),size = 6)

phi.plot.name = 'C:\\Liang\\Googlebox\\Research\\Project2\\BaleenWhales\\result_cluster\\results_ms_posterior\\biointer_alpha_gamma.png'
ggsave(phi.plot.name,width = 12,height = 8)

