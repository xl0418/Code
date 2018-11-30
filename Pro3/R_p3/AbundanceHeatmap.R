library(ggplot2)
library(ggthemes)
library(viridis)
Numspec = c()
for(i in c(1:25)){
  Numspec = c(Numspec, rnorm(1000,i*100,100))
}

psi = rep(c(0:4),each = 5000)
phi1 = rep(c(0:4),each = 1000)
phi = rep(phi1, 5)
phi = as.character(phi)
psi = as.character(psi)
abundid = rep(c(1:1000),25)
groupind = rep(c(1:25),each = 1000)
groupind = as.character(groupind)
htmaptest = data.frame(Numspec,psi,phi,abundid,groupind)
htmaptest$groupind =  factor(htmaptest$groupind, levels = c(1:25)) 


p <- ggplot(htmaptest, aes(x = abundid, y = groupind, fill = Numspec))

p_out <- p + geom_tile() +
  scale_fill_viridis_c(option = "A") +theme_hc()+
  labs(title = "Abundance distribution",
       x = "Abundance", y = "Group", fill = "Number of species") +
  theme(legend.position = "bottom", legend.box.just = "bottom")

p_out
