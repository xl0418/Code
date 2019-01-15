
# Option 2: density plot
for(i in c(1:25)){
  Richness = c(Richness, rnorm(1000,i*100,100))
}

# Richness = rnorm(25000,300,400)
psi = rep(c(0:4),each = 5000)
phi1 = rep(c(0:4),each = 1000)
phi = rep(phi1, 5)
phi = as.character(phi)
psi = as.character(psi)
colorvalue1 = rep(c(0,1,0,1,0),each = 1000)
colorvalue = rep(colorvalue1,5)
colorvalue = as.character(colorvalue)
dftest = data.frame(Richness,psi,phi,colorvalue)

ggplot(df, aes(x = Richness, y = phi,fill = ..x..)) + 
  geom_density_ridges_gradient(scale = 2) + theme_ridges()+theme_wsj()+
  facet_wrap(vars(psi),ncol = 6)+
  scale_fill_viridis(name = 'Richness', option = "C")+
  labs(title = 'Abundance distribution')+
  theme(legend.title = element_text(colour="black", size=10, face="bold"),
        legend.position = "bottom")


# Option 2-1: density plot
Richness = c()
for(i in c(1:25)){
  Richness = c(Richness, rnorm(1000,i*100,100))
}

psi = rep(c(0:4),each = 5000)
phi1 = rep(c(0:4),each = 1000)
phi = rep(phi1, 5)
phi = as.character(phi)
psi = as.character(psi)
colorvalue1 = rep(c(0,1,0,1,0),each = 1000)
colorvalue = rep(colorvalue1,5)
colorvalue = as.character(colorvalue)
groupind = rep(c(1:25),each = 1000)
groupind = as.character(groupind)
dftest = data.frame(Richness,psi,phi,colorvalue,groupind)
dftest$groupind =  factor(dftest$groupind, levels = c(1:25)) 

ggplot(dftest, aes(x = Richness, y =groupind,fill = ..x..)) + 
  geom_density_ridges_gradient(scale = 2) + theme_ridges()+theme_wsj()+
  # facet_wrap(vars(psi),ncol = 5)+
  scale_fill_viridis(name = 'Richness', option = "C")+
  labs(title = 'Abundance distribution')+
  theme(legend.title = element_text(colour="black", size=10, face="bold"),
        legend.position = "bottom")


# Option 3
ggplot(df, aes(x = Richness, y = phi)) + 
  geom_density_ridges_gradient(aes(fill = factor(as.integer(df$colorvalue) %% 2)),
                               scale = 2,size =0.25,alpha = .4, color = "white") + theme_ridges()+theme_wsj()+
  facet_wrap(vars(psi),ncol = 5)+
  scale_fill_manual(name = 'Richness',values = c('0' = '#2A7FFF', '1' = '#5599FF'))+
  labs(title = 'Abundance distribution',x = "psi",
       y = "phi")+
  theme(legend.position='none')

# Option 3-1
ggplot(dftest, aes(x = Richness, y = groupind)) + 
  geom_density_ridges_gradient(aes(fill = factor(as.integer(dftest$colorvalue) %% 2)),
                               scale = 2,size =0.25,alpha = .4, color = "white") + theme_ridges()+theme_wsj()+
  scale_fill_manual(name = 'Richness',values = c('0' = '#2A7FFF', '1' = '#5599FF'))+
  labs(title = 'Abundance distribution',x = "psi",
       y = "phi")+
  theme(legend.position='none')