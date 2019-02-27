library(ggplot2)
dir = 'C:/Liang/Googlebox/Research/Project3/simdata_1e+07newpara/1e+07'
jclabel = c(0,.2,.4,.6,.8,1)
plabel = c(0,.1,.3,1,3,10)
richnesstable=c()
for(phi_index in c(1:6)){
  for(psi_index in c(1:6)){
    filename = paste0('/Hevent',phi_index,psi_index,'.Rdata')
    fullfile = paste0(dir,filename)
    source(file = fullfile)
    speciesrichness = t(events)[,2]
    phi = plabel[phi_index]
    psi = jclabel[psi_index]
    tempt = cbind(speciesrichness,phi,psi)
    tempt = cbind(tempt,c(1:dim(events)[2]))
    richnesstable = rbind(richnesstable,tempt)     
  }
}

richdf = data.frame(richnesstable)
colnames(richdf) = c('Richness','phi','psi','Turnover')

richplot <- ggplot(data=richdf,aes(x=Turnover,y=Richness))+
  geom_line(linetype = "dashed")+
  geom_point()+
  facet_grid(phi~psi, labeller=label_both)+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))

richplot
