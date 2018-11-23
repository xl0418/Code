library(DDD)
library(ggplot2)
library(ggridges)

source('C:/Liang/Code/Pro3/R_p3/event2L.R', echo=TRUE)
source('C:/Liang/Code/Pro3/R_p3/multipleplot.R', echo=TRUE)
dir = 'C:/Liang/Code/Pro3/10_9simdata/'

r_df = NULL
for(i in c(0:4)){
  for(j in c(0:4)){
    comb=paste0(i,j)
    rname = paste0(dir,'R',comb,'.Rdata')
    source(rname)
    r_df_g = cbind(t(R),i,j)
    r_df = rbind(r_df,r_df_g)
    
  }
}
df = as.data.frame(r_df)
colnames(df) = c('Richness','psi','phi')
df$psi = as.character(df$psi)
df$phi = as.character(df$phi)

# useful
# ggplot(df, aes(x=Richness,y=psi, fill=phi)) + geom_density_ridges() +
#   theme_ridges() + 
#   theme(legend.position = "none")

dislist = list()
for(count in c(1:5)){
  dfpart = df[df['psi'] == count-1,]
  
  dislist[[count]] = ggplot(dfpart, aes(x=Richness, fill=phi)) + 
    geom_density(alpha=0.55)+
    theme_classic()+
    theme(legend.position = "None")
  # ggplot(dfpart, aes(x=Richness, fill=phi)) +
  #   geom_histogram(binwidth=.5, alpha=.5, position="identity")
}
multiplot(plotlist = dislist, cols = 1)
