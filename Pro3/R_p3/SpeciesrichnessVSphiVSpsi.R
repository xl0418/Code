library(ggplot2)
library("lattice")

dir = 'C:/Liang/Code/Pro3/10_9simdata/'
richness = matrix(0,5,5)
# species richness of the whole test.
for(i in c(0:4)){
  for(j in c(0:4)){
    comb=paste0(i,j)
    Rname =  paste0(dir,'R',comb,'.Rdata')
    source(Rname)
    no_species = ncol(R)
    richness[i+1,j+1]=no_species
    }
}

## Example data
colnames(richness)=paste0(expression(phi), "=", c(1:5) )
rownames(richness)=paste0(expression(psi), "=", c(1:5) )

## Try it out
par(mar=c(3,4,2,2))
levelplot(t(richness[c(nrow(richness):1) , ]))

# #heatmap
# ggplot(richness_data, aes(phi_axis, psi_axis )) +
#   geom_tile(aes(fill = spec11), color = "white") +
#   scale_fill_gradient(low = "red", high = "green",limits = range(0,1)) +
#   ylab("psi") +
#   xlab("phi") +
#   theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 12),
#         plot.title = element_text(size=16),
#         axis.title=element_text(size=14,face="bold"),
#         axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(fill = "per capita probability")

