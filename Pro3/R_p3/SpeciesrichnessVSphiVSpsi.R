library(ggplot2)
library("lattice")

dir = 'C:/Liang/Code/Pro3/data1_20181220/'
richness = matrix(0,6,6)
# species richness of the whole test.
for(i in c(1:6)){
  for(j in c(1:6)){
    comb=paste0(i,j)
    Rname =  paste0(dir,'HR',comb,'.Rdata')
    source(Rname)
    no_species = ncol(R)
    richness[i,j]=no_species
    }
}

## Example data
colnames(richness)=paste0(expression(phi), "=", c(1:6) )
rownames(richness)=paste0(expression(psi), "=", c(1:6) )

## Try it out
par(mar=c(3,4,2,2))
levelplot(t(richness[c(nrow(richness):1) , ]),
          xlab="Phylogenetic strength", ylab="JC strength", main="Species richness for mediate spatial parameter")

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

