library(POJCE)
library(plot3D)
library(rgl)
library(plot3Drgl)
library(DDD)
library(ggplot2)
source('C:/Liang/Googlebox/Research/Project3/Project3_R/event2L.R', echo=TRUE)

library(car)


phi = 0.1
psi = phi
sA = 0.5
sD = 0.5
ticks = 20000
seed = 29
L_size = 25
v = 0.001

result = pojce(phi=phi, psi=psi, sA=sA, sD=sD, ticks=ticks, randomSeed=5, L=L_size, v=v)
L = event2L(result = result)
phy = DDD::L2phylo(L,dropextinct =F)

plot(phy)
print(result$events)
print(length(which(L[,4]==-1)) == result$events$ns[length(result$events$ns)])

snapshot = 5
pop_vector = result$log$R[[snapshot]]
phy_dis_vector = result$log$D[[snapshot]]

#per capita probability of colonization
col_per_cap = function(phi,psi,pop_vector,D){
  num_spe = length(pop_vector)
  pop_matrix = matrix(nrow = num_spe,ncol = num_spe)
  for(i in 1:num_spe){
    pop_matrix[i,]=pop_vector/pop_vector[i]
  }
  p = pop_matrix^phi * D^psi
  per = rowSums(p/sum(p))
  return(per)
}

per_cap_vec = col_per_cap(phi = 1, psi = 1, pop_vector = pop_vector,D = phy_dis_vector)


#per species probability of colonization

col_per_spe = function(phi,psi,pop_vector,D){
  per_cap_vec = col_per_cap(phi = phi, psi = psi, pop_vector = pop_vector,D = D)
  per_vec = pop_vector*per_cap_vec
  per_spe_vec = per_vec/sum(per_vec)
  return(per_spe_vec)
}

per_spe_vec = col_per_spe(phi = 1, psi = 1, pop_vector = pop_vector,D = phy_dis_vector)


num_spe = length(pop_vector)

per_cap_col_dataframe = data.frame()


phi_range = seq(0,4,0.04)
psi_range = seq(0,4,0.04)
for(i in phi_range){
  for(j in psi_range){
    per_cap_vec = col_per_cap(phi = i, psi = j, pop_vector = pop_vector,D = phy_dis_vector)
    
    per_cap_col_dataframe = rbind(per_cap_col_dataframe,per_cap_vec)
  }
}
colnames(per_cap_col_dataframe) = paste0("spec",c(1:num_spe))
psi_axis = rep(psi_range,length(phi_range))
phi_axis = rep(phi_range,each = length(psi_range))
per_cap_col_dataframe = cbind(phi_axis,psi_axis,per_cap_col_dataframe)




## example data (on my env, plotrgl(some_graph) crushes Rstudio, so I used volcano)
volc <- reshape2::melt(volcano)  

with(per_cap_col_dataframe, scatter3D(x = phi_axis, y = psi_axis, z = spec12, ticktype="detailed", pch=16, 
                     xlab="phi", ylab="psi", zlab="per capita probability", main=""))
plotrgl(lighting = TRUE, smooth = TRUE, cex=2)
scatter3d(x = phi_axis, y = psi_axis, z = per_cap_col_dataframe$spec12,point.col = "blue", surface=FALSE)

## When graph is made, the left panel is focused
title3d(main = "Earthquakes off Fiji", line=4, cex=2)

## next3d() changes the focus into the right panel
next3d(clear = F)
title3d("Richter", cex = 2, line = 2)
title3d("Magnitude", cex = 2, line = 0.8)
# text3d(0, 1, 1.2, "Richter", cex=2)    # almost same
# text3d(0, 1, 1.1, "Magnitude", cex=2)

next3d(clear = F) # if you want to refocus the left panel


#heatmap
ggplot(per_cap_col_dataframe, aes(phi_axis, psi_axis )) +
  geom_tile(aes(fill = spec11), color = "white") +
  scale_fill_gradient(low = "red", high = "green",limits = range(0,1)) +
  ylab("psi") +
  xlab("phi") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "per capita probability")

