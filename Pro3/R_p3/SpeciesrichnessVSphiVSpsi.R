richness_data = data.frame()

for(i in phi_group){
  for(j in psi_group){
    richness_each = NULL
    for (seed in 1:100) {
      file = paste("C:/Liang/Googlebox/Research/Project3/simdata/",i,"phi",j,"psi","sim",seed,".Rdata",sep = "")
      load(file = file)
      richness_each = c(richness_each,result$events$ns[length(result$events$ns)])
      
    }
    rich_each_mean = mean(richness_each)
    rich_each_max = max(richness_each)
    rich_each_min = min(richness_each)
    richness_data0 = c(i,j,rich_each_mean,rich_each_max,rich_each_min)
    richness_data = rbind(richness_data,richness_data0)
    }
}
names(richness_data)=c("phi","psi","mean","max","min")
phi_axis = rep(c(0,0.001,0.01,0.1,1),5)
psi_axis =  rep(c(0,0.001,0.01,0.1,1),each = 5)


ggplot(data=richness_data, aes(x=phi, y=mean, group=1)) +
  geom_line(linetype = "dashed")+
  geom_point()

## example data (on my env, plotrgl(some_graph) crushes Rstudio, so I used volcano)

with(richness_data, scatter3D(x = phi_axis, y = psi_axis, z = mean, ticktype="detailed", pch=16, 
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

