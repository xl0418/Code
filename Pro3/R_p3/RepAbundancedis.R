library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
source(paste0(getwd(),'/g_legend.R'))
source('C:/Liang/Code/Pro3/R_p3/event2L.R', echo=TRUE)
source('C:/Liang/Code/Pro3/R_p3/barplot3d.R', echo=TRUE)
moviedir = 'C:/Liang/Googlebox/Research/Project3/batchsim_results/'
dir = 'C:/Liang/Googlebox/Research/Project3/batchsim_results/1e+07/'
scenario = c('LR','MR','HR')
sce.short = c('L','M','H')
jclabel = c(0,0.5,1)
plabel = c(0,1e2,1e4,1e6,1e8,-1)
diversity.upperlimit = 250
diversity.lowerlimit = 80
lowcol = "#F2F0F7"
highcol = '#54278F'
tas1 = list()
abund.df = NULL
max.logabund = 12
probs = c(0,.25,.5,.75,1)
for(i_n in c(1:3)){
  sce = scenario[i_n]
  f.name = sce.short[i_n]
  for(i in c(1:3)){
    for(j in c(1:6)){
      abund=NULL
      for(rep in c(61:100)){
        comb=paste0(i,j)
        rname = paste0(dir,'spatialpara1e+07',f.name,comb,'/',sce,comb,'rep',rep,'.csv')
        Rs = read.csv(rname,header = FALSE)
        log.Rs = log10(Rs)
        freq = hist(as.numeric(log.Rs))
        counts = freq$counts
        abund = rbind(abund, counts)
      }
      quantile1 = apply(abund, MARGIN = 2, FUN = quantile, probs = probs)
      col.quan = dim(quantile1)[2]
      if(col.quan<12){
        quantile1 <- cbind(quantile1,matrix(0,5,12-col.quan))
      }
      abund.df = rbind(abund.df,cbind(quantile1,i,j,i_n,probs))
    }
  }
}

colnames(abund.df) <- c(c(1:12),c('i','j','scenario','prob'))
mean.df <- abund.df[which(abund.df[,16]==0.5),1:12]
mean.list <- list()
group.index=1
for(i in c(1:nrow(mean.df))){
  mean.list <- rbind(mean.list,cbind(c(1:12),mean.df[i,],group.index))
  group.index = group.index+1
}
mean.dateframe <- as.data.frame(mean.list)
colnames(mean.dateframe) <- c('x','freq','group')


ids <- factor(c('0.25','0.5','0.75'))

values <- data.frame(
  id = ids,
  value = c(0.25,0.5,0.75)
)

positions <- data.frame(
  id = rep(ids, each = 12),
  x = rep(c(1:12),3),
  y = c(quantile1[2,],quantile1[3,],quantile1[4,])
)



# Currently we need to manually merge the two together
datapoly <- merge(values, positions, by = c("id"))

ggplot(data = datapoly, aes(x = x, y = y)) +
  geom_area(data = subset(positions[positions$id==0.75,], y >= 0), fill = "blue", alpha = 0.2) +
  geom_area(data = subset(positions[positions$id==0.5,], y >= 0), fill = "blue", alpha = 0.6) +
  geom_area(data = subset(positions[positions$id==0.25,], y >= 0), fill = "white", alpha = 0.2) +
  theme_bw()


library(viridis)
d <- data.frame(x = rep(1:5, 3) + c(rep(0, 5), rep(0.3, 5), rep(0.6, 5)),
                y = c(rep(0, 5), rep(1, 5), rep(3, 5)),
                height = c(0, 1, 3, 4, 0, 1, 2, 3, 5, 4, 0, 5, 4, 4, 1))
ggplot(d, aes(x, y, height = height, group = y, fill = factor(x+y))) +
  geom_ridgeline_gradient() +
  scale_fill_viridis(discrete = TRUE, direction = -1, guide = "none")


ggplot(mean.dateframe, aes(x = x, y =group,fill = ..x..)) + 
  geom_density_ridges_gradient(scale = 2) + theme_ridges()+theme_wsj()+
  # facet_wrap(vars(psi),ncol = 5)+
  scale_fill_viridis(name = 'Richness', option = "C")+
  labs(title = 'Abundance distribution')+
  theme(legend.title = element_text(colour="black", size=10, face="bold"),
        legend.position = "bottom")

