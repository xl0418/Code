library(ggplot2)
library(ggthemes)
dir = 'C:/Liang/Code/Pro3/data1_20181220/'

neufile =  paste0(dir,'neuRl.Rdata')
source(file = neufile)
norR = R/sum(R)

norRdf = as.data.frame( cbind(c(1:length(norR)),sort(c(norR),decreasing = TRUE)))
colnames(norRdf) = c('Species rank','Percentage of species abundance')
ggplot(data=norRdf, aes(x=norRdf$`Species rank`, y=norRdf$`Percentage of species abundance`)) +
  geom_line(linetype = "dashed")+ geom_rangeframe()+
  geom_point()+theme_tufte()+xlab('Species rank')+ylab('Percentage of species abundance')+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
