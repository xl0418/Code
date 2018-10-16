source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
emdatadir = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/aze_extant_renamed.tre'

fulltree = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/aze_full_renamed.tre'

emdata = read.tree(fulltree)
plot(emdata,show.tip.label = FALSE)
brt=branching.times(emdata)
if(min(brt)<0){
  brt = brt+abs(min(brt))
}
range(brt)


dropextinct = F
L_ext = phylo2L(emdata)
brt_preL = c(brt[emdata$edge[,1]-length(emdata$tip.label)])
range(brt_preL)


phylo_test = DDD::L2phylo(L_ext,dropextinct = dropextinct)
plot(phylo_test,show.tip.label = FALSE)
brt=branching.times(emdata)
if(min(brt)<0){
  brt = brt+abs(min(brt))
}
range(brt)