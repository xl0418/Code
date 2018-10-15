library(ape)
emdatadir = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/aze_extant_renamed.tre'

fulltree = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/aze_full_renamed.tre'

emdata = read.tree(fulltree)
plot(emdata,show.tip.label = FALSE)

phylo2L = function(emdata,dropextinct = TRUE){
brt=branching.times(emdata)
if(min(brt)<0){
  brt = brt+abs(min(brt))
}
brt_preL = c(brt[emdata$edge[,1]-length(emdata$tip.label)])
pre.Ltable = cbind(brt_preL,emdata$edge,emdata$edge.length,brt_preL-emdata$edge.length)
extantspecies.index = pre.Ltable[which(pre.Ltable[,5]<=0),3]
tipsindex = c(1:num.species)
extinct.index3 = subset(tipsindex,!(tipsindex %in% extantspecies.index))

eeindicator = matrix(0,length(emdata$edge.length),1)
eeindicator[match(extantspecies.index,pre.Ltable[,3])]=-1
ext.pos = match(extinct.index3,pre.Ltable[,3])
eeindicator[ext.pos]= pre.Ltable[ext.pos,5]
pre.Ltable=cbind(pre.Ltable,eeindicator)

sort.L = pre.Ltable[order(pre.Ltable[,1],decreasing = TRUE),]
num.species=emdata$Nnode+1
nodesindex = unique(emdata$edge[,1])
L = sort.L
realL=NULL
do=0
while(do == 0){
j = which.min(L[,3])
daughter = L[j,3]
# if(daughter %in% nodesindex){
#   do = 1
# }else{
parent = L[j,2]
if(parent %in% nodesindex){
L[which(L[,2]==parent),2] = daughter
if(length(which(L[,3]==parent))==0){
  realL = rbind(realL,L[j,],row.names = NULL)
  L = L[-j,,drop=FALSE]
}else{
L[which(L[,3]==parent),6] = L[j,6]
L[which(L[,3]==parent),3] = daughter
L = L[-j,,drop=FALSE]
}
}else{
  realL = rbind(realL,L[j,],row.names = NULL)
  L = L[-j,,drop=FALSE]
}

if(nrow(L)==0){
  do = 1
}
# }
}
realL = realL[order(realL[,1],decreasing = T),]
if(dropextinct == FALSE){
#get the exact branching times.
lastlength.index = which(realL[,6]==-1)[length(which(realL[,6]==-1))]
realL[,1]=realL[,1]+realL[lastlength.index,4]
}
L = realL[,c(1,2,3,6)]

daughter.index = L[,3]
daughter.realindex = c(1:nrow(L))
parent.index = L[,2]
parent.realindex = match(parent.index, daughter.index)

L[,2]=parent.realindex
L[,3]=daughter.realindex
L[1,2] = 0
L[1,3] = -1
L[2,2] = -1
for(i in c(2:nrow(L))){
  if(L[i-1,3]<0){
    mrows = which(L[,2]==abs(L[i-1,3]))
    L[mrows,2] = L[i-1,3]
    L[mrows,3] = -1* L[mrows,3]
  }
}
return(L)
}
L_ext = phylo2L(emdata,dropextinct = F)

phylo_test = DDD::L2phylo(L_ext,dropextinct = F)
plot(phylo_test,show.tip.label = FALSE)
