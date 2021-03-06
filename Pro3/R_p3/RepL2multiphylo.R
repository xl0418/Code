library(DDD)
library(ape)

dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces/'
sce.short = c('H','M','L')
scenario = NULL
sce.short.comb.vec = NULL
for(i.letter in sce.short){
  for(j.letter in sce.short){
    sce.folder = paste0('sce',i.letter,j.letter)
    scenario = c(scenario,sce.folder)
    sce.short.comb = paste0(i.letter,j.letter)
    sce.short.comb.vec = c(sce.short.comb.vec,sce.short.comb)
  }
}
ev2l <- function(ltable, cond){
  ancestor0 <- cond[1:2]
  ancestor1 <- cond[3:4]
  
  row1 <- c(10^7, 0.0, 1.0, ancestor0[2])
  row2 <- c(10^7, 1, 2, ancestor1[2])
  ltable[,2] <- ltable[,2] + 1
  ltable[,3] <- ltable[,3] + 1
  ltable <- rbind(row1, row2, ltable)
  phylo_L <- DDD::L2phylo(ltable, dropextinct = TRUE)
  return(phylo_L)
}

for(i_n in c(1:9)){
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]
  for(i in c(1:5)){
    for(j in c(1:6)){
      print(paste(scefolder,i,j))
      multitree=list()
      comb=paste0(i,j)
      
      for(rep in c(1:100)){

        rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'psi',i,'s_phi',j,'rep',rep,'.csv')
        L.table = read.csv(rname,header = FALSE)
        L.table <- as.matrix(L.table)
        ancestor_con <- L.table[1,]
        L.table <- L.table[2:nrow(L.table),]
        # print(nrow(L.table))
        single.phylo = ev2l(L.table, ancestor_con)
        multitree = c(multitree,list(single.phylo))
      }
      multitreefile <- paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/','multitreen',letter.comb,comb,'.tre')

      class(multitree) <- "multiPhylo"
      write.tree(multitree,file=multitreefile)
    }
  }
}

