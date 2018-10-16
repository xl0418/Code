source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
emdatadir = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/aze_extant_renamed.tre'

fulltree = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/aze_full_renamed.tre'

emdata = read.tree(emdatadir)
plot(emdata,show.tip.label = FALSE)


dropextinct = F
L_ext = phylo2L(emdata)
phylo_test = DDD::L2phylo(L_ext,dropextinct = dropextinct)
plot(phylo_test,show.tip.label = FALSE)

