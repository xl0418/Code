library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
library(ggtree)
library(ape)
source(paste0(getwd(),'/g_legend.R'))
source('C:/Liang/Code/Pro3/R_p3/multi3dbar.R', echo=TRUE)


# Load phylogenetic tree; Input 1
tree <- read.tree('C:/Liang/Googlebox/Research/Project3/BCI_data/control.tre')
# plot(tree)
# Read the BCI data
mydatapath = 'C:/Liang/Googlebox/Research/Project3/BCI_data/'
if (!exists('bci.full7')) {
  attach(paste(mydatapath, 'bci.full7.rdata', sep = ''))
}

# Extract the species mnemonic and the geographic info; Input 2
bci_alive_data <-
  subset(bci.full7,
         status = 'A',
         select = c('sp', 'gx', 'gy'))
sp_code <- unique(bci_alive_data[, 1])

# Get full species names
load('C:/Liang/Googlebox/Research/Project3/BCI_data/bci.spptable.Rdata')

sp_name <- NULL
genus_name <- NULL
for (i in c(1:length(sp_code))) {
  index <- which(bci.spptable$sp == sp_code[i])
  if (length(index) == 0) {
    sp_name <-
      c(sp_name, 0)
    genus_name <-
      c(genus_name, 0)
  } else {
    sp_name <-
      c(sp_name, bci.spptable$Species[index])
    genus_name <-
      c(genus_name, bci.spptable$Genus[index])
  }
}

species_name <- paste0(genus_name, '_', sp_name)
sp_code_name <- cbind(sp_code, species_name)

recognized_species <- which(species_name %in% tree$tip.label)
recognized_code_name <- sp_code_name[recognized_species,]

# Calculate the abundance of each species; Input 3
abundance <- NULL
for(i in c(1:nrow(recognized_code_name))) {
  abundance_species <- length(which(bci_alive_data[, 1] == recognized_code_name[i, 1]))
  abundance <- c(abundance, abundance_species)
}

# Wrap the sp code, species names and the abundance data into a data frame.
abundance_data <- data.frame(sp = recognized_code_name[,1], species = recognized_code_name[,2], abundance = abundance)

# Subsetting the species; Input 4
# not_in_species <- subset(species_name, !(species_name %in% tree$tip.label))
species_in <- subset(species_name, (species_name %in% tree$tip.label))

tree <- keep.tip(tree, species_in)
distance_matrix <- ape::cophenetic.phylo(tree)



interleave <- function(x,y){
  lx <- length(x)
  ly <- length(y)
  n <- max(lx,ly)
  as.vector(rbind(rep(x, length.out=n), rep(y, length.out=n)))
}

x_label_fontsize = 12
y_label_fontsize = 12
mu_title_fontsize = 16
mig_title_fontsize = 16
x_title_fontsize = 16
y_title_fontsize = 16


count1 = 1
p = list()
# AWMIPD mode
pdmode = 'exp'
a = 1

# BCI tree
p[[count1]] <- ggtree(tree,layout = "circular") #+xlim(0,age)

count1 = count1+1

# plot SAR
species.area = NULL
for(local.scale in c(1:19,seq(20,500,10))){
  mean.data = NULL

  ns.vec=NULL

  submatrix.vec = c(0:(500 %/% local.scale - 1))*local.scale
  for(row.num in submatrix.vec){
    local_row_x_half1 <- which(bci_alive_data$gx < (row.num+local.scale) & 
                          bci_alive_data$gx > row.num)
    local_row_x_half2 <- which(bci_alive_data$gx < (row.num+local.scale+500) & 
                                 bci_alive_data$gx > row.num+500)
    local_row_x <- c(local_row_x_half1, local_row_x_half2)
    local_row <- which(bci_alive_data$gy[local_row_x] < (row.num+local.scale) & 
                           bci_alive_data$gy[local_row_x] > row.num)
    
    local.richness = length(unique(bci_alive_data$sp[local_row]))
    ns.vec = c(ns.vec, local.richness)
  }
  quantile1 = quantile(ns.vec)
  species.area = rbind(species.area,c(quantile1,local.scale))
  
}

species.area.df <- as.data.frame(species.area)
colnames(species.area.df) <- c('0', '25', '50', '75', '100', 'area')
area = species.area.df$area
df_lineage_sar = species.area.df

df_lineage_sar$`0` = log(df_lineage_sar$`0`)
df_lineage_sar$`25` = log(df_lineage_sar$`25`)
df_lineage_sar$`50` = log(df_lineage_sar$`50`)
df_lineage_sar$`75` = log(df_lineage_sar$`75`)
df_lineage_sar$`100` = log(df_lineage_sar$`100`)
df_lineage_sar$area = log(df_lineage_sar$area)
df_min_max_sar = data.frame(id = "min_max", value = 1, x = c(df_lineage_sar$area,rev(df_lineage_sar$area)), y = c(df_lineage_sar$'0',rev(df_lineage_sar$'100')))
df_0025_sar = data.frame(id = "0025", value = 2, x = c(df_lineage_sar$area,rev(df_lineage_sar$area)), y = c(df_lineage_sar$'25',rev(df_lineage_sar$'75')))
df_mean_sar = data.frame(id = 'mean',y = df_lineage_sar$`50`,x = df_lineage_sar$area)
df_lineage_all_sar = rbind(df_min_max_sar,df_0025_sar)


p[[count1]] <- ggplot(df_mean_sar, aes(x = x, y = y)) +
  theme(legend.position="none",
        panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")
        ,plot.margin=unit(c(0,0,0,0),"cm"))+
  geom_polygon(data = df_min_max_sar, aes(  group = id),fill = "gray70", alpha = 0.8)+
  geom_polygon(data = df_0025_sar, aes( group = id),fill = "gray27", alpha = 0.8)+
  geom_line()+
  coord_cartesian(xlim=c(log(1),log(500)),ylim=c(0,log(1000))) + 
  scale_y_continuous(name = "No. of species",breaks = log(c(1,10,100,500)),labels = c('1','10','100', '500'))+
  scale_x_continuous(name = "Area",breaks = log(c(1,10,100, 500)),labels = c('1','10','100', '500'))+
  theme(axis.text.y=element_text(angle=90,size = y_label_fontsize),
        axis.text.x=element_text(size = x_label_fontsize))

count1 = count1+1

# SAD
x.breaks = seq(0,17,1)

Rs = abundance_data$abundance
log.Rs = log2(Rs)
freq = hist(as.numeric(log.Rs),plot=FALSE,breaks = x.breaks)
counts = freq$counts

mean.sim = counts

abund.df = cbind(mean.sim,c(1:(length(x.breaks)-1)))
colnames(abund.df) <- c('mean','species')
abund.df <- as.data.frame(abund.df)

my_labs <- interleave(seq(1,length(x.breaks)-1,2), "")
my_labs = my_labs[1:17]


p[[count1]]  <- ggplot(abund.df) +
  geom_bar( aes(x=species, y=mean),width = 0.6, stat="identity", fill="red", alpha=0.7) +
  geom_line(aes(species, mean),size=0.8,color="blue")+
  #theme_gdocs()+ #scale_color_calc()+
  scale_x_continuous(name="Abundance (log2)", breaks=seq(1,length(x.breaks)-1,1),labels = my_labs) +
  scale_y_continuous(name="Frequency",breaks=seq(0,60,20))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(),strip.text.x = element_text(size = 12, colour = "black"),
        strip.text.y = element_text(size = 12, colour = "black"))

count1 = count1+1

# AWMPID
species_name_order <- rownames(distance_matrix)
abundance_index <-
  match(species_name_order, abundance_data$species)
abundance <- abundance_data$abundance[abundance_index]
sp_order <- abundance_data$sp[abundance_index]

# phylogenetic distance
if (pdmode == 'inv') {
  ID = 1 / distance_matrix
  diag(ID) = 1
  AD.matrix = sweep(ID, MARGIN = 2, as.matrix(abundance), `*`)
} else if (pdmode == 'exp') {
  ID = exp(- a * distance_matrix)
  AD.matrix = sweep(ID, MARGIN = 2, as.matrix(abundance), `*`)
}

total.dvalues = rowSums(AD.matrix) * as.matrix(abundance)
awmipd_real_value = total.dvalues/sum(total.dvalues)
D.normalized = (total.dvalues - min(total.dvalues)) / (max(total.dvalues) -
                                                         min(total.dvalues))

geographic_dis <- bci_alive_data
order_in_geo_data <- match(geographic_dis$sp , sp_order)
geographic_dis['awmipd'] <- D.normalized[order_in_geo_data]

# low blue, high red
p[[count1]] <-
  ggplot2::ggplot(geographic_dis, aes(gx, gy)) + geom_point(aes(color = awmipd), size = 0.5, alpha = 0.8) +
  theme(
    legend.position = '',
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    panel.background = element_blank()
  ) +
  xlab("") + ylab("") + scale_color_gradient2(
    low = "#005CAF",
    mid = 'green',
    high = "#D0104C",
    midpoint = 0.5,
    name = "AWMIPD"
  )


count1 = count1+1

# LTT plot
age = 1
ltt_data <- ape::ltt.plot.coords(tree, tol = 1e-5)
ltt_data[,2] <- log(ltt_data[,2])
ltt_df <- data.frame(ltt_data)
colnames(ltt_df) <- c('time', 'species')

p[[count1]] <- ggplot(ltt_df, aes(x = time, y = species)) + geom_line() +
  theme(legend.position="none",
        panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")
        ,plot.margin=unit(c(0,0,0,0),"cm"))+
  coord_cartesian(xlim=c(-age,0),ylim=c(log(2),log(400))) + scale_y_continuous(name="No. of species",breaks = c(log(2),log(10),log(50),log(400)),labels = c(2,10,50,400))+
  scale_x_continuous(name="Time",breaks = -rev(seq(0,1,1/5)),labels = c('1','0.8','0.6','0.4','0.2','0'))+
  theme(axis.text.y=element_text(angle=90,size = y_label_fontsize),axis.text.x=element_text(size = x_label_fontsize))


  

m = matrix(1:5,ncol = 1)


col_labels = c('BCI data')
row_labels = c('Tree', 'SAR','SAD','SPD','LTT')

phi1 <- textGrob(row_labels[1], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi2 <- textGrob(row_labels[2], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi3 <- textGrob(row_labels[3], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi4 <- textGrob(row_labels[4], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi5 <- textGrob(row_labels[5], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))

psi1 <- textGrob(col_labels[1], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))


row_titles <- arrangeGrob(phi1,phi2,phi3,phi4,phi5,ncol = 1)
column_titles <- arrangeGrob(psi1,ncol = 1)


# label = textGrob("Number of lineages",gp=gpar(fontsize=y_title_fontsize), rot = 90)

g_ltt4 = arrangeGrob(grobs = p, layout_matrix = m)
g_ltt5 = textGrob("area", gp=gpar(fontsize=x_title_fontsize),rot = 0)
g_ltt1 = textGrob("")
g_a = textGrob("(a)")
g_b = textGrob("(b)")

ltt.sce <- grid.arrange(column_titles,g_ltt1,g_ltt4,row_titles,ncol = 2,widths = c(20,1),heights = c(2,25))

grid.arrange(column_titles, p[[1]], p[[2]], p[[3]], p[[5]], p[[4]], ncol = 2, 
             layout_matrix = rbind(c(1,1), c(2,3), c(4, 5), c(6,6)),
             heights = c(1,10,10,10))

dir_save <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces_results/'
savefilename <- paste0(dir_save,'mixing_result_BCI.pdf')
ggsave(savefilename,ltt.sce,width = 3,height = 10)




