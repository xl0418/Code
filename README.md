# Code
Under construction...  
This repository contains all the codes of my PhD projects.  
For the details of my PhD projects, please follow the links below: 

[Project 1: Detecting local diversity-dependence in diversification](https://github.com/xl0418/PhdIntroProject1)

[Project 2: Trait-population coevolution model](https://github.com/xl0418/PhdIntroProject2)

Project 3: under construction

So in this page, I am ganna post something interesting about coding that I developed in the projects. I am open for contributions but my code is a bit of mess that you may feel hard to read. There are also some hard codeings that the scripts are nested to each other. I will improve them slowly :-)  

## Plot cycles at the tips of a phylogenetic tree
In the first project, I was working on a two-location model. A phylogenetic tree that can imply where the species live is useful to model users. Thus I explored how to color the tree and how to plot cycles at the tips. Now this can be done by the function `plottree`

For example, we can generate a tree under our spatial model and save it.  
```
pars=c(0.8,0.2,20) #parameters of speciation rate, extinction rate and carrying capacity
result = sddsim(n=2,parsN=c(1,1),age=15,pars=pars , seed_fun = seed_fun, lambda_allo0 = 0.2, M0=1,K_fix = 1)
save(result,file = filename)
```  
Then, the function can be applied to the data and plot the tree like follows  
`
plottree(file = filename,dropextinct =T)
`
<div align=center><img width="550" height="550" src="https://github.com/xl0418/PhdIntroProject2/blob/master/Example/exampletree.png"/></div>  

The red and blue denote the emdemic species while the green denotes the widespread species. At tips, two cycles are implying two locations. Filled cycle means the species occupies that location. 

The ledgend's position can be adjusted but not from the arguments of the function. If you want to contribute, pls click [`plottree`](https://github.com/xl0418/Code/blob/0ebc0a244d3547d757503395f1cbedc9638b9261/Pro1/code_pro1/Plottree.R)

### phylo2L function 

I guess this function is specially useful to [our group](https://www.rug.nl/staff/r.s.etienne/) in which we play with L table. L table is an alternative way to a phylo class for phylogenetic information storage. A function called L2phylo has been implemented in the [DDD package](https://cran.r-project.org/web/packages/DDD/index.html) that converts an L table to a phylo class. This function [`phylo2L`](https://github.com/xl0418/Code/blob/99133e6e5744be7382c038edc5701cd494d8e76c/Pro2/R_p2/phylo2L.R) does the conversion the other way around. Thus, if you want to apply your model to an empirical data. This may be useful to you. 


### pruneL function

Actually, this toy function can be replaced by the previous funtion [`phylo2L`](https://github.com/xl0418/Code/blob/99133e6e5744be7382c038edc5701cd494d8e76c/Pro2/R_p2/phylo2L.R) combining with L2phylo when specifying drop extinction true. But this function [`pruneL`](https://github.com/xl0418/Code/blob/f4dfd4acc15af6855572fb4659f396cea14bb83b/Pro2/R_p2/pruneL.R) can prune an L table directly. Don't you think that will save you one second?  

### Animation of the trait-population coevolution model

### Animation of ABC-MCMC

### ABC-sequential Monte Calo algorithm
