# Code
Under construction...  
This repository contains all the codes of my PhD projects.
For the details of my PhD projects, please follow the links below: 

[Project 1: Detecting local diversity-dependence in diversification](https://github.com/xl0418/PhdIntroProject1)

[Project 2: Trait-population coevolution model](https://github.com/xl0418/PhdIntroProject2)

Project 3: under construction

So in this page, I am ganna post something interesting about coding that I developed in the projects. 

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
