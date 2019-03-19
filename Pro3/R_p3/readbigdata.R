fil = 'c:\\Liang\\Googlebox\\Research\\Project3\\simdata_190122\\10e6\\HD11.csv'
x = read.csv(fil,header = FALSE)
xm <- as.matrix(x)
hist(xm)
