library(rgl)
con_probability_per_capita <- function(N,D,phi,psi){
  if(length(N) == dim(D)[1]){
    dim.N = length(N)
    N.matrix = matrix(0,dim.N,dim.N)
    pi = c()
    for(i in c(1:dim.N)){
      N.matrix[i,]=1-N[i]/sum(N)
      pi[i]=sum(N.matrix[i,]^psi *(D[i,]/sum(D[i,]))^phi)
    }
    pi = pi/sum(pi)
    return(pi)
  }else{
    print('The length of N should equal the size of D')
  }
}

# phi vs. D
p.p.d <- function(phi,Dvec){
  psi <- 0.5
  N <- c(50,2,8,40)
  D <- matrix(c(0,5,10,20,5,0,5,10,10,5,0,5,20,10,5,0),4,4)
  D[1,4] = Dvec
  D[lower.tri(D)] <- t(D)[lower.tri(D)]
  pi <- con_probability_per_capita(N=N,D=D,phi=phi,psi=psi)
  return(pi[4])
}

phi = seq(0,3,length.out = 100)
Dvec = seq(0,1000,length.out = 100)
p.p.d.vec <- Vectorize(p.p.d)

z2 <- outer(phi, Dvec, p.n.d.vec)
nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(z2, nbcol)
persp3d(phi, Dvec, z2, theta=50, phi=25, expand=0.75, col=color[zcol],
        ticktype="detailed", xlab="phi", ylab="Distance", zlab="Probability",axes=TRUE)



# phi vs. psi
p.p.p <- function(phi,psi){
  N <- c(50,2,8,40)
  D <- matrix(c(0,5,10,20,5,0,5,10,10,5,0,5,20,10,5,0),4,4)
  D[lower.tri(D)] <- t(D)[lower.tri(D)]
  pi <- con_probability_per_capita(N=N,D=D,phi=phi,psi=psi)
  return(pi[4])
}

phi = seq(0,3,length.out = 100)
psi = phi
p.p.p.vec <- Vectorize(p.p.p)
z2 <- outer(phi, psi, p.p.p.vec)
nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(z2, nbcol)
persp3d(phi, psi, z2, theta=50, phi=25, expand=0.75, col=color[zcol],
        ticktype="detailed", xlab="phi", ylab="psi", zlab="Probability",axes=TRUE)




# phi vs. N
p.p.N <- function(phi,Nvec){
  N <- c(50,2,8,40)
  N[1] <- Nvec
  psi <- 0.4
  D <- matrix(c(0,5,10,20,5,0,5,10,10,5,0,5,20,10,5,0),4,4)
  D[lower.tri(D)] <- t(D)[lower.tri(D)]
  pi <- con_probability_per_capita(N=N,D=D,phi=phi,psi=psi)
  return(pi[1])
}

phi = seq(0,3,length.out = 100)
Nvec = seq(0,1000,length.out = 100)
p.p.N.vec <- Vectorize(p.p.N)
z2 <- outer(phi, Nvec, p.p.N.vec)
nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(z2, nbcol)
persp3d(phi, Nvec, z2, theta=50, phi=25, expand=0.75, col=color[zcol],
        ticktype="detailed", xlab="phi", ylab="N", zlab="Probability",axes=TRUE)





# N vs. psi
p.N.p <- function(psi,Nvec){
  N <- c(50,2,8,40)
  N[1] <- Nvec
  phi <- 0.4
  D <- matrix(c(0,5,10,20,5,0,5,10,10,5,0,5,20,10,5,0),4,4)
  D[lower.tri(D)] <- t(D)[lower.tri(D)]
  pi <- con_probability_per_capita(N=N,D=D,phi=phi,psi=psi)
  return(pi[1])
}

psi = seq(0,3,length.out = 100)
Nvec = seq(0,1000,length.out = 100)
p.N.p.vec <- Vectorize(p.N.p)
z2 <- outer( psi,Nvec, p.N.p.vec)
nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(z2, nbcol)
persp3d(psi, Nvec, z2, theta=50, phi=25, expand=0.75, col=color[zcol],
        ticktype="detailed", xlab="psi", ylab="N", zlab="Probability",axes=TRUE)



# D vs. psi
p.D.p <- function(psi,Dvec){
  phi <- 0.4
  N <- c(50,2,8,40)
  D <- matrix(c(0,5,10,20,5,0,5,10,10,5,0,5,20,10,5,0),4,4)
  D[1,4] = Dvec
  D[lower.tri(D)] <- t(D)[lower.tri(D)]
  pi <- con_probability_per_capita(N=N,D=D,phi=phi,psi=psi)
  return(pi[4])
}

psi = seq(0,3,length.out = 100)
Dvec = seq(0,1000,length.out = 100)
p.D.p.vec <- Vectorize(p.D.p)

z2 <- outer(psi, Dvec, p.D.p.vec)
nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(z2, nbcol)
persp3d(psi, Dvec, z2, theta=50, phi=25, expand=0.75, col=color[zcol],
        ticktype="detailed", xlab="psi", ylab="Distance", zlab="Probability",axes=TRUE)





# D vs. N
p.D.N <- function(Dvec,Nvec){
  phi <- psi <- 0.4
  N <- c(50,2,8,40)
  N[4] = Nvec
  D <- matrix(c(0,5,10,20,5,0,5,10,10,5,0,5,20,10,5,0),4,4)
  D[1,4] = Dvec
  D[lower.tri(D)] <- t(D)[lower.tri(D)]
  pi <- con_probability_per_capita(N=N,D=D,phi=phi,psi=psi)
  return(pi[4])
}

Nvec = seq(1,1000,length.out = 100)
Dvec = seq(0,1000,length.out = 100)
p.D.N.vec <- Vectorize(p.D.N)

z2 <- outer( Dvec,Nvec, p.D.N.vec)
nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(z2, nbcol)
persp3d(psi, Dvec, z2, theta=50, phi=25, expand=0.75, col=color[zcol],
        ticktype="detailed", xlab="Distance", ylab="N", zlab="Probability",axes=TRUE)


# spatial ability check
spatial_con_probability_per_capita <- function(N,D,phi,psi,ppjc,pd){
  if(length(N) == dim(D)[1]){
    dim.N = length(N)
    N.matrix = matrix(0,dim.N,dim.N)
    pi = c()
    for(i in c(1:dim.N)){
      N.matrix[i,]=N/N[i]
      pi[i]=sum(N.matrix[i,]^psi *D[i,]^phi)
    }
    pi = pi/sum(pi)*pd
    return(pi)
  }else{
    print('The length of N should equal the size of D')
  }
}



