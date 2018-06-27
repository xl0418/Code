text = NULL
phi = c(0,0.001,0.01,0.1,1)
psi = phi
for(i in c(1:5)){
  for(j in c(1:5)){
    text1 = paste0("phi=", phi[i]," psi=",psi[j])
    text = c(text,text1)
  }
}

fileConn<-file("c:/Liang/Code/Pro3/R_p3/parameterfile.txt")
writeLines(text, fileConn)
close(fileConn)