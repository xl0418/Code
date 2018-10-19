library(DDD)
# Arguments introduction
# n: number of locations in the model. Multiple locations model is already incorporated.
# parsN: a vector of number of species on each location at initial stage. In total, the number of species
#        in the system should be 2.
# For example: if set n = 2, parsN can be set as parsN = c(1,1) which means 1 species for each location
#              parsN = c(2,0) which means 2 species on location 1 and 0 species on location 2.
#              In this case, we can set migration rate to 0 to approach our spatial model to non-spatial model.
# age: crown age. As the simulation is conditioning on the survival of two ancestor specis, the crown age is fixed.
# pars: a vector of parameters. (lambda, mu, K'( or K)).
#              lambda: intrisic speciation rate
#              mu: extinction rate
#              K': ecological limit when argument K_fix = 0
#              K : carrying capacity when argument K_fix = 1
# seed_fun: seed for the simulation.
# lambda_allo0: intrisic allopatric speciation rate.
# M0: intrisic migration rate.
# K_fix: fix carrying capacity (K_fix = 1) or fix ecological limit (K_fix = 0)

sddsim<-function(n,parsN,age=20,lambda_n,mu_n,K_n , seed_fun = 29, lambda_allo0 = 0.2, M0=1,K_fix = 1){
  set.seed(seed_fun)
  if (length(parsN) == n){
    if(sum(parsN) ==2 | sum(parsN)==1){
      done = 0
      while(done == 0)
      {
        event_info_track = NULL
        brts = c(0,"start")
        lambda0 = lambda_n
        mu0 = mu_n

        #fix K, the carrying capacity
        if(K_fix == 1) {  K = K_n   # K is carrying capacity
        K= lambda0*K/(lambda0 - mu0) # here it is K'
        }
        #fix K', the limit to the diversity
        if(K_fix == 0){K = K_n}

        if(length(K) == 1){
          K_loc = rep(K_n,n)   # carrying capacity for each location
        }
        if(length(K) > 1){
          K_loc = K   # carrying capacity for each location
        }

        probs_record = NULL
        t <- 0
        # Na : number of species in location A
        # Nb : number of species in location B
        # Nab: number of species in location AB
        # age : time
        # Ki : carrying capacity of location i
        # lambda0 : initial speciation rate
        # mu0 : initial extinction rate


        # string of events
        #track locations as the maps of the "phylo" class
        track_loc = list()
        for(i in 1:n){
          track_loc[[i]]=0
          names(track_loc[[i]])=c(paste(i))
        }
        N=sum(parsN)
        Ntable_index = Nindex(n)
        Nlength = sum(Ntable_index)
        Ntable = matrix(0,nrow =1, ncol =Nlength)
        Ntable[1:n] = parsN
        i=0
        # L : Ltable used in L2phylo function of DDD package
        # L = data structure for lineages,
        # . L[,1] = branching times
        # . L[,2] = index of parent species
        # . L[,3] = index of daughter species
        # . L[,4] = time of extinction
        # . L[,5] = index of location
        # j = index running through L

        loctable= NULL

        L = NULL
        for(j in 1:N){
          loc = which(Ntable[1,]!=0)
          L = rbind(L, c(0, 1-j, (-1)^j*j, -1, loc[1]))
          L = matrix(L, ncol = 5)
          loc1= matrix(0,1,n)
          loc1[1,loc[1]] = 1
          loctable = rbind(loctable, loc1)
          linlist = cbind(L[,3], L[,5])
          Ntable[1,loc[1]] = Ntable[1,loc[1]] -1
          newL=j
        }

        Ntable = matrix(0,nrow =1, ncol =Nlength)
        Ntable[1:n] = parsN

        spec_num = sum(Nindex(n)[1,])
        N_loc = matrix(0,nrow = n,ncol = spec_num)
        N_loc_col = matrix(0,nrow = n,ncol = spec_num)
        Ndistribution = event_matrix(n)

        #index for sym speciation
        for(j in 1:n){
          N_loc_col[j,] = which(Ndistribution[j,] == 1,arr.ind = TRUE)#index of each loc in sym spec table
        }
        B_symspec = c(N_loc_col)

        #index for extinction
        B_ext = NULL
        x = Ndistribution
        y = split(x, rep(1:ncol(x), each = nrow(x)))
        for(j in 1:n){
          x1 = x
          x1[j,] = x1[j,]-1
          y1 = split(x1, rep(1:ncol(x1), each = nrow(x1)))
          z = match(y1,y)
          z = z[!is.na(z)]
          B_ext = rbind(B_ext, z)
        }
        B_ext = c(rep(0,n), B_ext)

        #index for allo speciation
        lambda_allo = NULL
        allo_dau1 = NULL
        allo_dau2 = NULL
        allo_par = NULL
        N_allo = colSums(Ndistribution)
        for(j in 1:(n-1)){
          allo_index = which(N_allo > 1 )
          allo_matrix = Ndistribution[,allo_index,drop = FALSE]
          allo_matrix_daughter = matrix(0,nrow = nrow(allo_matrix),ncol = ncol(allo_matrix))
          allo_matrix_daughter1 = allo_matrix_daughter
          allo_matrix_daughter1[1:j,] = allo_matrix[1:j,]
          allo_matrix_daughter2 = allo_matrix
          allo_matrix_daughter2[1:j,] = 0
          allo_col_cut1 = which(colSums(allo_matrix_daughter1) == 0)
          allo_col_cut2 = which(colSums(allo_matrix_daughter2) == 0)
          allo_col_cut = c(allo_col_cut1,allo_col_cut2)

          if(length(allo_col_cut) != 0) {
            allo_matrix_daughter1_rest = allo_matrix_daughter1[,-allo_col_cut]
            allo_matrix_daughter2_rest = allo_matrix_daughter2[,-allo_col_cut]
            allo_index_final = allo_index[-allo_col_cut]
          }
          else {
            allo_matrix_daughter1_rest = allo_matrix_daughter1
            allo_matrix_daughter2_rest = allo_matrix_daughter2
            allo_index_final = allo_index

          }
          allo_dau1 = cbind(allo_dau1,allo_matrix_daughter1_rest)
          allo_dau2 = cbind(allo_dau2,allo_matrix_daughter2_rest)
          allo_par = c(allo_par, allo_index_final)
        }
        B_allodau1 = match(data.frame(allo_dau1),data.frame(Ndistribution))
        B_allodau2 = match(data.frame(allo_dau2),data.frame(Ndistribution))
        allo_index = allo_par


        # ignore the linear topology of locations
        matrix_allo = rbind(B_allodau1,B_allodau2,allo_index)
        duplicated.columns <- duplicated(t(matrix_allo))
        new.matrix <- matrix_allo[, !duplicated.columns,drop = FALSE]
        B_allodau1 = new.matrix[1,]
        B_allodau2 = new.matrix[2,]
        allo_index = new.matrix[3,]


        #index for migration
        B_mig_target = which(Ndistribution == 0, arr.ind = TRUE)
        B_mig_from = c(B_mig_target[,2])
        B_mig_to = c(B_mig_target[,1])

        B_mig_bec = NULL
        for(j in 1:length(B_mig_from)){
          Ndis_aftermig = Ndistribution
          Ndis_aftermig[B_mig_to[j],B_mig_from[j]] = 1
          loc_aftermig = as.vector(Ndis_aftermig[,B_mig_from[j]])
          B_mig_bec = c(B_mig_bec,which(sapply( y ,function(x)all(x==loc_aftermig))))
        }
        B_mig_bec = as.numeric(B_mig_bec)

        #number of events
        #sym spec
        num_ss = 0
        for(i in 1:(n-1))
        {
          num_ss = num_ss+choose(n-1,i)
        }
        num_ss = n*(1+num_ss)

        #ext
        num_ext = num_ss

        #allo spec
        num_as = length(allo_index)

        #mig
        num_mig = 0
        for(i in 1:(n-1)){
          num_mig = num_mig + choose(n,i)*(n-i)
        }

        #number of events
        num_event = num_ss + num_ext + num_as + num_mig
        probs_part1 = num_ss
        probs_part2 = num_ext
        probs_part3 = num_as
        probs_part4 = num_mig


        B<- c(1:num_event)

        i = 0
        log_likelihood = 0
        while(t[i+1]< age){
          i<-i+1

          # speciation event & extinction event
          lambda_sym = rep(0,n)

          if(length(lambda0) == 1) lambda0 = rep(lambda0,n)

          if(length(mu0) == 1){ mu = rep(mu0,n)
          }else mu = mu0

          sym_spec_event = matrix(0,nrow = n,ncol = spec_num)
          ext_event = matrix(0,nrow = n,ncol = spec_num)

          for(j in 1:n){
            N_loc[j,] = Ntable[i,which(Ndistribution[j,]==1)]  #number of each loc in sym spec table
            lambda_sym[j]=max(lambda0[j]*(1-sum(N_loc[j,])/K_loc[j]),0)
            sym_spec_event[j,] = lambda_sym[j]*N_loc[j,]
            ext_event[j,] = mu[j] *N_loc[j,]
          }
          prob_spec_sym = c(sym_spec_event)
          prob_ext = c(ext_event)


          # Migration
          prob_mig = NULL
          Mig_dir = rep(0,n)
          for(j in 1:n){
            Mig_dir[j] = max(M0*(1-sum(N_loc[j,])/K_loc[j]),0)
          }
          for(j in 1:(ncol(Ntable)-1)){
            tar = which(Ndistribution[,j]==0)
            prob_mig_each = Ntable[i,j]*Mig_dir[tar]
            prob_mig_each = matrix(prob_mig_each,ncol = length(prob_mig_each))
            prob_mig = cbind(prob_mig,prob_mig_each)
          }
          prob_mig = c(prob_mig)


          # Allopatric speciation event
          lambda_allo = NULL
          for(j in allo_index){
            Mig_base =  M0
            if(Mig_base == 0) lambda_allo_each = 1*Ntable[i,j]
            else lambda_allo_each = max(lambda_allo0/Mig_base, 0 )*Ntable[i,j]
            lambda_allo = cbind(lambda_allo,lambda_allo_each)
          }
          prob_spec_allo = c(lambda_allo)

          #probs of all events
          probs= c(prob_spec_sym,prob_ext,prob_spec_allo,prob_mig)
          probs_record = rbind(probs_record,probs)
          #Total rate
          TR =sum(probs)

          if(TR ==0) break
          else{
            A<-DDD::sample2(B,1, prob = probs)
            t_wait = rexp(1,rate=TR)
            t[i+1]=t[i]+t_wait

            prob_waiting_time = 1-exp(-TR*t_wait)

            if(t[i+1]>age) break
            loc2 = matrix(0,1,n)
            log_likelihood = log_likelihood + log(prob_waiting_time) + log(probs[A]/TR)

            # Sympatric speciation
            if (is.element(A,1:probs_part1)){

              b1<-A%%n
              if(b1 == 0) b1 = n
              Ntable=rbind(Ntable,Ntable[i,])
              Ntable[i+1,b1] = Ntable[i,b1]+1
              newL = newL + 1;
              list0 = matrix(linlist,ncol = 2)
              b3 <- B_symspec[A]

              list1 = linlist[list0[,2]== b3]

              list2 = matrix(list1, ncol = 2)

              linlist1 = list2[,1]

              ranL= DDD::sample2(linlist1,1)
              L = rbind(L,c(t[i+1],ranL,sign(ranL) * newL,-1,b1))
              linlist = rbind(linlist,c(sign(ranL) * newL,b1))
              loc2[1,b1] = 1
              loctable = rbind(loctable, loc2)
              #change of loc info
              track_loc[[abs(ranL)]]=c(track_loc[[abs(ranL)]],t[i+1])
              names(track_loc[[abs(ranL)]])[length(names(track_loc[[abs(ranL)]]))]=c(paste(b3))
              track_loc[[newL]]=track_loc[[abs(ranL)]]
              names(track_loc[[newL]])[length(names(track_loc[[newL]]))]=c(paste(b1))

              event_info_track = c(event_info_track, c(paste0("t",abs(ranL)," gave birth to ",newL)))
              brts = rbind(brts, c(t[i+1],paste0("t",abs(ranL)," gave birth to ","t",newL)))
            }

            #Extinction
            else if(is.element(A,(probs_part1+1):(probs_part1+probs_part2))) {


              b1<-A%%n
              if(b1 == 0) b1 = n
              Ntable=rbind(Ntable,Ntable[i,])
              b3 <- B_symspec[A-probs_part1]
              list0 = matrix(linlist,ncol = 2)
              list1 = linlist[list0[,2]== b3]
              list2 = matrix(list1, ncol = 2)
              linlist1 = list2[,1]
              ranL= DDD::sample2(linlist1,1)
              loctable[abs(ranL),b1] = 0
              # Go extinct
              if(is.element(A,(probs_part1+1):(probs_part1+n))){


                Ntable[i+1,b3] = max(Ntable[i,b3]-1,0)
                L[abs(ranL),4] = t[i+1]
                if(length(L[L[,4]== -1]) == 0) break
                else {
                  v = which(linlist[,1] == ranL)
                  linlist = linlist[-v,,drop=FALSE]
                  linlist = linlist[order(linlist[,1]),]
                  linlist = matrix(linlist,ncol=2)

                  track_loc[[abs(ranL)]]=c(track_loc[[abs(ranL)]],t[i+1])
                  names(track_loc[[abs(ranL)]])[length(names(track_loc[[abs(ranL)]]))]=c("0")
                  event_info_track = c(event_info_track, c(paste0("t",abs(ranL)," went extinct")))

                }
              }
              # Contraction
              else{

                b2 <- B_ext[A-probs_part1]
                Ntable[i+1,b3] = Ntable[i,b3]-1
                Ntable[i+1,b2] = Ntable[i,b2]+1
                L[abs(ranL),5] <- b2
                v = which(linlist[,1] == ranL)
                linlist[v,2] = b2
                linlist = linlist[order(linlist[,1]),]
                linlist = matrix(linlist,ncol=2)


                track_loc[[abs(ranL)]]=c(track_loc[[abs(ranL)]],t[i+1])
                names(track_loc[[abs(ranL)]])[length(names(track_loc[[abs(ranL)]]))]=c(paste(b2))

                event_info_track = c(event_info_track, c(paste0("t",abs(ranL)," contracted to ",b2)))

              }

            }


            # Allopatric speciation
            else if(is.element(A,((probs_part1+probs_part2)+1):(probs_part1+probs_part2+probs_part3))) {

              A1 = A - (probs_part1+probs_part2)

              Ntable=rbind(Ntable,Ntable[i,])
              Ntable[i+1,B_allodau1[A1]] = Ntable[i,B_allodau1[A1]]+1
              Ntable[i+1,B_allodau2[A1]] = Ntable[i,B_allodau2[A1]]+1
              Ntable[i+1,allo_index[A1]] = Ntable[i,allo_index[A1]]-1
              newL = newL + 1
              list0 = matrix(linlist,ncol = 2)
              list1 = linlist[list0[,2]==allo_index[A1]]
              list2 = matrix(list1, ncol = 2)
              linlist1 = list2[,1]
              ranL= DDD::sample2(linlist1,1)
              L[abs(ranL),5] <- B_allodau1[A1]
              loctable[abs(ranL), ] = Ndistribution[,B_allodau1[A1]]
              v = which(linlist[,1] == ranL)
              linlist[v,2] = B_allodau1[A1]
              L = rbind(L,c(t[i+1],ranL,sign(ranL) * newL,-1,B_allodau2[A1]))
              linlist = rbind(linlist,c(sign(ranL) * newL,B_allodau2[A1]))
              loc2[1,] = Ndistribution[,B_allodau2[A1]]
              loctable = rbind(loctable,loc2)
              #change of loc info
              track_loc[[abs(ranL)]]=c(track_loc[[abs(ranL)]],t[i+1])
              names(track_loc[[abs(ranL)]])[length(names(track_loc[[abs(ranL)]]))]=c(paste(B_allodau1[A1]))
              track_loc[[newL]]=track_loc[[abs(ranL)]]
              names(track_loc[[newL]])[length(names(track_loc[[newL]]))]=c(paste(B_allodau2[A1]))
              event_info_track = c(event_info_track, c(paste("t",abs(ranL)," gave birth to ",B_allodau1[A1]," & ",B_allodau2[A1])))
              brts = rbind(brts, c(t[i+1],paste0("t",abs(ranL)," gave birth to ",B_allodau1[A1]," & ",B_allodau2[A1])))

            }

            #Migration
            else {


              Mig1 = B_mig_from
              Mig2 =  B_mig_bec
              Mig3 = B_mig_to
              Am = A - (probs_part1+probs_part2+probs_part3)
              b2 = Mig1[Am]
              b3 = Mig2[Am]
              Ntable=rbind(Ntable,Ntable[i,])
              Ntable[i+1,b2]=Ntable[i,b2]-1
              Ntable[i+1,b3]=Ntable[i,b3]+1
              linlist = matrix(linlist,ncol = 2)
              list1 = linlist[linlist[,2]==b2]
              list2 = matrix(list1, ncol = 2)
              linlist1 = list2[,1]
              ranL1 = DDD::sample2(linlist1,1)
              L[abs(ranL1),5] <- b3
              v = which(linlist[,1] == ranL1)
              linlist[v,2] = b3
              linlist = linlist[order(linlist[,1]),]
              linlist = matrix(linlist,ncol=2)
              loctable[abs(ranL1),Mig3[Am]] = 1
              #change of loc info
              track_loc[[abs(ranL1)]]=c(track_loc[[abs(ranL1)]],t[i+1])
              names(track_loc[[abs(ranL1)]])[length(names(track_loc[[abs(ranL1)]]))]=c(paste(b3))
              event_info_track = c(event_info_track, c(paste(abs(ranL1)," dispersed and becomes",b3)))

            }

            if(sum(linlist[,1] < 0) == 0 | sum(linlist[,1] > 0) == 0)
            {
              break
            }

            N[i+1]=sum(Ntable[i+1,])

          }
        }

        if(sum(linlist[,1] < 0) == 0 | sum(linlist[,1] > 0) == 0)
        {
          done = 0

        } else {
          done = 1
        }
      }

      if(dim(L)[1]==1)
        print(paste("It goes extinct at the begining!"))
      else{
        L[,1]=age-L[,1]
        notmin1 = which(L[,4] != -1)
        L[notmin1,4] = age - c(L[notmin1,4])
        L[which(L[,4] == age + 1),4] = -1

        tes = L2phylo(L,dropextinct = T)
        tas = L2phylo(L,dropextinct = F)
        out = list(tes = tes,tas = tas,L = L, loctable = loctable, Ntable = Ntable,t = t,log_likelihood = log_likelihood,probs_record = probs_record,track_loc = track_loc,event_info_track=event_info_track,brts = brts)

        return(out)
      }
    }
    else print(paste("please input 1 or 2 species in total!"))
  }
  else print(paste("please match the length of parsN to n"))

}
print(paste("Please input initial number of lineages of each loc with no more than 2 species in total! Match number n and dimentions of parsN "))

