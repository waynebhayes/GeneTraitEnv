GetpBad<-function(G, E, X_tr, Y_tr, GR_tr, SA_T, threshold){
  
  N<-10000
  counter<-0
  counter_bad<-0
  pbad_sum<-0
  
  current<--1
  obj<-1e30
  num_bad<-0
  GR_mean_diff<-0
  
  genes<-dim(G)[2]
  envi<-dim(E)[2]
  traits<-dim(G)[1]
  prev_best <- 2

  L<-list(G,E)
  GR_guess_tr<-rep(0, length(GR_tr))
  TG_tr<-G%*%X_tr
  TE_tr<-E%*%Y_tr
  
  while (counter<N){
    
    whichRow<-sample(1:traits, 1)
    if (current%%2==0){
      select<-1
      whichCol<-sample(1:genes, 1)
    }
    if (current%%2==1){
      select<-2
      whichCol<-sample(1:envi, 1)  
    }
    
    delta<-runif(1, -sqrt(threshold)/1.2, sqrt(threshold)/1.2)  
    L[[select]][whichRow, whichCol]<-L[[select]][whichRow, whichCol]+delta  
    
    if (select==1){
      TG_tr[whichRow,]<-TG_tr[whichRow,]+delta*X_tr[whichCol,]
    }
    else {
      TE_tr[whichRow,]<-TE_tr[whichRow,]+delta*Y_tr[whichCol,]
    }
  
    # Possible improvement
    TG_tr<-L[[1]]%*%X_tr
    TE_tr<-L[[2]]%*%Y_tr
    
    for (i in 1:length(GR_tr)){
      GR_guess_tr[i]<-TG_tr[,i]%*%TE_tr[,i]
    }
    
    rho_sp <- cor(GR_tr, GR_guess_tr, method = "spearman")
    rho_pr <- cor(GR_tr, GR_guess_tr, method = "pearson")
    
    GR_mean_diff<-mean(abs(GR_guess_tr-GR_tr))
    
    rho <- (0*rho_sp + 2*rho_pr)/2
    
    obj_diff<-abs(GR_mean_diff)
    obj_new = abs(GR_mean_diff*(1 - rho))
    SA_pAcc<-exp(-(obj_new-obj)/SA_T)
    SA_r <- runif(1)
    
    if (SA_r < SA_pAcc) {
      obj<-obj_new
    }
    else {
      L[[select]][whichRow,whichCol]<-L[[select]][whichRow,whichCol]-delta  
      counter_bad<-counter_bad+1
      pbad_sum<-pbad_sum+SA_pAcc
    }
   
    current<-current+1
    counter<-counter+1
  }
  
  if (counter_bad==0){
    return(1)
  }
  else {
    return(pbad_sum/counter_bad)    
  }

  
  
}