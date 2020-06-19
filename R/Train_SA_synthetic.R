Train_SA_synthetic-function(num_exp, G, E, X_tr, Y_tr, GR_tr, obj_best, T_range){
  
  
  cat("Cycle ", num_cycle, "Round", num_round, "\n", file="run_report.txt", append=T)
  cat("Cycle ", num_cycle, "Round", num_round, "\n")
  
  
  #current<--1
  obj<-1e30
  num_bad<-0
  num_print<-1
  norm_temp<-2
  GR_mean_diff<-0
  traits<-dim(G)[1]
  genes<-dim(G)[2]
  envi<-dim(E)[2]
  
  ### Critical bug fixed: envi<-dim(E)[1] changed to envi<-dim(E)[2]
  
  prev_best <- 2
  finished <- 0
  numCirc <- 1000
  circBuf<-rep(1, numCirc)
  sumCirc <- numCirc
  iCirc <- 0
  # record<-rep(0, 100000)
  # record[1]<-5
  
  pBad<-0.1
  
  SA_Tinit<-T_range[1]
  SA_Tfinal<-T_range[2]
  SA_numIters <- 2000000
  SA_iter <- 0
  
  SA_T <- SA_Tinit
  SA_lambda = -log(SA_Tfinal/SA_Tinit)
  
  
  
  
  L<-list(G,E)
  small_progress<-FALSE
  iter<-0
  GR_guess_tr<-rep(0, length(GR_tr))
  TG_tr<-G%*%X_tr
  TE_tr<-E%*%Y_tr
  TG[TG!=1]<-1
  #TG_tr_n<-TG_tr
  #TE_tr_n<-TE_tr
  part<-genes/(genes+envi)
  G_best<-G
  E_best<-E
  best<-c(0,0,0,0)
  names(best)<-c("iter", "obj", "rho", "GRdif")
  
  
  while (TRUE){
    current<-runif(1)
    whichRow<-sample(1:traits, 1)
    #if (current%%2==0){
    if (current>1){  
      select<-1
      whichCol<-sample(1:genes, 1)
    }
    #if (current%%2==1){
    if (current<=1){
      select<-2
      whichCol<-sample(1:envi, 1)  
    }
    if (pBad>1e-11){
      delta<-runif(1, -sqrt(pBad), sqrt(pBad))  
    }
    # Note: for this version, only change E.
    
    #if (current%%3==0){
    #    select<-3
    #    delta_B<-runif(1, -0.1, 0.1)
    #}
    
    L[[select]][whichRow, whichCol]<-L[[select]][whichRow, whichCol]+delta  
    
    
    
    #if (select==1){
    #  TG_tr_n[whichRow,]<-TG_tr[whichRow,]+delta*X_tr[whichCol,]
    #}
    #else {
    #  TE_tr_n[whichRow,]<-TE_tr[whichRow,]+delta*Y_tr[whichCol,]
    #}
    
    if (select==1){
      TG_tr[whichRow,]<-TG_tr[whichRow,]+delta*X_tr[whichCol,]
    }
    else {
      TE_tr[whichRow,]<-TE_tr[whichRow,]+delta*Y_tr[whichCol,]
    }
    
    
    # Possible improvement
    #TG_tr<-L[[1]]%*%X_tr
    TE_tr<-L[[2]]%*%Y_tr
    
    #if (SA_iter<2){
    #  if (select==1){
    #    print(which(TG_tr!=TG_tr_n, arr.ind=T))  
    #    cat("TG", whichrow, SA_iter, "\n")
    #    print(TG_tr_n[whichRow,1:10])
    #    print(TG_tr[whichRow,1:10])
    #    print(all(TG_tr_n[whichRow,1:10]==TG_tr[whichRow,1:10]))
    #  }
    #  else{
    #    print(which(TE_tr!=TE_tr_n, arr.ind=T))  
    #    cat("TE", whichRow, SA_iter, "\n")
    #   print(TE_tr_n[whichRow,1:10])
    #    print(TE_tr[whichRow,1:10])
    #    print(all(TE_tr_n==TE_tr))
    #  }
    
    #}
    
    
    for (i in 1:length(GR_tr)){
      # GR_guess[i]<-(TG[,i]+L[[3]])%*%TE[,i]
      GR_guess_tr[i]<-TG_tr[,i]%*%TE_tr[,i]
    }
    
    # GR_guess<-GR_guess/max(GR_guess)
    # norm_max<- max(abs(GR-GR_guess/max(abs(GR_guess))))
    # norm_mean<- mean(abs(GR-GR_guess/max(abs(GR_guess))))
    rho_sp <- cor(GR_tr, GR_guess_tr, method = "spearman")
    rho_pr <- cor(GR_tr, GR_guess_tr, method = "pearson")
    
    GR_mean_diff<-mean(abs(GR_guess_tr-GR_tr))
    
    rho <- (0*rho_sp + 2*rho_pr)/2
    # norm <- (4*norm_mean + 0*norm_max)/4
    
    
    # obj_diff<-abs(GR_mean_diff)
    # obj_rho<-abs(1-rho)
    # obj_df<-abs(num_exp-sum(sqrt(abs(L[[1]])))-sum(sqrt(abs(L[[2]]))))/num_exp
    # obj_new<-obj_diff + obj_rho + obj_df
    obj_new <- abs(GR_mean_diff^1.5*(1 - rho))
    # obj_new = (0*(norm -rho) + 4*norm*(1 - rho))/4
    
    
    SA_s <- SA_iter / SA_numIters
    SA_iter <- SA_iter + 1
    SA_T <- SA_Tinit * exp (-SA_lambda * SA_s)
    SA_pAcc <- exp(-(obj_new-obj)/SA_T)
    SA_r <- runif(1)
    
    if(is.na(SA_r < SA_pAcc)){
      cat("iter",c(as.integer(iter+2)), "SA_r", SA_r, "SA_s", SA_s, "SA_lambda", SA_lambda, "\n")
    }
    
    if (SA_r < SA_pAcc) {
      # norm_temp<-norm
      if (obj>obj_new){
        G_best<-L[[1]]
        E_best<-L[[2]]
        best["iter"]<-iter
        best["obj"]<-obj_new
        best["rho"]<-rho
        best["GRdif"]<-GR_mean_diff
      }
      
      obj<-obj_new
    }
    else {
      L[[select]][whichRow,whichCol]<-L[[select]][whichRow,whichCol]-delta  
      iCirc <- iCirc + 1
      if(iCirc>numCirc){
        iCirc<-0
      }
      sumCirc <- sumCirc - circBuf[iCirc %% numCirc + 1]
      circBuf[iCirc %% numCirc + 1] <- SA_pAcc
      sumCirc <- sumCirc + SA_pAcc
    }
    
    pBad<-sumCirc/numCirc
    
    status_freq <- 10000
    
    num_print<-num_print+1
    
    if(num_print %% status_freq == 0) {
      sumCirc <- 0
      sumCirc<-sum(circBuf)
      cat("iter",c(as.integer(iter+2)), "pBad", pBad, "T", SA_T, "obj", obj, "rho", rho, "GRdif", GR_mean_diff, "\n")
      cat("iter",c(as.integer(iter+2)), "pBad", pBad, "T", SA_T, "obj", obj, "rho", rho, "GRdif", GR_mean_diff, "\n", file="run_report.txt", append=T)
      # record[num_print/status_freq+1]<-obj
      if (SA_iter >= SA_numIters) { #} abs(record[num_print/status_freq+1]-record[num_print/status_freq])<0.0001){
        print("Cannot move further")
        cat("Cannot move further\n", file="run_report.txt", append=T)
        print(best)
        cat("Train:\n","Lab Growth Rates:\n")
        print(GR_tr)
        cat("Dot Production Prediction of Lab Growth Rates:\n")
        print(apply((G_best%*%X_tr)*(E_best%*%Y_tr), 2, sum))
        cat("iter", best["iter"], "obj", best["obj"], "rho", best["rho"], "GRdif", best["GRdif"], "\n", file="run_report.txt", append=T) 
        
        cat("Dot Production Prediction of Lab Growth Rates:\n", apply((G_best%*%X_tr)*(E_best%*%Y_tr), 2, sum), "\n", file="run_report.txt", append=T)
        cat("Lab Growth Rates:\n", GR_tr, "\n", file="run_report.txt", append=T)
        print(best["obj"])
        cat("Best obj is larger than the goal\n")
        num_cycle<<-num_cycle+1
        
        if (num_cycle>3){
          
          num_cycle<<-1
          num_round<<-num_round+1
          
          if (num_round>3){
            
            print("Cannot get a good solution. Raise the goal of obj by 100%")
            obj_best<-obj_best*2
            num_round<<-1
            cat("Now goal of obj is:", obj_best, "\n")
            cat("Now goal of obj is:", obj_best, "\n", file="run_report.txt", append=T)
            
            cat("Round ", num_round, "\n", "Goal of obj ", obj_best, "\n", file="run_report.txt", append=T)
            cat("Round ", num_round, "\n", "Goal of obj ", obj_best, "\n")
            
          } 
          
          cat("G, E redefined. Anneal again with new G and E as input\n", file="run_report.txt", append=T)
          print("G, E redefined. Anneal again with new G and E as input")
          
          G_initial<<-matrix(rnorm(traits*genes), traits, genes)
          E_initial<<-matrix(rnorm(traits*envi), traits, envi)
          
          G<-G_initial
          E<-E_initial
          
          cat("Initial G:\n", G, "\n", file="run_result.txt", append=T)
          cat("Initial E:\n", E, "\n", file="run_result.txt", append=T)
          
          run_again<-Train_SA(num_exp, G, E, X_tr, Y_tr, GR_tr, obj_best, T_range)
        }
        
        else {
          cat("Anneal again with the best G and E so far as input\n", file="run_report.txt", append=T)
          print("Anneal again with the best G and E so far as input")
          run_again<-Train_SA(num_exp, G_best, E_best, X_tr, Y_tr, GR_tr, obj_best, T_range)  
        }
        
        return (run_again) 
      }
    }  
    
    
    iter<-iter+1
    current<-current+1
    
    if (obj < obj_best)# (num_bad > 200 | obj < obj_best) #rho>rho_best & norm_temp<norm_temp)
    {
      print("Good enough")
      cat("Good enough\n", file="run_report.txt", append=T)
      print(c(as.integer(iter), obj, rho, GR_mean_diff))
      cat("iter", c(as.integer(iter), "obj", obj, "rho", rho, "GRdif", GR_mean_diff, "\n"), file="run_report.txt", append=T)
      cat("Train:\n","Lab Growth Rates:\n")
      cat("Train:\n","Lab Growth Rates:\n", file="run_report.txt", append=T)
      print(GR_tr)
      cat(GR_tr, "\n", file="run_report.txt", append=T)
      cat("Dot Product Prediction of Growth Rates:\n")
      cat("Dot Product Prediction of Growth Rates:\n", file="run_report.txt", append=T)
      cat(GR_guess_tr, "\n", file="run_report.txt", append=T)
      print(GR_guess_tr)
      return (list(obj_new, G, E, L, GR_guess_tr, iter, GR_mean_diff)) 
    }
    
  }
}
