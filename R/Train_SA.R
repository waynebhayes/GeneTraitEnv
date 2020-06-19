# Note: R has no 0 index. All vector/matrix indices start from 1

Train_SA<-function(num_exp, G, E, X_tr, Y_tr, GR_tr, obj_best, T_range){
  
  # num_exp (scalar): total number of data points
  # G (matrix): trait-mutation matrix
  # E (matrix): trait-environmental conditionsconditionCall
  # X_tr (matrix): training matrix formed by vectors of mutation. Each column is a mutaion vector of a strain.
  # Y_tr (matrix): training matrix formed by vectors of environmental conditions. Each column is an environmental condition of a strain.
  # GR_tr (vector): growth rates of all the strains. These are "true" values.
  # obj_best (scalar): the threshold of objective function. Once reached, the training ends.
  # T_range (vector): the range of initial temperature and ending temperature.
  
  cat("Cycle ", num_cycle, "Round", num_round, "\n", file="run_report.txt", append=T)
  cat("Cycle ", num_cycle, "Round", num_round, "\n")
  
  
  cat("Cycle ", num_cycle, "Round", num_round, "\n", file="run_report.txt", append=T)
  cat("Cycle ", num_cycle, "Round", num_round, "\n")
  
  # This variable stores the best (lowest) objetive value so far 
  obj<-1e30
  
  
  num_bad<-0
  num_print<-1
  norm_temp<-2
  GR_mean_diff<-0
  traits<-dim(G)[1]
  genes<-dim(G)[2]
  envi<-dim(E)[2]
  
  finished <- 0
  
  # The bad results will be cached in a 1000-element vector
  numCirc <- 1000
  circBuf<-rep(1, numCirc)
  sumCirc <- numCirc
  iCirc <- 0
  
  pBad<-0.1
  
  # Take the initial and ending temperatures.
  SA_Tinit<-T_range[1]
  SA_Tfinal<-T_range[2]
  
  # Each annealing has 2 million steps
  SA_numIters <- 2000000
  
  # Initialize the counter of iterations
  SA_iter <- 0

  SA_T <- SA_Tinit
  SA_lambda = -log(SA_Tfinal/SA_Tinit)
  
  # Merge G and E into a list.
  L<-list(G,E)

  iter<-0
  
  # The vector of simulated growth rates
  GR_guess_tr<-rep(0, length(GR_tr))
  
  # Multiply G and X_tr to generate a matrix representing the genetic strengths of each traits of each strain.
  TG_tr<-G%*%X_tr
  # Multiply E and Y_tr to generate a matrix representing the environmental influences on each traits of each strain.
  TE_tr<-E%*%Y_tr
  
  # The possibilities of G and E being selected to change the element value is correlated to their sizes. 
  part<-genes/(genes+envi)
  
  # Record the G and E with the best performance, and the corresponding correlation coefficient and difference between lab and simulated growth rates

  G_best<-G
  E_best<-E
  best<-c(0,0,0,0)
  names(best)<-c("iter", "obj", "rho", "GRdif")
  

  while (TRUE){
    
    # Randomly determine which element of which matrix will be changed
    current<-runif(1)
    whichRow<-sample(1:traits, 1)
    
    if (current<=part){  
      select<-1
      whichCol<-sample(1:genes, 1)
    }
  
    if (current>part){
        select<-2
        whichCol<-sample(1:envi, 1)  
    }
    
    # Determine the value to be added to the selected matrix element. 
    if (pBad>1e-11){
        delta<-runif(1, -sqrt(pBad), sqrt(pBad))  
    }

    L[[select]][whichRow, whichCol]<-L[[select]][whichRow, whichCol]+delta  
    
    # Possible improvement
    TG_tr<-L[[1]]%*%X_tr
    TE_tr<-L[[2]]%*%Y_tr
    
    # Below is the simpler version for  TG_tr<-L[[1]]%*%X_tr, TE_tr<-L[[2]]%*%Y_tr. Not sure if that's correct.
    
    # if (select==1){
    #   TG_tr[whichRow,]<-TG_tr[whichRow,]+delta*X_tr[whichCol,]
    # }
    # else {
    #   TE_tr[whichRow,]<-TE_tr[whichRow,]+delta*Y_tr[whichCol,]
    # }
    
    # The simulated growth rate of each strain is the dot product of its trait_genetic and its trait_environment vector
    for (i in 1:length(GR_tr)){
      GR_guess_tr[i]<-TG_tr[,i]%*%TE_tr[,i]
    }
    
    # Calculate the pearson correlation coefficient
    rho_pr <- cor(GR_tr, GR_guess_tr, method = "pearson")
    # Calculate the average difference between the simulated 
    GR_mean_diff<-mean(abs(GR_guess_tr-GR_tr))
    
    rho <- rho_pr
    
    # The new objective is calculated
    obj_new <- abs(GR_mean_diff*(1 - rho))
    
  
    # Upadate current annealing temperature according to the number of iterations
    SA_s <- SA_iter / SA_numIters
    SA_iter <- SA_iter + 1
    SA_T <- SA_Tinit * exp (-SA_lambda * SA_s)
    
    # Calculate the possibility of accuracy, and compare it with a random number between 0 and 1
    SA_pAcc <- exp(-(obj_new-obj)/SA_T)
    SA_r <- runif(1)
    
    # A good performance is defined by SA_r<SA_pACC, then the obj variable is updated. If that is "real" improvement by smaller new objective, save the G, E and other relative parameters 
    
    if (SA_r < SA_pAcc) {
        
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
    
    # If not, return the G or E matrix to the original form and insert the corresponding pAcc to the 1000-element queue
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
    
    # The pBad is defined by the average of all the pAcc stored in the queue
    pBad<-sumCirc/numCirc
    
    status_freq <- 10000

    num_print<-num_print+1
    
    # Report the parameters after every 10000 steps
    if(num_print %% status_freq == 0) {
      sumCirc <- 0
      sumCirc<-sum(circBuf)
      cat("iter",c(as.integer(iter+2)), "pBad", pBad, "T", SA_T, "obj", obj, "rho", rho, "GRdif", GR_mean_diff, "\n")
      cat("iter",c(as.integer(iter+2)), "pBad", pBad, "T", SA_T, "obj", obj, "rho", rho, "GRdif", GR_mean_diff, "\n", file="run_report.txt", append=T)
      
    # If the iteration reaches the limit, but the objective is still worse than the threshold, then report the best performance so far
    # Then, use the G and E that generate the best fitting so far as input to restart another annealing process. Try this in 3 cycles. If still fails, restart a new round of training by initializing with new randomly generated G and E matrices. If the falilure lasts after 3 rounds (a total of 9 cycles), raise the threshold of objective by 100% for more tolerant restriction.     

      if (SA_iter >= SA_numIters){
        print("Cannot move further")
        cat("Cannot move further\n", file="run_report.txt", append=T)
        print(best)
        cat("Train:\n","Lab Growth Rates:\n")
        print(GR_tr)
        cat("Dot Production Prediction of Lab Growth Rates:\n")
        print(GR_guess_tr)
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
    
    # If the training reaches the threshold, report and record all important parameters
    if (obj < obj_best)
    {
        print("Good enough")
        cat("Good enough\n", file="run_report.txt", append=T)
        print(c(as.integer(iter), obj, rho, GR_mean_diff))
        cat("iter", c(as.integer(iter), "obj", obj, "rho", rho, "GRdif", GR_mean_diff, "\n"), file="run_report.txt", append=T)
        cat("Train:\n","Lab Growth Rates:\n")
        cat("Train:\n","Lab Growth Rates:\n", file="run_report.txt", append=T)
        print(GR_tr)
        cat(GR_tr, "\n", file="run_report.txt", append=T)
        cat("Dot Product Prediction of Growth Rates:\n", GR_guess_tr, "\n")
        cat("Dot Product Prediction of Growth Rates:\n", GR_guess_tr, "\n", file="run_report.txt", append=T)
        return (list(obj_new, G, E, L, GR_guess_tr, iter, GR_mean_diff)) 
    }

  }
}
