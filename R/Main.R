Main<-function(input, basal, seed, traits, K, obj_best){
  
  DB<-Preprocessing(input)
  
  
  val_result<-list("Lab"=c(), "Simu"=c())
  tr_result<-list("Lab"=c(), "Simu"=c())
  G_each_K<-c()
  E_each_K<-c()
  
  colnames_total<-colnames(DB)
  col_gene<-(grep("Gene_start", colnames_total)+1):(grep("Growth_Rate", colnames_total)-1)
  col_gene_start<-grep("Gene_start", colnames_total)+1
  col_env<-(grep("Env_start", colnames_total)+1):(grep("Gene_start", colnames_total)[1]-1)
  
  num_exp<-nrow(DB)
  
  sample_ID<-DB[, "sample_ID"]
  DB_env<-DB[ , col_env]
  DB_gene<<-DB[ , col_gene]
  GR_lab<-DB[, "Growth_Rate"]
  GR_basal<-DB[, "Basal_Growth_Rate"]
  X<-t(DB_gene)
  X[is.na(X)]<-0
  Y<-t(DB_env)
  Y[is.na(Y)]<-0

  if (basal){
    GR<-GR_lab-GR_basal  
  }
  else{
    GR<-GR_lab
  }
  
  paperdiv<-TRUE
  
  output<-Kfold(seed, K, paperdiv, sample_ID)
  tr_whichrow<-output[[1]]
  va_whichrow<-output[[2]]
  exp_index_shuffle<-output[[3]]
  
  result<-list()
  
  norm_best<-0.03
  rho_best<-0.8
  
  
  score<-vector("numeric")
  validation_score<-vector("numeric")
  
  GR_va_guess<-c()
  GR_va_all<-c()
  result_train<- vector(mode = "list", length = K)
  cat("seed=", seed, "traits=", traits, "\n", file="run_report.txt")
  cat("seed=", seed, "traits=", traits, "\n", file="run_result.txt")
  
  for (i in 1:K){
    
    G_initial<<-matrix(rnorm(traits*nrow(X)), traits, nrow(X))
    E_initial<<-matrix(rnorm(traits*nrow(Y)), traits, nrow(Y))
    
    G<-G_initial
    E<-E_initial
    
    cat("K=", i)
    cat("K=", i, "\n", file="run_report.txt", append=T)
    
    
    cat("Initial G:\n", G, "\n", file="run_result.txt", append=T)
    cat("Initial E:\n", E, "\n", file="run_result.txt", append=T)
    
    
    cat("K=", i, "\n", file="run_result.txt", append=T)
    cat("K=", i, "\n", file="run_result_top.txt", append=T)
    
    
    tr_index<-exp_index_shuffle[tr_whichrow[[i]]]
    va_index<-exp_index_shuffle[va_whichrow[[i]]]
    X_tr<-X[, tr_index]
    Y_tr<-Y[, tr_index]
    X_va<-X[, va_index]
    Y_va<-Y[, va_index]
    GR_tr<-GR[tr_index]
    GR_va<-GR[va_index]

    SA_Tinit<-1
    threshold<-0.5
    
    num_cycle<<-1
    num_round<<-1
    
    cat("Round ", num_round, "\n", "Goal of obj ", obj_best, "\n", file="run_report.txt", append=T)
    cat("Round ", num_round, "\n", "Goal of obj ", obj_best, "\n")
    
    
    
    while(GetpBad(G, E, X_tr, Y_tr, GR_tr, SA_Tinit, threshold)<threshold){
      SA_Tinit<-SA_Tinit*10
    }
    SA_Tinit<-SA_Tinit/2
    while(GetpBad(G, E, X_tr, Y_tr, GR_tr, SA_Tinit, threshold)>threshold){
      SA_Tinit<-SA_Tinit/2
    }

    SA_Tfinal<-SA_Tinit/10
    
    threshold<-1e-6
    while(GetpBad(G, E, X_tr, Y_tr, GR_tr, SA_Tfinal, threshold)>threshold){
      SA_Tfinal<-SA_Tfinal/10
    }
    SA_Tfinal<-SA_Tfinal*2
    while(GetpBad(G, E, X_tr, Y_tr, GR_tr, SA_Tfinal, threshold)<threshold){
      SA_Tfinal<-SA_Tfinal*2
    }

    # T_range<-c(5e-1, 1e-11)
    T_range<-c(SA_Tinit, SA_Tfinal)
    print(T_range)  
    
    cat("SA_init=", SA_Tinit, "SA_Tfinal", SA_Tfinal, "\n", file="run_report.txt", append=T)
     
    result_train[[i]]<-Train_SA(num_exp, G, E, X_tr, Y_tr, GR_tr, obj_best, T_range)
    
    score<-append(score, result_train[[i]][[1]])
    
    G<-result_train[[i]][[4]][[1]]
    cat("G=", t(G), "\n", file="run_result.txt", append=T)
    
    
    #G_prime<-G
    #G_prime[abs(G_prime)<sort(abs(G_prime), decreasing=TRUE)[length(sample_ID)]]<-0
    #cat("G_prime:\n", t(G_prime), "\n", file="run_result_top.txt", append=T)
    
    #if (i>1){
    #  G_prime_same_idx<-None_zero_idx_same(G_prime_prev, G_prime)
    #  cat("Non-zero element positions where G_prime shares with the previous one:\n", file="run_result_top.txt", append=T)
    #  write.table(G_prime_same_idx,  file="run_result_top.txt", append=T)
      
    #}
    #G_prime_prev<-G_prime
    
    
    
    E<-result_train[[i]][[4]][[2]]
    cat("E=", t(E), "\n", file="run_result.txt", append=T)
    #E_prime<-E
    # if (length(sample_ID)<nrow(E_prime)*ncol(E_prime)){
    #  E_prime[abs(E_prime)<sort(abs(E_prime), decreasing=TRUE)[length(sample_ID)]]<-0
    # }
    
    
    # cat("E_prime:\n", E_prime, "\n", file="run_result_top.txt", append=T)
    # if (i>1){
    #   E_prime_same_idx<-None_zero_idx_same(E_prime_prev, E_prime)
    #   cat("Non-zero element positions where E_prime shares with the previous one:\n", file="run_result_top.txt", append=T)
    #   write.table(E_prime_same_idx,  file="run_result_top.txt", append=T)
    #   
    # }
    # E_prime_prev<-E_prime
    # 
    GR_tr_guess<-result_train[[i]][[5]]
    
    tr_result$Lab<-append(tr_result$Lab, GR_tr)
    tr_result$Simu<-append(tr_result$Simu, GR_tr_guess)
    
    G_each_K<-append(G_each_K, G)
    E_each_K<-append(E_each_K, E)
    
    result_validate<-Validate(G, E, X_va, Y_va, GR_va)
    
    validation_score<-append(validation_score, result_validate[[1]])
    
    val_result$Lab<-append(val_result$Lab, GR_va)
    val_result$Simu<-append(val_result$Simu, result_validate[[2]])
    
    cat("Test:\n", "Lab Growth Rates:\n")
    cat("Test:\n", "Lab Growth Rates:\n", file="run_report.txt", append=T)
    print(GR_va)
    cat(GR_va, file="run_report.txt", append=T)
    cat("Dot Product Prediction of Growth Rate:\n")
    cat("Dot Product Prediction of Growth Rate:\n", file="run_report.txt", append=T)
    print(result_validate[[2]])
    cat(result_validate[[2]], file="run_report.txt", append=T)
    
  }
  Draw(K, tr_result, "training")
  Draw(K, val_result, "testing")
  
  return(list(score, validation_score, tr_result, val_result, G_each_K, E_each_K))
}