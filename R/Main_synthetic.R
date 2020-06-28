Main_synthetic<-function(input, basal, seed, traits, K, obj_best){
  
  DB_wild<-Preprocessing(input, synthesis=T)
  DB_wild<-Preprocessing(DB_wild, synthesis=F)
  
  
  val_result<-list("Lab"=c(), "Simu"=c())
  tr_result<-list("Lab"=c(), "Simu"=c())
  G_each_K<-c()
  E_each_K<-c()
  
  colnames_total<<-colnames(DB_wild)
  col_env<<-c(colnames_total[(grep("Env_start", colnames_total)+1):
                         (grep("Gene_start", colnames_total)[1]-1)])
  col_gene<<-c(colnames_total[(grep("Gene_start", colnames_total)+1):
                          (grep("Growth_Rate", colnames_total)-1)])
  col_gene_start<<-grep("Gene_start", colnames_total)+1
  
  
  
  sample_ID<-DB_wild[, "sample_ID"]
  
  # DB_syn<<-Add_synthesis(DB_wild, "env", 20)
  DB_syn<-Add_synthesis_2(DB_wild, 20)
  col_env<-c("Temperature", "Carbon_Glycerol",	"Carbon_Glucose",	"Carbon_L.Lactate", "Growth_Rate")
  
  
  
  if ("X.1" %in% colnames(DB_syn)){
     DB_syn<<-DB_syn[,-c("X.1")]
   }
   
  
  # DB_syn<<-DB_syn[,-c("Growth_Rate", "Basal_Growth_Rate")]
  
  # GR_syn<-DB_syn[, GR_synthesis]
  GR_syn<-DB_syn[, "GR_synthesis"]
  shuffle_idx<-sample(1:nrow(DB_syn))
  

  
  
  
  X<-t(DB_syn[,col_gene])
  # X<-t(DB_syn[shuffle_idx,col_gene,with=F])
  X[is.na(X)]<-0
  Y<-t(DB_syn[,col_env[which(col_env!="Growth_Rate")]])
  # Y<-t(DB_syn[shuffle_idx,col_env,with=F])
  Y[is.na(Y)]<-0
  
  result<-list()
  norm_best<-0.03
  rho_best<-0.8
  
  #print(tr_whichrow)
  
  score<-vector("numeric")
  validation_score<-vector("numeric")
  
  GR_va_guess<-c()
  GR_va_all<-c()
  result_train<- vector(mode = "list", length = K)
  num_exp<-nrow(DB_syn)
  # num_exp<-nrow(DB_syn)
  round<-floor(num_exp/K)
  cat("Synthetic experiment.\n", "seed=", seed, "traits=", traits, "\n", file="run_report.txt")
  cat("Synthetic experiment.\n", "seed=", seed, "traits=", traits, "\n", file="run_result.txt")
  cat("Synthetic experiment.\n", "seed=", seed, "traits=", traits, "\n", file="run_result_top.txt")
  
  write.csv(DB_syn, "Data_Synthesized.csv")
  
  for (i in 1:K){
    
    G_initial<-matrix(rnorm(traits*nrow(X)), traits, nrow(X))
    E_initial<-matrix(rnorm(traits*nrow(Y)), traits, nrow(Y))
    
    G<-G_initial
    E<-E_initial
    cat("Initial G:\n", G, "\n", file="run_result.txt", append=T)
    cat("Initial E:\n", E, "\n", file="run_result.txt", append=T)
    
  
    
    if (i!=K){
      va_index<-shuffle_idx[(round*(i-1)+1):(round*i)]  
    } else{
      va_index<-shuffle_idx[(round*(i-1)+1):num_exp]
    }
    
    tr_index<-setdiff(shuffle_idx, va_index)
    cat("K=", i, "\n", file="run_result.txt", append=T)
    cat("K=", i, "\n", file="run_result_top.txt", append=T)
    
    
    X_tr<-X[, tr_index]
    Y_tr<-Y[, tr_index]
    X_va<-X[, va_index]
    Y_va<-Y[, va_index]
    GR_tr<-GR_syn[tr_index]
    GR_va<-GR_syn[va_index]
    
    # return(list(X_tr, Y_tr, GR_tr))
    
    SA_Tinit<-1
    threshold<-0.5
    
    num_cycle<<-1
    num_round<<-1
    
    

    cat("Round ", num_round, "\n", "Goal of obj ", obj_best, "\n", file="run_report.txt", append=T)
    cat("Round ", num_round, "\n", "Goal of obj ", obj_best, "\n")
    
    
              
    while(GetpBad_synthetic(G, E, X_tr, Y_tr, GR_tr, SA_Tinit, threshold)<threshold){
      print(GetpBad(G, E, X_tr, Y_tr, GR_tr, SA_Tinit, threshold))
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
    cat("K=", i, "\n", file="run_report.txt", append=T)
    
    
    result_train[[i]]<-Train_SA_synthetic(num_exp, G, E, X_tr, Y_tr, GR_tr, obj_best, T_range)
    
    score<-append(score, result_train[[i]][[1]])
    
    G<-result_train[[i]][[4]][[1]]
    cat("G=", t(G), "\n", file="run_result.txt", append=T)
    
    
    
    E<-result_train[[i]][[4]][[2]]
    cat("E=", t(E), "\n", file="run_result.txt", append=T)
    
    GR_tr_guess<-result_train[[i]][[5]]
    
    tr_result$Lab<-append(tr_result$Lab, GR_tr)
    tr_result$Simu<-append(tr_result$Simu, GR_tr_guess)
    
    G_each_K<-append(G_each_K, G)
    E_each_K<-append(E_each_K, E)
    
    result_validate<-Validate(G, E, X_va, Y_va, GR_va)
    
    validation_score<-append(validation_score, result_validate[[1]])
    
    val_result$Lab<-append(val_result$Lab, GR_va)
    val_result$Simu<-append(val_result$Simu, result_validate[[2]])
    
    print(GR_va)
    print(result_validate[[2]])
    
  }
  Draw(K, tr_result, "training")
  Draw(K, val_result, "testing")
  
  return(list(score, validation_score, tr_result, val_result, G_each_K, E_each_K))
}