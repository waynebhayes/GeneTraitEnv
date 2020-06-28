Validate<-function(G, E, X_va, Y_va, GR_va){
  TG_va<-G%*%X_va
  TE_va<-E%*%Y_va
  
  GR_guess_va<-rep(0, length(GR_va))
  
  for (i in 1:length(GR_guess_va)){
    # GR_guess[i]<-(TG[,i]+L[[3]])%*%TE[,i]
    GR_guess_va[i]<-TG_va[,i]%*%TE_va[,i]
  }
  
  rho_sp <- cor(GR_va, GR_guess_va, method = "spearman")
  rho_pr <- cor(GR_va, GR_guess_va, method = "pearson")
  
  GR_mean_diff<-mean(abs(GR_guess_va-GR_va))
  
  rho <- (0*rho_sp + 2*rho_pr)/2

  obj<-abs(GR_mean_diff*(1-rho))
  
  return(list(obj, GR_guess_va))
  
}