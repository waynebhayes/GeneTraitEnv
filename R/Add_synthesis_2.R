Add_synthesis_2<-function(DB_wild, num){
  
  DB_sim<-DB_wild[, c("Temperature", "Carbon_Glycerol",	"Carbon_Glucose",	"Carbon_L.Lactate", "Growth_Rate")]
  DB_sim[is.na(DB_sim)]<-0
  DB_remain<-DB_wild[, col_gene]
  DB_sim$Fluc<-rep(0, nrow(DB_sim))
  DB_sim$GR_synthesis<-rep(0, nrow(DB_sim))
  DB_sim$GR_synthesis[1:nrow(DB_wild)]<-DB_wild[,"Growth_Rate"]
  
  
  DB_syn<-DB_sim

  
  
  carbon_source<-c("Carbon_Glycerol",	"Carbon_Glucose",	"Carbon_L.Lactate")
  carbon_eff<-c(1,1.17,0.66)
  names(carbon_eff)<-carbon_source
  
  for (i in 1:nrow(DB_sim)){
    carbon_old<-carbon_source[which(DB_sim[i, carbon_source]>0)]
    data_temp<-DB_sim[rep(i,num), ]
    T_base<-DB_sim[i,"Temperature"]
    conc_base<-DB_sim[i,carbon_old]
    GR_real<-DB_sim[i,"Growth_Rate"]
    
    
    for (j in 1:num){
      carbon_new<-carbon_source[sample(1:length(carbon_source),1)]
      data_temp[j, carbon_old]<-0
      data_temp[j, carbon_new]<-conc_base*runif(1, 0.7, 1.3)
      data_temp[j, "Temperature"]<-data_temp[j, "Temperature"]*runif(1, 0.7, 1.3)
      data_temp[j, "Fluc"]<-0.005*rnorm(1)
      data_temp[j, "GR_synthesis"]<-GR_real/conc_base*carbon_eff[carbon_old]/carbon_eff[carbon_old]*data_temp[j, carbon_new]/T_base^2*data_temp[j, "Temperature"]^2+data_temp[j, "Fluc"]
      data_temp[j, "Growth_Rate"]<-0
      if (data_temp[j, "GR_synthesis"]>1.5)
        cat("Original rate:", GR_real, "  Synthesized rate:", data_temp[j, "GR_synthesis"], "\n")
    }
    
    DB_syn<-rbind(DB_syn, data_temp)
    
  }
  
  DB_syn<-cbind(DB_syn, DB_remain)
  
  return(DB_syn)
  
  
}