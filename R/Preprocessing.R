Preprocessing<-function(input, synthesis = F){
  
  if (is.character(input)){
    DB<-read.csv(input, stringsAsFactors = FALSE)[,-1]  
  } else {
    DB<-input
  }
  
  
  colnames_total<-colnames(DB)
  col_gene<-(grep("Gene_start", colnames_total)+1):(grep("Growth_Rate", colnames_total)[1]-1)
  col_gene_start<-grep("Gene_start", colnames_total)+1
  
  if (synthesis==T){
    DB_wild<-data.frame()
    
    for (i in 1:nrow(DB)){
      if (sum(DB[i, col_gene])==0){
        DB_wild<-rbind(DB_wild, DB[i,])
      }
    }
    
    return(DB_wild)
  }
  
  i<-col_gene_start
  empty_col<-rep(0, nrow(DB))
  while(i<grep("Growth_Rate", colnames(DB))[1]){
    
    if (grepl("_Deletion", colnames(DB)[i])){
      underscore_pos<-regexpr("_Deletion", colnames(DB)[i])
      attributes(underscore_pos)<-NULL
      gene_name<-substr(colnames(DB)[i], 1, underscore_pos-1)
      gene_all_con<-grep(gene_name, colnames(DB))
      
      new_col<-empty_col
      
      for (j in 1:nrow(DB)){
        if (sum(DB[j, gene_all_con])==0){
          new_col[j]<-1
        }
      }
      
      DB<-Ins_right(DB, i, paste(gene_name, "_Normal", sep=""))
      DB[ ,paste(gene_name, "_Normal", sep="")]<-new_col
      
      if (sum(DB[, i])==0){
        DB<-Delete(DB, colnames(DB)[i])
        i<-i+1
      }
      else{
        i<-i+2
      }
      next;
    }
    
    if (sum(DB[, i])==0){
      DB<-Delete(DB, colnames(DB)[i])
    }
    else{
      i<-i+1
    }
  }
  
  if (is.character(input)){
    write.csv(DB, "Data_Matrix_Reduced.csv")  
  } else {
    write.csv(DB, "Data_Matrix_Reduced_wild.csv")
  }
  
  return(DB)
}