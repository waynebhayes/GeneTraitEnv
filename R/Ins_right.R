Ins_right<-function(DB, Num, NewName){
  
  New<-DB[, Num]
  
  if (Num < ncol(DB)){
    table_secondhalf<-DB[,(Num+1):ncol(DB)]
    name_firsthalf<-colnames(DB)[1:Num]
    name_secondhalf<-colnames(DB)[(Num+1):ncol(DB)]
    
    DB<-cbind(DB, New)
    
    DB[, Num+1]<-New
    DB[, (Num+2):ncol(DB)]<-table_secondhalf
    colnames(DB)<-c(name_firsthalf, NewName, name_secondhalf)
  }
  
  else {
    name_new<-c(colnames(DB), NewName)
    DB<-cbind(DB, New)
    colnames(DB)<-name_new
    }
  return(DB)
}