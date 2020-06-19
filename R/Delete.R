Delete<-function(DB, Colname){
  
  Num<-grep(Colname, colnames(DB))
  newtable_secondhalf<-DB[,(Num+1):ncol(DB)]
  newtable_firsthalf<-DB[,1:(Num-1)]
  
  name_secondhalf<-colnames(DB)[(Num+1):ncol(DB)]
  name_firsthalf<-colnames(DB)[1:(Num-1)]
  
  newDB<-cbind(newtable_firsthalf, newtable_secondhalf)
  colnames(newDB)<-c(name_firsthalf, name_secondhalf)
  
  return(newDB)
}