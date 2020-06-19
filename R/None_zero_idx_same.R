None_zero_idx_same<-function(m1, m2){
  
  
  product_multi<-m2*m1
  product_minus<-m2-m1
  
  result<-which(product_multi!=0, arr.ind=T)
  
  value_change<-rep(0, nrow(result))
  sign_change<-rep("?", nrow(result))
  
  for (i in 1:nrow(result)){
    idx<-as.vector(result[i, ])
    r<-idx[1]
    c<-idx[2]
    value_change[i]<-product_minus[r,c]
    sign_change[i]=case_when(
      (sign(m1[r,c])==-1 && sign(m2[r,c])==-1) ~ "-",
      (sign(m1[r,c])==1 && sign(m2[r,c])==1) ~ "+",
      (sign(m1[r,c])==-1 && sign(m2[r,c])==1) ~ "- -> +",
      (sign(m1[r,c])==1 && sign(m2[r,c])==-1) ~ "+ -> -"
    )
    if (is.na(sign_change[i])) {
      print(m1[idx])
      print(m2[idx])
    }
  }
  
  
  result<-cbind(result, sign_change, value_change)
  
  colnames(result)<-c("Row", "Col", "Sign Change", "Value Change")
  
  return(result)
}