Shuffle<-function(exp_index, paperdiv){
  
  exp_index_shuffle<-exp_index
  num_papers<-length(unique(exp_index[ ,2]))
    
  if (paperdiv){
    
    for (i in 1:num_papers){
      exp_paper<-matrix(exp_index[exp_index[ ,2]==i, ], ncol=2)
      exp_paper_shuffle<-exp_paper[sample(nrow(exp_paper)), ]
      exp_index_shuffle[exp_index_shuffle[ ,2]==i, ]<-exp_paper_shuffle
    }
  }
  
  else {
    exp_index_shuffle<-exp_index_shuffle[sample(exp_index_shuffle[ ,1]), ]
  }
  
  return (exp_index_shuffle)
}