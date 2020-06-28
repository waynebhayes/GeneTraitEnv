Kfold<-function(seed, K, paperdiv, sample_ID){
  
  set.seed(seed)
  
  num_exp<-length(sample_ID)
  
  dash_pos<-regexpr("-", sample_ID)
  attributes(dash_pos)<-NULL
  
  
  
  paper_ID<-rep(0, num_exp)
  for (i in 1:num_exp){
    paper_ID[i]<-substr(sample_ID[i], 1, dash_pos[i]-1)
  }
  
  papers<-unique(paper_ID)
  
  num_papers<-length(papers)
  
  from_paper<-rep(0, num_exp)
  
  for (i in 1:num_exp){
    from_paper[i]<-grep(paper_ID[i], papers)
  }
  
  exp_index<-cbind(1:num_exp, from_paper)

  exp_index_shuffle<-Shuffle(exp_index, paperdiv)
  
  num_papers<-length(unique(exp_index_shuffle[ ,2]))
  num_validation_samples<-vector("integer")
  
  tr_row_index<-list()
  va_row_index<-list()
  
  for (i in 1:K){
    tr_row_index[[i]]<-vector("integer")
    va_row_index[[i]]<-vector("integer")
  }
  
  
  if (paperdiv){
    rows_acc<-0
    for (i in 1:num_papers){
      all_exp<-matrix(exp_index[exp_index[,2]==i, ], ncol=2)
      num_selected<-floor(nrow(all_exp)/K)
      num_validation_samples<-append(num_validation_samples, num_selected)
      
      for (j in 1:K){
        if (num_selected!=0){
          row_index<-(num_selected*(j-1)+1):(num_selected*j)
          va_row_index[[j]]<-append(va_row_index[[j]], row_index+rows_acc)
        }
        else{
          row_index<-NULL
        }
        
        tr_row_index[[j]]<-append(tr_row_index[[j]], setdiff(1:nrow(all_exp), row_index)+rows_acc)
      }
      rows_acc<-rows_acc+nrow(all_exp)
    }
  }  
  else {
      for (j in 1:K){
        num_validation_samples<-floor(length(sample_ID/K))
        row_index<-(num_validation_samples*(j-1)+1):(num_validation_samples*j)
        va_row_index[[j]]<-append(va_row_index[[j]], row_index)
        tr_row_index[[j]]<-append(tr_row_index[[j]], setdiff(1:length(sample_ID), row_index))
      }
  }
  
  return (list(tr_row_index, va_row_index, exp_index_shuffle[ ,1]))
}