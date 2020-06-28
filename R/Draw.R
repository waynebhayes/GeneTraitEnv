Draw<-function(K, result_input, keyword){
  
  progress<-c("0%", paste(round(100*(1:K)/K), "%", sep=""))
  
  filename<-paste(keyword, ".jpeg", sep="")
  
  jpeg(file=filename, width=240*2*K, height=240*2*K)
  layout(matrix(c(1:K, 1:K),K,2))
  step<-length(result_input[[1]])/K
  j<-0
  
  for (i in 1:K){
    
    lab_data<-result_input$Lab[((i-1)*step+1):(i*step)]
    simu_data<-result_input$Simu[((i-1)*step+1):(i*step)]
    
    order_data<-order(lab_data)
    lab_data_ordered<-lab_data[order_data]
    simu_data_rearr<-simu_data[order_data]
    
    plot(lab_data_ordered, main=paste(progress[i], progress[i+1], sep="-"), ylim=c(min(c(lab_data, simu_data)), max(c(lab_data, simu_data))), cex.main=K/2, cex.lab=K/2, cex.axis=K/2, cex=K/2)
    points(simu_data_rearr, col=2, cex=K/2)
    legend("bottomright", legend=c("Real Data", "Simulation"), pch=c(1,1), col=c("black", "red"))
  }
  
  dev.off()
  
  return()
}
