setwd("~/Desktop/school/bayesian-data-analysis/hw9")

## q2: switches in binary data
outcome <- c(1, 1, 0, 0, 0, 0, 0,1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)

# posterior predictive check
simY <- function(){
  # simulate y
  theta <- rbeta(1,8,14)
  y <- c()
  while(sum(y==0)<13){
    y <- c(y,rbinom(1,1,theta))
  }
  # count number of switches
  T_ <- 0
  for(i in seq_along(y)){
    if(i==1){
      next
    }
    if(y[i]!=y[i-1]){
      T_ <- T_ + 1
    }
  }
  return(T_)
}

sims <- data.frame(T_=unlist(lapply(1:10000,function(i) simY())))
require(ggplot2)
switchesPlot <- ggplot(sims,aes(x=T_)) + 
  geom_histogram(binwidth=1,fill="#56B4E9") + 
  geom_vline(xintercept=3) +
  xlab("Number of Switches") + ylab("") +
  ggtitle("Posterior Predictive Check\nof Number of Switches in Binary Simulation")
saveRDS(switchesPlot,"switchesPlot.rds")

freqs <- as.data.frame(table(sims$T_))
sum(freqs[as.numeric(freqs$Var1)<=3,2])/sum(freqs$Freq)
