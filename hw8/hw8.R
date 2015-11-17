setwd("~/Desktop/school/bayesian-data-analysis/hw8")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

survey.resp <- read.table("http://www.stat.columbia.edu/~gelman/bda.course/incentives_data_clean.txt", skip = 12)

# indicator for incentive occurance
survey.resp$incentive <- ifelse(survey.resp$f==0,0,1)

idx <- data.frame(idx=1:39,sid=unique(survey.resp$sid))
survey.resp <- merge(survey.resp,idx)
survey.resp$one <- 1
survey.resp$r.dif <- NULL
survey.resp$v.dif <- NULL

stanData <- as.list(survey.resp)
stanData$X <- survey.resp[,c('one','incentive','m','b')]
stanData$N <- nrow(survey.resp)
stanData$M <- ncol(stanData$X)

#fit1 <- stan(file='model1.stan',data=stanData,iter=100)
#saveRDS(fit1,'fit1.rds')
fit1result <- stan(fit=fit1,data=stanData,iter=1000)


getInference <- function(stanData){
  fit1result <- stan(fit=fit1,data=stanData,iter=1000)
  results <- extract(fit1result)
  params <- cbind(results$betas,results$sigma,results$tau)
  colnames(params) <- c(names(stanData$X),'sigma','tau')
  
  means <- colMeans(params)
  
  intervals   <- do.call("rbind",(lapply(1:ncol(params),function(x) quantile(params[,x],c(.25,.75)))))
  halfWidth   <- (intervals[,2] - intervals[,1])/2
  return(100*data.frame(param=means,se=halfWidth))
}

trial1 <- getInference(stanData)


interaction1 <- stanData$X
interaction1 <- within(interaction1,{
  mb <- m*b
  incV <- incentive*survey.resp$v
  incT <- incentive*survey.resp$t
  incF <- incentive*survey.resp$f
  incM <- incentive*m
  incB <- incentive*b
})
  
stanData2 <- stanData
stanData2$X <- interaction1
stanData2$M <- ncol(interaction1)
trial2 <- getInference(stanData2)  


interaction2 <- within(interaction1,{
  incVT <- incV*survey.resp$t
  incVB <- incV*b
  incTB <- incT*b
  incVF <- incV*survey.resp$v
  incVM <- incV*m
  incTF <- incT*survey.resp$f
  incTM <- incT*m
  incFM <- incF*m
  incFB <- incF*b
  incMB <- incM*b
})

stanData3 <- stanData2
stanData3$X <- interaction2
stanData3$M <- ncol(interaction2)

trial3 <- getInference(stanData3)


clean3 <- data.frame(parameter=row.names(trial3),trial3=with(trial3, paste0(round(param,1)," (",round(se,1),")")),stringsAsFactors = F)
clean2 <- data.frame(parameter=row.names(trial2),trial2=with(trial2, paste0(round(param,1)," (",round(se,1),")")),stringsAsFactors = F)
clean1 <- data.frame(parameter=row.names(trial1),trial1=with(trial1, paste0(round(param,1)," (",round(se,1),")")),stringsAsFactors = F)
merged <- merge(clean3,merge(clean1,clean2,by="parameter",all.y=T),by="parameter",all.x=T)
merged[is.na(merged)]<-''

merged$parameter
merged$Parameter <- c("Burden","Incentive X Burden","Incentive","Incentive X Form","Incentive X Form X Burden","Incentive X Form X Mode",
  "Incentive X Mode","Incetive X Mode X Burden","Incetive X Time","Incetive X Time X Burden","Incentive X Time X Form",
  "Incentive X Time X Mode","Incentive X Value","Incentive X Value X Burden","Incentive X Value X Form",
  "Incentive X Value X Mode","Incentive X Value X Time","Mode","Mode X Burden","Constant","Sigma","Tau")

ordered <- merged[c(20,3,18,1,19,13,4,7,9,2,5,6,8,10,11,12,14,15,16,17,21,22),]
final <- ordered[,c("Parameter","trial1","trial2","trial3")]
saveRDS(final,"final.rds")

