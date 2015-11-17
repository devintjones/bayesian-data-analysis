setwd("~/Desktop/school/bayesian-data-analysis/hw7")
library(rstan)
library(ggplot2)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


survey <- read.csv("http://www.stat.columbia.edu/~gelman/bda.course/naes04.csv")
survey$X <- NULL

# bucket ages into 5 year bands
survey <- within(survey, ageGroups <- 5*round(survey$age/5))

modelData <- survey[,c('ageGroups','gender','gayKnowSomeone')]
modelData <- modelData[complete.cases(modelData),]
modelData <- within(modelData,{
  gender    <- ifelse(modelData$gender=='Male',1,0)
})

aggData <- modelData %>% 
  group_by(ageGroups,gender) %>%
  summarize(numKnown=sum(gayKnowSomeone=='Yes'),
            total=length(gender),
            prop=numKnown/total)

library(reshape2)
ageGroups <- model.matrix(~as.factor(ageGroups),aggData)
ageGroups <- ageGroups[,-1]

gender    <- ifelse(aggData$gender=="Male",1,0)

forStan <- as.list(aggData)
forStan$ageGroups <- ageGroups
forStan$gender <- gender
forStan$N <- nrow(aggData)
forStan$M <- ncol(ageGroups)
forStan$age <- aggData$ageGroups
fit1 <- stan('gp.stan',data=forStan,iter=100)
saveRDS(fit1,'fit1.rds')
#fit1 <- stan(fit1,data=forStan,iter=100)
probs <- extract(fit1)$prob
print(fit1)
plot(fit1)
head(aggData)

# GP of % people who know a gay person by age and gender
# model: y(x) = f_age(x) + f_gender(x) + eps

# f_age(x) ~ GP(0,k1)  
k1 <- function(x,x_prime,tau=1,l=1){
  tau**2 * exp(-sum((x-x_prime)**2)/l**2)
}

# f_gender(x) ~ GP(0,k1)
k2 <- function(x,x_prime,tau=1,l=1){
  tau**2 * exp(-(x-x_prime)**2/l**2)
}

# compute covariance matrix over aggregated data
makeK <- function(ageGroups,gender,sigma2=.1){
  N <- length(gender)
  K <- matrix(nrow=N,ncol=N)
  for(i in 1:N){
    for(j in 1:N){
      sigma[i,j] <- k1(ageGroups[i,],ageGroups[j,]) + k2(gender[i],gender[j]) + ifelse(i==j,.1,0)
    }
  }
  K
}
y <- forStan$numKnown
N <- forStan$total
theta <- rnorm(length(y))
log.post <- function(theta,y,N,ageGroups,gender){
  # log binomial
  sum(log(mapply(dbinom,x=y,size=N,prob=inv.logit(theta)))) +
    -.5*t(theta)%*%solve(K)%*%theta
}

# plot stan results against data

estMean <- colMeans(probs)
intervals   <- do.call("rbind",(lapply(1:ncol(probs),function(x) quantile(probs[,x],c(.025,.975)))))
colnames(intervals) <- c("lower","upper")

plotData <- within(aggData,group<- paste0(gender,ageGroups))
plotData <- cbind(prePlot,intervals,estMean)

ggplot(data=plotData,aes(x=group,y=prop)) + 
  geom_point(color="#56B4E9",size=2) +
  geom_point(aes(x=group,y=estMean),color="black") +
  geom_errorbar(aes(ymax=upper,ymin=lower),width=.1)

head(prePlot)

subset(prePlot,prop>upper | prop<lower)



