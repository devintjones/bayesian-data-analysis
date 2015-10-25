setwd('~/Desktop/school/bayesian-data-analysis/hw3')

library(rstan)
library(ggplot2)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

y <- c(-2, -1, 0, 1.5, 2.5)
data1 <- list(y=y,N=length(y))
fit1 <- stan(file='hw3.1.stan',data=data1,iter=2000)
plot(fit1)

# 1.a) plot the posterior density
curve(dcauchy(x)/(.01*sum(dcauchy(x))),0,1,ylim=c(0,1.3),ylab="P(theta|y)",xlab="theta")

curve(log(dcauchy(x)/(.01*sum(dcauchy(x)))),0,1,ylab="P(theta|y)",xlab="theta")

# 1.d)
curve(dnorm(x,0,(1-mean(y)**2)/2),0,1,ylim=c(0,1),ylab="P(theta|y",xlab="Theta")


# Q2

chickens <- read.table('http://www.stat.columbia.edu/~gelman/bda.course/chickens.txt',header=T,skip=1)
names(chickens) <- c('freq','N_s','mean_s','se_s','N_e','mean_e','se_e')
head(chickens)
x
library(ggplot2)
ggplot(data=chickens,aes(x=freq,y=mean_e-mean_s)) + geom_point() + stat_smooth()

data2   <- as.list(chickens)
data2$N <- nrow(chickens)
model2t  <- stan('hw3.2.t.stan',data=data2,iter=2000)
print(model2t)
plot(model2t)

freq  <- extract(model2t)$freq_effect
sham  <- extract(model2t)$sham_effect
ratio <- extract(model2t)$ratio
freqStats <- data.frame(
  freqMeans = apply(freq,2,mean),
  freqLower = apply(freq,2,function(x) quantile(x,prob=.025)),
  freqUpper = apply(freq,2,function(x) quantile(x,prob=.975)),
  shamMeans = apply(sham,2,mean),
  shamLower = apply(sham,2,function(x) quantile(x,prob=.025)),
  shamUpper = apply(sham,2,function(x) quantile(x,prob=.975)),
  ratioMeans = apply(ratio,2,mean),
  ratioLower = apply(ratio,2,function(x) quantile(x,prob=.025)),
  ratioUpper = apply(ratio,2,function(x) quantile(x,prob=.975))
  )
freqStats <- within(freqStats, actualRatio<- chickens$mean_e/chickens$mean_s)
saveRDS(freqStats,"freqStats.rds")

ratioPlot <- ggplot(data=freqStats,aes(x=1:nrow(freqStats),y=actualRatio)) + 
  geom_errorbar(aes(ymin=ratioLower,ymax=ratioUpper),width=.1) +
  geom_point() +
  xlab("Observation") +
  ylab("Ratio: Observed & estimated interval") +
  ggtitle("Ratio of Exposed ratio\n to Shame Ratio")
saveRDS(ratioPlot,"ratioPlot.rds")

meanRatio <- data.frame(meanRatio=extract(model2t)$mean_ratio) 
meanPlot <- ggplot(meanRatio,aes(x=meanRatio)) + 
  geom_histogram() +
  xlab("Mean Effect of Treatment") +
  ylab("") +
  ggtitle("Posterior of the Mean Effect of Treatment")
saveRDS(meanPlot,"meanPlot.rds")


