setwd('~/Desktop/school/bayesian-data-analysis/hw2')

library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data <- read.table("http://www.stat.columbia.edu/~gelman/bda.course/electric.txt",skip=2)
names(data) <- c("city","grade","ePre","ePost","cPre","cPost","replace.suppliment")

# compare discriptive stats
library(dplyr)
summarized <- data %>% group_by(grade) %>% summarise(e.premean=mean(e.pretest),
                                       e.postmean=mean(e.posttest),
                                       c.premean=mean(c.pretest),
                                       c.postmean=mean(c.posttest))

summarized <- within(summarized,{e.diff <-e.postmean-e.premean
                 c.diff <-c.postmean-c.premean})



# this takes a subset of the original data and simulates the stan model
runStan <- function(df,plotTitle){

  # unpack df columns
  for(i in 1:ncol(df)){
    assign(colnames(df)[i],df[,i])
  }
  
  N <- nrow(df)
  
  fit <- stan(file="model1_try2.stan",
              data=c("cPost","cPre","ePost","ePre","N"),chains=4,iter=100)
  saveRDS(fit,paste0(gsub(" ","",plotTitle),".rds"))

  paramsTry2 <- extract(fit)
  
  tdiff <- paramsTry2$tdiff
  dim(tdiff)<-NULL
  names(tdiff) <- gsub(" ","",plotTitle)
  
  plot <- ggplot(data=as.data.frame(tdiff)) + 
    geom_histogram(binwidth=1,aes(x=tdiff)) +
    ggtitle(plotTitle)
  
  stats <- list(interval=quantile(paramsTry2$tdiff,c(.25,.975)),
                mean=mean(paramsTry2$tdiff),
                plot=plot,
                data=tdiff)
  return(stats)
}

# test on everything
stats <- runStan(data,"All Students")

# subset interesting segements
grades        <- lapply(unique(data$grade), function(x) subset(data,grade==x))
cities        <- lapply(unique(data$city ), function(x) subset(data,city==x))
names(grades) <- paste0("grade",unique(data$grade))
names(cities) <- paste0("city", unique(data$city))

# run stan model on all subsets
allSubsets    <- c(grades,cities)
allResults    <- lapply(seq_along(allSubsets),function(i) runStan(allSubsets[[i]],names(allSubsets)[i]))
names(allResults) <- names(allSubsets)

# plot distribution of all posteriors
dataResults <- lapply(seq_along(allResults),function(i) 
  data.frame(subset=names(allResults)[i],
             tdiff=allResults[[i]]$data) )
dataResults <- do.call('rbind',dataResults)

ggplot(data=dataResults,aes(tdiff)) + 
  geom_histogram(binwidth=1) +
  facet_wrap(~ subset) +
  ggtitle("Effect of Test Controlling for Dimensions")

allResults[[1]]$interval[1]
# all posterior means & 95% quantiles
summaryStats <- lapply(seq_along(allResults),function(i) 
  data.frame(subset=names(allResults)[i],
             mean=allResults[[i]]$mean,
             int025=allResults[[i]]$interval[1],
             int975=allResults[[i]]$interval[2]) )
summaryStats <- do.call('rbind',summaryStats)
row.names(summaryStats) <- NULL

#########
# part 2
#########
library(ggplot2)

# make data for logistic regression from Table 3.1
x    <- c(-0.86,-0.30,-0.05,0.73)
n    <- c( 5,5,5,5)
y    <- c( 0,1,3,5)

x_trials <- sort(rep(x,5))
y_trials <- unlist(lapply(y,function(x) c(rep(0,5-x),rep(1,x))))
J        <- length(x_trials)

# fit stan model to plot prior simulation
fit <- stan(file="model2pre.stan",data=c("x_trials","y_trials","J"),chains=4,iter=1000)
saveRDS(fit,'model2pre.rds')

# post processing
fitPrior             <- readRDS('model2pre.rds')
simsPrior            <- extract(fitPrior)
simsPlotPrior        <- as.data.frame(simsPrior$b_temp)
names(simsPlotPrior) <- c('alpha','beta')

# visualition of prior simulation
ggplot(data=simsPlotPrior,aes(x=alpha,y=beta)) + geom_point()
cor(simsPlotPrior$alpha,simsPlotPrior$beta)


# Tune posterior on logistic regression model
fit <- stan(file="model2post.stan",data=c("x_trials","y_trials","J"),chains=4,iter=1000)
saveRDS(fit,'model2Post.rds')

# post process
plot(fit)
sims <- extract(fit)
simsPlot <- as.data.frame(sims$b_prior)
names(simsPlot) <- c('alpha_post','beta_post')

# compute expected LD50
simsPlot <- within(simsPlot,LD50 <- -alpha_post/beta_post)

# posterior distribution plot of LD50
# note: similar to results found in book
ggplot(data=simsPlot,aes(x=LD50)) + geom_histogram(binwidth=.02) +xlim(-.5,.4)
mean(simsPlot$LD50)

# posterior distribution of beta.
# note: all > 0 indicating increase in dosage -> increase in death likelihood
ggplot(data=simsPlot,aes(x=beta_post)) + geom_histogram()
