setwd('~/Desktop/school/bayesian-data-analysis/hw2')

library(rstan)
library(ggplot2)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data <- read.table("http://www.stat.columbia.edu/~gelman/bda.course/electric.txt",skip=2)
names(data) <- c("city","grade","ePre","ePost","cPre","cPost","replace.suppliment")




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
allResults    <- lapply(seq_along(allSubsets),function(i) 
  runStan(allSubsets[[i]],names(allSubsets)[i]))
names(allResults) <- names(allSubsets)

# add total average effect
allResults$allStudents <- stats

# plot distribution of all posteriors
dataResults <- lapply(seq_along(allResults),function(i) 
  data.frame(subset=names(allResults)[i],
             tdiff=allResults[[i]]$data) )
dataResults <- do.call('rbind',dataResults)
dataResults <- subset(dataResults,subset!='allStudents')

facetPlot <- ggplot(data=dataResults,aes(tdiff)) + 
  geom_histogram(binwidth=1) +
  facet_wrap(~ subset) +
  ylab('') +
  xlab('Simulated Effect of TV Viewer Program') +
  ggtitle("Effect of Test Controlling for Interesting Factors\nPosterior Simulation")
saveRDS(facetPlot,"facetPlot.rds")

# all posterior means & 95% quantiles
summaryStats <- lapply(seq_along(allResults),function(i) 
  data.frame(subset=names(allResults)[i],
             mean=allResults[[i]]$mean,
             int025=allResults[[i]]$interval[1],
             int975=allResults[[i]]$interval[2]) )
summaryStats <- do.call('rbind',summaryStats)
row.names(summaryStats) <- NULL
saveRDS(summaryStats,"summaryStats.rds")


# part b
# compare discriptive stats
summarizedGrade <- data %>% 
  mutate(ediff=ePost-ePre,
         cdiff=cPost-cPre,
         tdiff=ediff-cdiff) %>%
  group_by(grade,replace.suppliment) %>% 
  summarise(count=n(),
            tdiffmean = mean(tdiff),
            tdiffse   = sd(tdiff)/sqrt(n()),
            tdiffLower= tdiffmean + qt(.975,n()-1)*tdiffse,
            tdiffUpper= tdiffmean - qt(.975,n()-1)*tdiffse)
summarizedGrade

summarizedCity <- data %>% 
  mutate(ediff=ePost-ePre,
         cdiff=cPost-cPre,
         tdiff=ediff-cdiff) %>%
  group_by(city,replace.suppliment) %>% 
  summarise(count=n(),
            tdiffmean = mean(tdiff),
            tdiffse   = sd(tdiff)/sqrt(n()),
            tdiffLower= tdiffmean + qt(.975,n()-1)*tdiffse,
            tdiffUpper= tdiffmean - qt(.975,n()-1)*tdiffse)
summarizedCity

suppliment <-  data %>% 
  group_by(replace.suppliment) %>% 
  summarize(count=n(),
            tdiffmean = mean(tdiff),
            tdiffse   = sd(tdiff)/sqrt(n()),
            tdiffLower= tdiffmean + qt(.975,n()-1)*tdiffse,
            tdiffUpper= tdiffmean - qt(.975,n()-1)*tdiffse)
kable(suppliment)

##########
# part 2
#########

# make data for logistic regression from Table 3.1
x    <- c(-0.86,-0.30,-0.05,0.73)
n    <- c( 5,5,5,5)
y    <- c( 0,1,3,5)
x_trials <- sort(rep(x,5))
y_trials <- unlist(lapply(y,function(x) c(rep(0,5-x),rep(1,x))))
J        <- length(x_trials)

# fit stan model to plot prior simulation
fit <- stan(file="model2pre.stan",data=c("x_trials","y_trials","J"),chains=4,iter=1000)
saveRDS(fit,'model2.rds')

# post processing
fitPrior             <- readRDS('model2.rds')
simsPrior            <- extract(fitPrior)
simsPlotPrior        <- as.data.frame(simsPrior$b_temp)
names(simsPlotPrior) <- c('alpha','beta')

# visualition of prior simulation
priorPlot <- ggplot(data=simsPlotPrior,aes(x=alpha,y=beta)) + geom_point() + ggtitle("Prior simulation of\nLogistic Regression Parameters")
saveRDS(priorPlot,"priorPlot.rds")
cor(simsPlotPrior$alpha,simsPlotPrior$beta)


# Tune posterior on logistic regression model
setwd("~/Desktop/school/bayesian-data-analysis/hw2")
x_trials <- x
y_trials <- y
fit <- stan(file="model2post.stan",data=c("x_trials","y_trials","J"),chains=4,iter=1000)
saveRDS(fit,'model2_post.rds')
library(shinystan)
print(fit)
shinystan1 <- launch_shinystan(fit)
rhats <- stan_get(fit, what = "rhat")

# post process
fit <- readRDS('model2_post.rds')
sims <- extract(fit)
simsPlot <- as.data.frame(sims$b_prior)
names(simsPlot) <- c('alpha_post','beta_post')

# compute expected LD50
simsPlot <- within(simsPlot,LD50 <- -alpha_post/beta_post)

# posterior distribution plot of LD50
# note: similar to results found in book
ld50Plot <- ggplot(data=simsPlot,aes(x=LD50)) + 
  geom_histogram(binwidth=.02) +xlim(-.5,.4) + 
  ggtitle("Posterior Simulation\nof LD50") +
  ylab('') +
  xlab('LD50 Dosage (log scale)')
saveRDS(ld50Plot,'ld50Plot.rds')
mean(simsPlot$LD50)

# posterior distribution of beta.
# note: all > 0 indicating increase in dosage -> increase in death likelihood
betaPlot <- ggplot(data=simsPlot,aes(x=beta_post)) + 
  geom_histogram() + 
  ggtitle("Posterior Simulation\nBeta, the effect of Dosage") +
  ylab('') +
  xlab('Simulated Values of Beta\nThe inverse logit of the effect of log(dosage)\non probality of response')
saveRDS(betaPlot,"betaPlot.rds")


