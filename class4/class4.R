## fitting quadratic models

setwd('~/Desktop/school/bayesian-data-analysis/class4')

N <- 10
x <- 1:N
a <- 2
b <- -1
y <- a + b*x + rnorm(N,0,1)
data <- list(N=N,x=x,y=y)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit<- stan(file='class4.stan',data=data,iter=100)

a_hat <- mean(extract(fit)$a)
b_hat <- mean(extract(fit)$b)
sigma_hat <- mean(extract(fit)$sigma_y)

plot(x,y)
abline(a_hat,b_hat,col="red")
abline(a_hat+sigma_hat,b_hat,col="red",lty=2)
abline(a_hat-sigma_hat,b_hat,col="red",lty=2)
