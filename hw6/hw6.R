setwd("~/Desktop/school/bayesian-data-analysis/hw6")
library(boot)

## data
x.init   <- c(-0.86,-0.30,-0.05,0.73)
n        <- c( 5,5,5,5)
y.init   <- c( 0,1,3,5)
x        <- sort(rep(x.init,5))
y        <- unlist(lapply(y.init,function(x) c(rep(0,5-x),rep(1,x))))
J        <- length(x)

# log probability of y given params
log.prob.y <- function(x,y,alpha,beta){
    sum(-log(1+exp(alpha+x*beta))) + sum(y*(alpha+x*beta))
}

# derivative of log prob wrt alpha
d.alpha.log.prob.y <- function(x,y,alpha,beta){
  sum(y-inv.logit(alpha+beta*x))
}

# derivative of log prob wrt beta
d.beta.log.prob.y <- function(x,y,alpha,beta){
  sum((y-inv.logit(alpha+beta*x))*x)
}

# normal. used for alpha, beta
log.prob.normal <- function(mu.star,mu.prior,sigma.prior){
  log(dnorm(mu.star,mu.prior,sigma.prior))
}
d.mu.log.prob.normal <- function(mu.star,mu.prior,sigma.prior){
  (1/(sigma.prior**2))*(mu.prior-mu.star)
}
d.sigma.log.prob.normal <- function(mu.star,mu.prior,sigma.prior){
  -1/(sigma.prior)+(1/(sigma.prior**3))*(mu.star-mu.prior)**2
}

mu.star <- .5
mu.prior <- 0
sigma.prior <-.5
# test normal drivative wrt mean
(log.prob.normal(mu.star,mu.prior,sigma.prior) - log.prob.normal(mu.star,mu.prior+eps,sigma.prior))/eps
d.mu.log.prob.normal(mu.star,mu.prior+eps/2,sigma.prior)

# test normal derivative wrt sigma
(log.prob.normal(mu.star,mu.prior,sigma.prior) - log.prob.normal(mu.star,mu.prior,sqrt(sigma.prior**2-eps)))/eps
d.sigma.log.prob.normal(mu.star,mu.prior,sqrt(sigma.prior**2-eps/2))


# test that derivative formulas are correct
alpha <- .5
beta <- 1
eps <- 0.0001

## test alpha
numerical.alpha  <- (log.prob.y(x,y,alpha,beta) - log.prob.y(x,y,alpha-eps,beta))/eps
analytical.alpha <- d.alpha.log.prob.y(x,y,alpha-eps/2,beta)
if(abs(numerical.alpha-analytical.alpha)>.001) stop("alpha derivative is wrong")

## test beta
numerical.beta  <- (log.prob.y(x,y,alpha,beta) - log.prob.y(x,y,alpha,beta-eps))/eps
analytical.beta <- d.beta.log.prob.y(x,y,alpha,beta-eps/2)
if(abs(numerical.beta-analytical.beta)>.001) stop("beta derivative is wrong")


# total posterior
alpha.prior.mu    <- .8
alpha.prior.sigma <- 1
beta.prior.mu     <- 7.7
beta.prior.sigma  <- 4.9

posterior <- function(x,y,alpha,beta){
  log.prob.y(x,y,alpha,beta) + 
    log.prob.normal(alpha,alpha.prior.mu,alpha.prior.sigma) +
    log.prob.normal(beta,beta.prior.mu,beta.prior.sigma)
}

# total posterior gradient wrt alpha
d.post.alpha <-function(x,y,alpha,beta){
  d.alpha.log.prob.y(x,y,alpha,beta) +
    d.mu.log.prob.normal(alpha,alpha.prior.mu,alpha.prior.sigma)
}

# test
(posterior(x,y,alpha,beta) - posterior(x,y,alpha-eps,beta))/eps
d.post.alpha(x,y,alpha-eps/2,beta)

# total posterior gradient wrt beta
d.post.beta <-function(x,y,alpha,beta){
  d.beta.log.prob.y(x,y,alpha,beta) +
    d.mu.log.prob.normal(beta,beta.prior.mu,beta.prior.sigma)
}

# test
(posterior(x,y,alpha,beta) - posterior(x,y,alpha,beta-eps))/eps
d.post.beta(x,y,alpha,beta-eps/2)


# gradient used in HMC
d.log.p.theta <- function(x,y,alpha,beta){
  c(d.post.alpha(x,y,alpha,beta),d.post.beta(x,y,alpha,beta))
}



getDiag <- function(mat){
  unlist(lapply(1:ncol(mat), function(i) mat[i,i]))
}

# init M and phi
M <- diag(c(alpha.prior.sigma,beta.prior.sigma))
#M <- diag(c(1,2))
init.phi <- function(M){
  unlist(lapply(1:ncol(M),function(i) rnorm(1,0,M[i,i])))
}
init.theta <- function(alpha.prior.mu,alpha.prior.sigma,beta.prior.mu,beta.prior.sigma){
  c(runif(1,alpha.prior.mu-.5,alpha.prior.mu +.5),runif(1,beta.prior.mu-2,beta.prior.mu+2))
}

post.theta.phi <- function(x,y,theta,phi,vars){
  posterior(x,y,theta[1],theta[2]) + 
    log(dnorm(phi[1],0,vars[1])) + 
    log(dnorm(phi[2],0,vars[2]))
}


HMCiter <- function(theta,M,e,L){
  theta.init <- theta
  phi        <- phi.init <- init.phi(M)
  vars       <- getDiag(M)
  
  i <- 1
  while(i<=L){
    phi   <- phi   + .5*e*d.log.p.theta(x,y,theta[1],theta[2])
    theta <- theta + e*phi/vars
    phi   <- phi   + .5*e*d.log.p.theta(x,y,theta[1],theta[2])
    i <- i+1
  }

  # diff of log posterior probs
  r <- exp(post.theta.phi(x,y,theta,phi,vars) -
             post.theta.phi(x,y,theta.init,phi.init,vars))
  if(min(r,1) > runif(1)){
    theta.new <- theta
    accept <- 1
  }else{
    theta.new <- theta.init
    accept <- 0
  }
  return(list(theta=theta.new,r=r,accept=accept,phi=phi))  
}

HMC.chain <- function(M,e,L,iter=1000){
  theta <- init.theta(alpha.prior.mu,alpha.prior.sigma,beta.prior.mu,beta.prior.sigma)
  thetas <- c()
  accepts <- 0
  for(i in 1:iter){
    e.star <- runif(1,0,2*e)
    L.star <- runif(1,1,round(2*L))
    params <- HMCiter(theta,M,e.star,L.star)
    theta  <- params$theta
    
    if(i > iter/2){
      thetas <- rbind(thetas,theta)
      accepts <- accepts + params$accept
    }
  }
  acceptance <- accepts/(iter/2)
  return(list(thetas=thetas,avg.acceptance=acceptance))
}

e   <- .85
L   <- round(1/e)

M <- diag(x=c(2,10))
alpha.prior.sigma <- 1
beta.prior.sigma  <- 10

library(parallel)
results <- mclapply(1:4,function(x) HMC.chain(M,e,L,iter=1500))
printResults(results)
plotResults(results)

#saveRDS(results,"resultsMany.rds")

printResults <- function(results){
  cat('rhat for alpha:',rhat(results,1),
      '\nrhat for beta:',rhat(results,2),
      '\nneff alpha:',neff(results,1),
      '\nneff beta:',neff(results,2),
      '\nacceptance rate:',paste(lapply(1:4,function(i) results[[i]]$avg.acceptance),collapse=' '))
}

plotResults <- function(results){
  plotData <- as.data.frame(do.call('rbind',lapply(1:length(results), function(i) results[[i]]$thetas)))
  names(plotData) <- c('alpha','beta')
  require(reshape2)
  plotData <- melt(plotData)
  head(plotData)
  ggplot(plotData,aes(x=value)) +
    geom_histogram()+
    facet_wrap(~variable,nrow = 2) +
    ggtitle("Posterior Distribution of alpha, beta")
}


# compute rhat
rhat <- function(results,beta.idx=1){
  m <- length(results)
  n <- nrow(results[[1]]$thetas)
  chain.sims <- function(i) results[[i]]$thetas[,beta.idx]
  
  mean.sims  <- mean(unlist(lapply(1:m,function(i) mean(chain.sims(i)))))
  
  # deviance of mean of chain from mean of means of chains
  B <- n/(m-1)*sum(unlist(lapply(1:m, function(i) (mean(chain.sims(i)) - mean.sims)**2)))
  
  # mean standard deviation of chains
  W <- mean(unlist(lapply(1:m,function(i) sd(chain.sims(i)))))
  
  var.phi <- (n-1)/n*W + (1/n)*B
  
  rhat <- sqrt(var.phi/W)
  return(rhat)
}


# for use in autocorelation calcs
lag_ <- function(x,k=0) c(rep(NA, k), x)[1 : length(x)] 

# effective sample size for a cahing given m chains, n iters
neff <- function(results,beta.idx){
  m <- length(results)
  n <- nrow(results[[1]]$thetas)
  chain.sims <- function(i) results[[i]]$thetas[,beta.idx]
  
  mean.sims  <- mean(unlist(lapply(1:m,function(i) mean(chain.sims(i)))))
  
  # deviance of mean of chain from mean of means of chains
  B <- n/(m-1)*sum(unlist(lapply(1:m, function(i) (mean(chain.sims(i)) - mean.sims)**2)))
  
  # mean standard deviation of chains
  W <- mean(unlist(lapply(1:m,function(i) sd(chain.sims(i)))))
  
  # var beta|y
  var.psi <- (n-1)/n*W + (1/n)*B
  
  # variogram t
  v.t <- function(t){
    (1/(m*(n-1)))*sum(unlist(lapply(1:m, function(m){ 
      sum(((lag_(chain.sims(m),t)-chain.sims(m))**2)[(t+1):n])
    })))
  }
  
  # correlation at t
  rho.t <- function(t) 1 - (v.t(t)/(2*var.psi))
  
  # correaltions over all t
  rhos <- unlist(lapply(1:(n-1), function(t) rho.t(t)))
  
  # effective sample size
  neff <- m*n/(1 + 2*sum(rhos[1:which.max(lag_(rhos)+rhos<0)]))
  
  return(neff)
}
