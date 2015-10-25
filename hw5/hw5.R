setwd('~/Desktop/school/bayesian-data-analysis/hw5')




## 1.a
J <- 10
x <- runif(J,0,1)
n <- rpois(J,5)

deg.f <- 4

alpha <- rt(1,4)*2/2
beta  <- rt(1,4)*1/2

library(boot)
theta <- inv.logit(alpha + beta%*%x)

y <- rbinom(n,J,theta)


## 1.b

curve(dt(x,4),from=-10,to=10)
curve(dt.scaled(x,0,2,4),add=T,col="blue")

# add scale parameter to distirubtion function
dt.scaled <- function(x, m, s, df) dt((x-m)/s, df)/s

sample.post <- function(alpha.star,beta.star){
  p.alpha <- dt.scaled(alpha.star,0,2,4)
  p.beta  <- dt.scaled(alpha.star,0,1,4)
  theta   <- inv.logit( alpha.star +  beta.star * x)
  p.y     <- dbinom(y,J,theta)
  return (prod(p.y) * p.alpha * p.beta)
}

runs <- 200000
rejection.samps <- data.frame(alpha.star = runif(runs,-10,10),
                              beta.star  = runif(runs,-5,5))

# sample from uniform selection of alpha and beta
rejection.samps <- within(rejection.samps,prob <- mapply(sample.post,
                                                         alpha.star,beta.star))
# boolean
accepted<- with(rejection.samps, runif(runs,0,1)*max(prob,na.rm = T) < prob)

keepers <- rejection.samps[accepted,]
nrow(keepers)
nrow(rejection.samps)
plot(keepers$alpha.star,keepers$beta.star)

## 1.c

# function to bin and find mode
getMode <- function(RejectionSamples,cuts=100){
  # returns string: (a,b]
  modeRange <- names(sort(-table(cut(RejectionSamples,
                                     breaks=seq(min(RejectionSamples),
                                                max(RejectionSamples),
                                                length.out = cuts)
                                     )))[1])
  # clean string and find midpoint
  modeVal   <- sum(as.numeric(strsplit(gsub('\\(|]','',modeRange),',')[[1]]))/2
  return(modeVal)
}

# find mode and covar
alpha_mode <- getMode(keepers$alpha.star)
beta_mode  <- getMode(keepers$beta.star)
mu         <- c(alpha_mode,beta_mode)
covar      <- cov(keepers[,c('alpha.star','beta.star')])

# sample from bivariate normal
library(mnormt)
normalsim <- rmnorm(1000,mu,covar)
plot(normalsim)

## 1.d use importance sampling to estimate E(alpha|y) & E(beta|y)

### sample from multi variate t using covar matrix
tsim2 <- rmt(1000,mu,covar,4)
plot(tsim2)

### importance sampling
importance.sample <- function(vec){
  ## weight function: normal density
  w <- dunif(vec,min(vec),max(vec))/dnorm(vec,mean(vec),sd=(max(vec)-min(vec))/4)
  ## sampling function: samples from t distribution
  h <- vec
  return(mean(w*h)/mean(w))
}

importance.sample(alpha.tsim)
importance.sample(beta.tsim)
summary(tsim2)
curve(dunif(x,-10,10)/dnorm(x,mean=-.5,sd=5),from=-10,10)
dunif(5,-10,10)
## 2.1

metropolis <- function(X,y,iter=500,chains=4){
  
  N <- nrow(X)
  M <- ncol(X)

  computePbeta.y <- function(beta.val){
    # calc p(y|beta) & p(beta)
    py.beta <- dpois(y, exp(X %*% beta.val) )
    #py.beta[which(py.beta==0)] <- NA
    p.beta  <- dcauchy(beta.val,0,2.5)
    # p(beta|y) propto the product:
    return(sum(log(py.beta),na.rm = T)+sum(log(p.beta)))
  }

 compute.new.beta <- function(beta.init){
    
    dens.init <- computePbeta.y(beta.init)
    
    # random walk to next potential beta
    beta.star <- beta.init + rnorm(M,0,.1)
    # compute new density
    dens.star <- computePbeta.y(beta.star)
    
    r <- exp(dens.star - dens.init)
    
    if( runif(1) < min(1,r) ){
      beta.new <- beta.star 
    }else{
      beta.new <- beta.init
    }
    return(list(beta=beta.new,r=r))
  }

  chain <- function(beta){
    results <- list(beta=beta,r=NA)
    new.beta <- list(beta=beta)
    for(i in 1:iter){
      new.beta     <- compute.new.beta(new.beta$beta)
      results$beta <- rbind(results$beta,new.beta$beta)
      results$r    <- rbind(results$r,new.beta$r)
      
    }
    return(results)
  }
  
  init.betas <- lapply(1:chains,function(x) runif(M,-.5,.5))
  
  require(parallel)
  results <- mclapply(init.betas,chain)
  
  return(results)
}

# b
# b
X <- cbind(rnorm(50),rnorm(50),rnorm(50))
y <- rpois(nrow(X),exp(X%*%c(1,.5,-1)))

results <- metropolis(X,y,iter=500)
rs <- lapply(1:4,function(i)results[[i]]$r)
lapply(rs,function(x) mean(x,na.rm=T))
plotdata <- lapply(seq_along(results),function(x){
  df <- as.data.frame(results[[x]]$beta)
  names(df) <- lapply(1:ncol(X),function(i) paste0("beta",i))
  df$chain <- paste0("chain",x)
  return(df)
})

require(reshape2)
plotdata <- do.call("rbind",plotdata)
plotdata.melt <- melt(plotdata,id.vars = "chain")

plotdata.melt$idx <- rep(1:(nrow(plotdata.melt)/(12)),12)

metsims <- ggplot(plotdata.melt,aes(x=idx,y=value,color=chain)) + 
  geom_line() + 
  facet_grid(variable~.) +
  xlab("Iteration") +
  ylab("Parameter Value") +
  ggtitle("Metropolis Simulation")
saveRDS(metsims,"metsims.RDS")


fit3data <- list(N=nrow(X),M=ncol(X),X=X,y=y)
#fit3 <- stan(file='model.2.1.stan',data=fit3data,iter=500)  
fit3 <- stan(fit = fit3,data=fit3data,iter=500)  

results3 <- extract(fit3)
head(results3$betas)
betas <- as.data.frame(results3$betas)
names(betas) <- paste0("beta",1:3)
betas.plot <- melt(betas)
betas.plot$Iteration <- rep(1:(nrow(betas.plot))/3)
stansims <- ggplot(betas.plot,aes(x=value)) + geom_histogram() +   facet_grid(variable~.) +
  xlab("Value") +
  ylab("Frequency") +
  ggtitle("Stan Simulation")
saveRDS(stansims,"stansims.RDS")

plot(results3$betas)
plot(fit3)
print(fit3)
