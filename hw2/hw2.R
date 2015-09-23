setwd('~/Desktop/school/bayesian-data-analysis/hw2')

data <- read.table("http://www.stat.columbia.edu/~gelman/bda.course/electric.txt",skip=2)
names(data) <- c("city","grade","ePre","ePost","cPre","cPost","replace.suppliment")



# compare discriptive stats
library(dplyr)
summarized <- data %>% group_by(grade) %>% summarise(e.premean=mean(e.pretest),
                                       e.postmean=mean(e.posttest),
                                       c.premean=mean(c.pretest),
                                       c.postmean=mean(c.posttest))
summarized <- within(summarized,{e.diff <-e.postmean-e.premean
                 c.diff <-c.postmean-c.premean
})


# prep for stan
for(i in 1:ncol(data)){
  assign(colnames(data)[i],data[,i])
}
data$grade <- as.factor(data$grade)
grade <- model.matrix(~grade,data=data)
N <- nrow(data)
K <- length(unique(data$grade))
city <- as.vector(model.matrix(~city,data)[,2])
ones <- rep(1,N)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fit <- stan(file="model1.stan",
            data=c("cPost","cPre","ePost","ePre","grade","city","ones","N","K"),chains=4,iter=100)
saveRDS(fit,"fit1.rds")
plot(fit)
print(fit)
params <- extract(fit)

mean(params$tdiff)



#########
# part 2
#########


