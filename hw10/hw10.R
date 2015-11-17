setwd("~/Desktop/school/bayesian-data-analysis/hw10")

mortRate <- read.table('cdc/white_nonhisp_death_rates_from_1999_to_2013_by_sex.txt',header=T)

library(dplyr)
meanRatesAgeYear <- mortRate %>% 
  mutate(ageBin=round(Age/10)*10) %>%
  group_by(ageBin,Age,Year) %>% 
  summarize(avgRate=1000*sum(Deaths)/sum(Population),agePop=sum(Population)) %>%
  group_by(ageBin,Year) %>%
  mutate(weight=agePop/sum(agePop))

adjustedRates <- meanRatesAgeYear %>%
  group_by(ageBin,Year) %>%
  summarize(adjRate=sum(weight*avgRate),
            avgRate=mean(avgRate)) %>%
  group_by(ageBin) %>%
  mutate(adjStartVal=adjRate[which.min(Year)],
         avgStartVal=avgRate[which.min(Year)],
         adjIdx=adjRate/adjStartVal,
         avgIdx=avgRate/avgStartVal)

library(reshape2)
ggAdjRates <- melt(adjustedRates[,c('ageBin','Year','adjIdx','avgIdx')],id.vars = c("ageBin","Year"))

library(ggplot2)
ggplot(ggAdjRates,aes(x=Year,y=value,color=variable)) +
  geom_line() +
  geom_hline(y=1) +
  facet_wrap(~ageBin,nrow=3) +
  ggtitle("Change in Mortality Rates from 1999-2011\nBy Age Group")