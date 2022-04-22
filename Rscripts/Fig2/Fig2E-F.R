####################################
## Code to produce Fig. 2E-F      ##
####################################

# load libraries 
library(tidyverse)
library(readxl)
library(trend)

# load data 
PassageTimes = read_xlsx("PDTO_timeInCulture.xlsx")
## format data
PassageTimes = PassageTimes %>% pivot_longer(-c(Experiment,Grade),names_to = "Passage",values_to = "Date")
PassageTimes = PassageTimes %>% group_by(Experiment) %>% mutate(Time.secs = Date-min(Date,na.rm=T)) %>% ungroup()

## correct LNET24 manually because thawed so dates do not accurately account for passage time
PassageTimes[PassageTimes$Experiment=="LNET 24" & as.numeric(PassageTimes$Passage)>=2,]$Time.secs = PassageTimes[PassageTimes$Experiment=="LNET 24" & as.numeric(PassageTimes$Passage)>=2,]$Time.secs - (377-152)*60*60*24

# Code for Figure 2E
ggplot(PassageTimes, aes(x=Time.secs/60/60/24,y=Experiment,col=Experiment)) + geom_point() + geom_line()+ theme_bw() + 
  geom_vline(xintercept = 365) + xlab("Cumulative days after isolation") + facet_grid(Grade~., scales = "free_y") + guides(col=F)


## test acceleration
mkl  = lapply(unique(PassageTimes$Experiment) , function(x) mk.test(PassageTimes %>% filter(Experiment==x,!is.na(Time.secs)) %>% pull(Time.secs) %>% as.numeric() %>% diff()))
pvals = sapply(mkl,function(x) x$p.val)
qvals = p.adjust(pvals,method = "BH")

## plot p-value distribution
hist(pvals,nclass = 20) # relatively uniform

## find minimal q-value
min(qvals)

## write results 
mk.tab = sapply(mkl,function(x) unlist(x[c(3,5:6,2)]))
colnames(mk.tab) = unique(PassageTimes$Experiment)
mk.tab = rbind(mk.tab,q.value=qvals)

write.table(mk.tab,file = "Fig2_MKtests.tsv",sep = "\t") #included in TableS1

# Code for Figure 2F
P5.times = PassageTimes %>% filter(Passage=="5")
ggplot(P5.times, aes(x=Grade,y=Time.secs/60/60/24,col=Grade)) + geom_point()  + coord_cartesian(ylim=c(0,600)) + theme_bw()

## test significance of difference, reported in the text
anova(lm(as.numeric(P5.times$Time.secs)/60/60/24~P5.times$Grade)) 