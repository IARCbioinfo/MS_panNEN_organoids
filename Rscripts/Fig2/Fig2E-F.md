---
title: "Fig2E-F"
author: "N. Alcala"
date: "4/22/2022"
---



# Code to produce Fig. 2E-F from Dayton et al. (Submitted)

## load libraries 

```r
library(tidyverse)
```

```
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
```

```
## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
## ✓ tibble  3.1.6     ✓ dplyr   1.0.8
## ✓ tidyr   1.2.0     ✓ stringr 1.4.0
## ✓ readr   2.1.2     ✓ forcats 0.5.1
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```

```r
library(readxl)
library(trend)
```

## Load and format data 

```r
PassageTimes = read_xlsx("TableS1.xlsx",sheet = 3,skip = 2)
PassageTimes = PassageTimes %>% pivot_longer(-c(Experiment,Grade),names_to = "Passage",values_to = "Date")
PassageTimes = PassageTimes %>% group_by(Experiment) %>% mutate(Time.secs = Date-min(Date,na.rm=T)) %>% ungroup()

## correct LNET24 manually because thawed so dates do not accurately account for passage time
PassageTimes[PassageTimes$Experiment=="LNET 24" & PassageTimes$Passage %in% paste0("P",2:20),]$Time.secs = PassageTimes[PassageTimes$Experiment=="LNET 24" & PassageTimes$Passage %in% paste0("P",2:20),]$Time.secs - (377-152)*60*60*24
```

## Plot Figure 2E

```r
ggplot(PassageTimes, aes(x=Time.secs/60/60/24,y=Experiment,col=Experiment)) + geom_point() + geom_line()+ theme_bw() + 
  geom_vline(xintercept = 365) + xlab("Cumulative days after isolation") + facet_grid(Grade~., scales = "free_y") + guides(col=F)
```

```
## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
## "none")` instead.
```

```
## Don't know how to automatically pick scale for object of type difftime. Defaulting to continuous.
```

```
## Warning: Removed 106 rows containing missing values (geom_point).
```

```
## Warning: Removed 106 row(s) containing missing values (geom_path).
```

![](Fig2E-F_files/figure-html/Figure2E-1.png)<!-- -->

## Test temporal trend in passage times
We use the Mann-Kendall trend test for each experiment

```r
mkl  = lapply(unique(PassageTimes$Experiment) , function(x) mk.test(PassageTimes %>% filter(Experiment==x,!is.na(Time.secs)) %>% pull(Time.secs) %>% as.numeric() %>% diff()))
pvals = sapply(mkl,function(x) x$p.val)
qvals = p.adjust(pvals,method = "BH")
```

We plot the p-value distribution to manually check its uniformity


```r
hist(pvals,nclass = 20) # relatively uniform
```

![](Fig2E-F_files/figure-html/trend_pvals-1.png)<!-- -->

We find the minimal q-value (reported in text)


```r
min(qvals)
```

```
## [1] 0.1809287
```


## Write results 
These results are part of Table S1

```r
mk.tab = sapply(mkl,function(x) unlist(x[c(3,5:6,2)]))
colnames(mk.tab) = unique(PassageTimes$Experiment)
mk.tab = rbind(mk.tab,q.value=qvals)

write.table(mk.tab,file = "TableS1_MKtests.tsv",sep = "\t")
```


# Plot Figure 2F
Plot passage 5 times (note: SINETs did not have P5 times)

```r
P5.times = PassageTimes %>% filter(Passage=="P5", !str_detect(Experiment,"SINET"))
ggplot(P5.times, aes(x=paste0("G",Grade),y=Time.secs/60/60/24,col=Grade)) + geom_point()  + coord_cartesian(ylim=c(0,600)) + theme_bw() + ylab("Time in Days")+
  xlab("Grade")
```

```
## Don't know how to automatically pick scale for object of type difftime. Defaulting to continuous.
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

![](Fig2E-F_files/figure-html/Figure2F-1.png)<!-- -->

## Test difference in passage 5 times 

```r
anova(lm(as.numeric(P5.times$Time.secs)/60/60/24~P5.times$Grade)) 
```

```
## Analysis of Variance Table
## 
## Response: as.numeric(P5.times$Time.secs)/60/60/24
##                Df Sum Sq Mean Sq F value  Pr(>F)  
## P5.times$Grade  2  86567   43284    6.46 0.01821 *
## Residuals       9  60302    6700                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(lm(as.numeric(P5.times$Time.secs)/60/60/24~P5.times$Grade)) 
```

```
## 
## Call:
## lm(formula = as.numeric(P5.times$Time.secs)/60/60/24 ~ P5.times$Grade)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -98.000 -55.938   3.625  27.812 156.500 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)         358.50      40.93   8.759 1.07e-05 ***
## P5.times$Grade2.0  -113.50      57.88  -1.961  0.08152 .  
## P5.times$Grade3.0  -207.75      57.88  -3.589  0.00584 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 81.85 on 9 degrees of freedom
##   (1 observation deleted due to missingness)
## Multiple R-squared:  0.5894,	Adjusted R-squared:  0.4982 
## F-statistic:  6.46 on 2 and 9 DF,  p-value: 0.01821
```

```r
(kruskal.test(as.numeric(P5.times$Time.secs)/60/60/24~P5.times$Grade)) 
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  as.numeric(P5.times$Time.secs)/60/60/24 by P5.times$Grade
## Kruskal-Wallis chi-squared = 7.5385, df = 2, p-value = 0.02307
```

## Session Info 

```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.3.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] trend_1.1.4     readxl_1.3.1    forcats_0.5.1   stringr_1.4.0  
##  [5] dplyr_1.0.8     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0    
##  [9] tibble_3.1.6    ggplot2_3.3.5   tidyverse_1.3.1
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.1.1 xfun_0.29        bslib_0.3.1      haven_2.4.3     
##  [5] colorspace_2.0-2 vctrs_0.3.8      generics_0.1.2   htmltools_0.5.2 
##  [9] yaml_2.2.2       utf8_1.2.2       rlang_1.0.1      jquerylib_0.1.4 
## [13] pillar_1.7.0     withr_2.4.3      glue_1.6.1       DBI_1.1.2       
## [17] dbplyr_2.1.1     modelr_0.1.8     lifecycle_1.0.1  cellranger_1.1.0
## [21] munsell_0.5.0    gtable_0.3.0     rvest_1.0.2      evaluate_0.15   
## [25] labeling_0.4.2   knitr_1.38       tzdb_0.2.0       fastmap_1.1.0   
## [29] fansi_1.0.2      highr_0.9        broom_0.7.12     Rcpp_1.0.8.3    
## [33] backports_1.4.1  scales_1.1.1     jsonlite_1.7.3   farver_2.1.0    
## [37] fs_1.5.2         hms_1.1.1        digest_0.6.29    stringi_1.7.6   
## [41] grid_4.1.2       cli_3.1.1        tools_4.1.2      magrittr_2.0.2  
## [45] sass_0.4.0       crayon_1.4.2     pkgconfig_2.0.3  ellipsis_0.3.2  
## [49] xml2_1.3.3       reprex_2.0.1     lubridate_1.8.0  extraDistr_1.9.1
## [53] assertthat_0.2.1 rmarkdown_2.11   httr_1.4.2       rstudioapi_0.13 
## [57] R6_2.5.1         compiler_4.1.2
```
