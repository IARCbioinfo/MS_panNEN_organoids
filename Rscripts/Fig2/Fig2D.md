---
title: "Fig2D"
author: "N. Alcala"
date: "4/22/2022"
output:
  html_document: 
    keep_md: yes
  pdf_document: default
---


# Code to produce Fig. 2D from Dayton et al. (Submitted)

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
library(ggridges)
```



## Define colors
Colors for each experiment (organoid and parental tumor families), type, and grades:

```r
colors_org = c(LNET2="#aade87ff",LNET6="#5fd38dff",LNET13="#16502dff",LNET14="#6f917cff",
               LNET5="#e6a73cff",LNET10="#ff9955ff",LNET15="#ffd42aff", LNET16 = "#ff6600ff", LNET18= "#d0742fff", 
               LNET19="#2aff80ff", 
               LNET20 = "#f6e62bff", 
               LCNEC3="#ff8080ff",LCNEC4="#d35f5fff", LCNEC23 = "#ff5555ff", 
               LCNEC11="#ff5599ff",PANEC1="#8d5fd3ff",
               SINET7="#2ad4ffff",SINET8="#80b3ffff",SINET9="#5f8dd3ff",SINET12="#5fbcd3ff", SINET21="#0066ffff", SINET22="#2c5aa0ff")

colors_types = c(SI= "#5bc2f0ff", Lung = "#9b9972ff", Pancreas = "#8d5fd3ff")

colors_grades = c(G1= "#58b873f9", G2 = "#ff9955ff", G3 = "#f0677dff", "G1/G2" = "#58c1f0ff")
```

## Load data 

```r
expr_genes_KI67.tib = read_xlsx("/data/lungNENomics/work/organoids/SI/TableS1.xlsx",sheet = 2,skip=2,col_types = c("text","numeric","text","text","text"))
```

## Plot Figure 2D
Figure presenting the expression of *MKI67* in organoids and reference tumors:

```r
ggplot( expr_genes_KI67.tib %>% filter(Experiment=="Reference",Type!="SCLC",Grade!="NA") , aes(y=Grade,x=Expression.TPM,fill=Grade) ) + geom_density_ridges2(scale = 1.5,col=NA) + 
  geom_point(data = expr_genes_KI67.tib %>% filter(Experiment!="Reference",str_detect(Sample,"p[0-9.]*$")) , size=4, pch=16, col="white",fill="black") +
  geom_point(data = expr_genes_KI67.tib %>% filter(Experiment!="Reference",str_detect(Sample,"p[0-9.]*$")) , size=2.7, pch=16, col="black",fill="black") + 
  theme_classic()  + labs(y="Histopathological type",x=expression(italic(MKI67)~" Expression (TPM)") ) + 
  geom_vline(xintercept = 1,linetype="dashed") + 
  scale_x_log10(breaks=c(0.01,0.1,1,10,10**2),limits=c(0.01,200),labels=c("≤0.01",0.1,1,10,100)) + scale_fill_manual(values=alpha(colors_grades,0.5))
```

![](Fig2D_files/figure-html/Figure2D-1.png)<!-- -->

## compute a few stats
Compute mean expression levels per grade and tumor type, in organoids and reference samples

```r
expr_genes_KI67.tib %>% group_by(Experiment=="Reference",Type,Grade) %>% summarize(mean(Expression.TPM))
```

```
## # A tibble: 10 × 4
## # Groups:   Experiment == "Reference", Type [7]
##    `Experiment == "Reference"` Type  Grade `mean(Expression.TPM)`
##    <lgl>                       <chr> <chr>                  <dbl>
##  1 FALSE                       LCNEC G3                     39.0 
##  2 FALSE                       LNET  G1                      1.15
##  3 FALSE                       LNET  G2                      5.51
##  4 FALSE                       SINET G1/G2                   1.02
##  5 TRUE                        LCNEC G3                     37.8 
##  6 TRUE                        LNET  G1                      1.10
##  7 TRUE                        LNET  G2                   4268.  
##  8 TRUE                        LNET  NA                      2.05
##  9 TRUE                        SCLC  G3                     71.4 
## 10 TRUE                        SINET G1/G2                   3.54
```

```r
expr_genes_KI67.tib %>% filter(str_detect(Sample,"p[0-9]")) %>% group_by(Experiment=="Reference",Type,Grade) %>% summarize(mean(Expression.TPM))
```

```
## # A tibble: 4 × 4
## # Groups:   Experiment == "Reference", Type [3]
##   `Experiment == "Reference"` Type  Grade `mean(Expression.TPM)`
##   <lgl>                       <chr> <chr>                  <dbl>
## 1 FALSE                       LCNEC G3                    38.8  
## 2 FALSE                       LNET  G1                     1.45 
## 3 FALSE                       LNET  G2                     4.28 
## 4 FALSE                       SINET G1/G2                  0.691
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
##  [1] ggridges_0.5.3  readxl_1.3.1    forcats_0.5.1   stringr_1.4.0  
##  [5] dplyr_1.0.8     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0    
##  [9] tibble_3.1.6    ggplot2_3.3.5   tidyverse_1.3.1
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.1.1 xfun_0.29        bslib_0.3.1      haven_2.4.3     
##  [5] colorspace_2.0-2 vctrs_0.3.8      generics_0.1.2   htmltools_0.5.2 
##  [9] yaml_2.2.2       utf8_1.2.2       rlang_1.0.1      jquerylib_0.1.4 
## [13] pillar_1.7.0     withr_2.4.3      glue_1.6.1       DBI_1.1.2       
## [17] dbplyr_2.1.1     modelr_0.1.8     plyr_1.8.6       lifecycle_1.0.1 
## [21] cellranger_1.1.0 munsell_0.5.0    gtable_0.3.0     rvest_1.0.2     
## [25] evaluate_0.15    knitr_1.38       tzdb_0.2.0       fastmap_1.1.0   
## [29] fansi_1.0.2      highr_0.9        broom_0.7.12     Rcpp_1.0.8.3    
## [33] backports_1.4.1  scales_1.1.1     jsonlite_1.7.3   farver_2.1.0    
## [37] fs_1.5.2         hms_1.1.1        digest_0.6.29    stringi_1.7.6   
## [41] grid_4.1.2       cli_3.1.1        tools_4.1.2      magrittr_2.0.2  
## [45] sass_0.4.0       crayon_1.4.2     pkgconfig_2.0.3  ellipsis_0.3.2  
## [49] xml2_1.3.3       reprex_2.0.1     lubridate_1.8.0  assertthat_0.2.1
## [53] rmarkdown_2.11   httr_1.4.2       rstudioapi_0.13  R6_2.5.1        
## [57] compiler_4.1.2
```

