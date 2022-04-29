# PanNEN manuscript scripts
This repository contains the scripts used to produce Figures in panNEN organoids manuscript Dayton et al. (Submitted).

## Rscripts
Major package dependencies are mentioned below; see a list of all dependencies at the beginning of each script.

### Figure 2
#### Figure2D.R
- requires [ggbeeswarm R package](https://cran.r-project.org/web/packages/ggbeeswarm/index.html)

Produces density plots of expression levels.

#### [Figure2E-F.md](Rscripts/Fig2/Fig2E-F.md)
- requires [trend R package](https://cran.r-project.org/web/packages/trend/index.html)

Produces plots of passage times, computes temporal trend tests and writes results (Table S1).

### Figure 3
#### Figure3B.R
Produces violin plots of gene expression (using package [ggbeeswarm](https://github.com/eclarke/ggbeeswarm)).

#### Figure3CD.R
Produces UMAP representations of lung and pancreatic NENs and small intestine NETs (using umap).

### Figure 4
#### Figure4BC.R
Produces oncoplots of somatic alterations from whole-genome and RNA-seq data (package maftools). 

#### Figure4D.R
Produces circos plots of copy number variants and structural variants (package circos). 

## Citation
Dayton*, Alcala, ... , Foll, Fernandez-Cuesta*, Clevers*. 2022 Submitted.
