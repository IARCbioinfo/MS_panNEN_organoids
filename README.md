# PanNEN manuscript scripts
This repository contains the scripts used to produce Figures in panNEN organoids manuscript Dayton et al. (Submitted).

## Rscripts
Major package dependencies are mentioned below; see a list of all dependencies at the beginning of each script.

### Figure 2. NET and LCNEC PDTOs retain histologic features and relative growth-rate of parental tumor subtypes
#### [Figure2D.md](Rscripts/Fig2/Fig2D.md)
- requires [ggridges R package](https://cran.r-project.org/web/packages/ggridges/)

Produces density plots of expression levels of *MKI67* in organoids and reference tumors (data in Table S1).

#### [Figure2E-F.md](Rscripts/Fig2/Fig2E-F.md)
- requires [trend R package](https://cran.r-project.org/web/packages/trend/index.html)

Produces plots of passage times, computes temporal trend tests and writes results (Table S1).

### Figure 3. High-purity NEN PDTOs recapitulate the gene expression of original tumors
#### [Figure3B.md](Rscripts/Fig3/Fig3B_S3BCE.md)
- requires [ggbeeswarm R package](https://cran.r-project.org/web/packages/ggbeeswarm/index.html)

Produces violin plots of gene expression for various markers from Table S2, producing Figure 3B, S3B, C, and E.

#### Figure3CD
Produces UMAP representations of lung and pancreatic NENs and small intestine NETs (using umap).

### Figure 4. NEN PDTOs retain genomic features of parental tumors
#### Figure4BC
Produces oncoplots of somatic alterations from whole-genome and RNA-seq data (package maftools). 

#### Figure4D
Produces circos plots of copy number variants and structural variants (package circos). 

### Figure 5. NEN PDTOs recapitulate the intra-tumor heterogeneity of the parental tumor
#### Figure5A
...

## Citation
Dayton*, Alcala, ... , Foll, Fernandez-Cuesta*, Clevers*. 2022 Submitted.
