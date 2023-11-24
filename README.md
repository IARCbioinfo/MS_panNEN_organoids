# PanNEN manuscript scripts
This repository contains the data and scripts used to produce the genomic Figures in the panNEN organoids manuscript Dayton et al. (In press; https://www.biorxiv.org/content/10.1101/2022.10.31.514549v1) and its associated data note Alcala et al. biorxiv (2023; https://www.biorxiv.org/content/10.1101/2023.08.31.555732v1).

## Data pre-processing
We provide here the command lines used to process the RNA-sequencing and whole-genome sequencing (WGS) data. All processing was performed using the Nextflow pipelines freely available on this IARCbioinfo github account.

### RNA-seq

Step 1: mapping

Step 2: post-processing

Step 3: expression quantification

### WGS

Step 1: mapping

Step 2: post-processing

Step 3: variant calling

## Data
The data folder contains tab-separated files with processed data:
- *gene_expression_PDTOs_parents.tsv* contains the final gene expression matrix, with read counts for the 59607 features and 47 samples sequenced in the Dayton et al. study
- *gene_expression_PDTOs_parents.tsv* contains the annotation of the 59607 features (gene name, ensembl ID, protein ID, gene type, etc)

## Rscripts
The Rscripts folder contains markdown files detailing R commands used to produce the figures. Major package dependencies are mentioned below; see a list of all dependencies at the beginning of each script.

### Figure 2. NET and LCNEC PDTOs retain histologic features and relative growth-rate of parental tumor subtypes
#### [Figure2D.md](Rscripts/Fig2/Fig2D.md)
- requires [ggridges R package](https://cran.r-project.org/web/packages/ggridges/)

Produces density plots of expression levels of *MKI67* in organoids and reference tumors (data in Table S1).

![](Rscripts/Fig2/Fig2D_files/figure-html/Figure2D-1.png)

#### [Figure2E-F.md](Rscripts/Fig2/Fig2E-F.md)
- requires [trend R package](https://cran.r-project.org/web/packages/trend/index.html)

Produces plots of passage times, computes temporal trend tests and writes results (Table S1).
![](Rscripts/Fig2/Fig2E-F_files/figure-html/Figure2E-1.png)
![](Rscripts/Fig2/Fig2E-F_files/figure-html/Figure2F-1.png)

### Figure 3. High-purity NEN PDTOs recapitulate the gene expression of original tumors
#### [Figure3B_S3BCE.md](Rscripts/Fig3/Fig3B_S3BCE.md)
- requires [ggbeeswarm R package](https://cran.r-project.org/web/packages/ggbeeswarm/index.html)

Produces violin plots of gene expression for various markers from Table S2, producing Figure 3B, S3B, C, and E.

![](Rscripts/Fig3/Fig3B_S3BCE_files/figure-html/Figure3B-1.png)

#### [Figure3CD_S3FGHI.md](Rscripts/Fig3/Fig3CD_S3FGHI.md)
- requires [umap R package](https://cran.r-project.org/web/packages/umap/vignettes/umap.html)
- requires [mixOmics R package](http://mixomics.org/)

Produces the unsupervised (UMAP representations) and supervised (PLS) analyses of lung and pancreatic NENs and small intestine NETs presented in Figures 3CD and S3F-I.

![](Rscripts/Fig3/Fig3CD_S3FGH_files/figure-html/Figure3C-1.png)
![](Rscripts/Fig3/Fig3CD_S3FGH_files/figure-html/Figure3D-1.png)

### Figure 4. NEN PDTOs retain genomic features of parental tumors
#### [Figure4BC_S4BC.md](Rscripts/Fig4/Fig4BC_S4BC.md)
- requires [maftools R package](https://https://bioconductor.org/packages/release/bioc/html/maftools.html)

Produces oncoplots of somatic alterations from whole-genome and RNA-seq data and tumor mutational burden plots. 
![](Rscripts/Fig4/Fig4BC_S4BC_files/figure-html/Fig4B-1.png)
![](Rscripts/Fig4/Fig4BC_S4BC_files/figure-html/Fig4Bbot-1.png)
![](Rscripts/Fig4/Fig4BC_S4BC_files/figure-html/Fig4C-1.png)

#### [Figure4D_S4D.md](Rscripts/Fig4/Fig4D_S4D.md)
- requires [circlize R package](https://jokergoo.github.io/circlize_book/book/)

Produces circos plots of copy number variants and structural variants in high-purity (Fig. 4D) and mixed (Fig. S4D) samples. 

![](Rscripts/Fig4/Fig4D_S4D_files/figure-html/Fig4D-4.png)

### Figure 5. NEN PDTOs recapitulate the intra-tumor heterogeneity of the parental tumor
#### [Fig5_S5.md](Rscripts/Fig5/Fig5_S5.md)
- requires [DPclust R package](https://github.com/Wedge-lab/dpclust)

Produces Venn-Euler diagrams of shared variants (Fig. 5A and S5B) and joint plots of cancer cell fractions (CCFs).

![](Rscripts/Fig5/Fig5_S5_files/figure-html/fig5AS5B-4.png)
![](Rscripts/Fig5/Fig5_S5_files/figure-html/S5A-7.png)
![](Rscripts/Fig5/Fig5_S5_files/figure-html/evolplots-7.png)
![](Rscripts/Fig5/Fig5_S5_files/figure-html/Fig5D-1.png)
![](Rscripts/Fig5/Fig5_S5_files/figure-html/Fig5E-1.png)

### Figure 7. LNETs express EGFR
#### [FigureS7B.md](Rscripts/Fig7/FigS7B.md)
- requires [ggbeeswarm R package](https://cran.r-project.org/web/packages/ggbeeswarm/index.html)

Produces violin plots of EGFR gene expression from Table S2, producing Figure S7B.

![](Rscripts/Fig7/FigS7B_files/figure-html/FigureS7B-1.png)

### Data Note Figure 6. RNA-seq variant classification using a random forest algorithm
#### [DataNote_Fig6.md](Rscripts/DataNoteFig6/DataNote_Fig6.md)
- requires [caret R package](https://cran.r-project.org/web/packages/caret/index.html)

Produces ROC curve, confusion matrices, and feature importance metrics in Figure 6.

![](Rscripts/DataNoteFig6/DataNote_Fig6_files/figure-html/ROC-1.png)


## Citations
Dayton*, Alcala, ... , Foll, Fernandez-Cuesta*, Clevers*. Druggable Growth Dependencies and Tumor Evolution Analysis in Patient-Derived Organoids of Neuroendocrine Cancer. biorxiv 2022. doi: https://doi.org/10.1101/2022.10.31.514549

Alcala*, Voegele, ..., Dayton, Foll*. Multi-omic dataset of patient-derived tumor organoids of neuroendocrine neoplasms. biorxiv 2023. doi: https://doi.org/10.1101/2023.08.31.555732

