# load packages
library(DESeq2)
library(umap)
library(ggnewscale)
library(gghighlight)
library(ggalt)
library(RColorBrewer)

# list high purity samples
high_pur_samples = c("LNET13T","LNET13Tp1","LNET14T","LNET14Tp1","LNET6T","LNET6Tp1","LNET10T","LNET10Tp11","LNET10Tp4","LNET15M","LNET15Mp2","LNET16M","LNET16Mp1",
                     "PANEC1T","PANEC1Tp4","PANEC1Tp14", 
                     "LCNEC11M", "LCNEC11Mp3", "LCNEC23Mp3", "LCNEC3T", "LCNEC3Tp17.2", "LCNEC3Tp24", "LCNEC4T", "LCNEC4Tp24", "LCNEC4Tp7",
                     "SINET12M","SINET12Mp1.1", "SINET12Mp1.3", "SINET21M", "SINET21Mp2", "SINET22M", "SINET22Mp2", 
                     "SINET8M", "SINET8Mp2")

# low pur : "SINET7M", "SINET7Mp2", 

colors_clusters = c("Carcinoid-A1"="#E41A1C", "Carcinoid-A2"= "#377EB8", "Carcinoid-B"= "#4DAF4A", "Supra_carcinoid"="#984EA3",  
                    "LCNEC-like"="#80e5ff80" , #"LCNEC/TypeI"="#FF7F00", "LCNEC/TypeII"="#FFFF33", "LCNEC/NA"="#A65628", "SCLC/LCNEC-like"="#F781BF",
                    "SCLC-like"="#e9c6af80") #"LCNEC/SCLC-like"="#999999", "SCLC/SCLC-like"="gray")

colors_SIgroups = c("liver metastasis"="#E41A1C", "lymph node metastasis"="#377EB8", "mesenteric metastasis"="#4DAF4A","primary"="#984EA3")  
                    
## panel C: NET/NEC umap
### genes
require(readxl)
NatCommunDE.genes = read_xlsx("/data/lungNENomics/work/organoids/metadata/41467_2019_11276_MOESM15_ESM.xlsx",skip = 40)
coreGenesNCDE = NatCommunDE.genes$Gene[which( NatCommunDE.genes$core.genes.A+NatCommunDE.genes$core.genes.B+NatCommunDE.genes$core.genes.LCNEC>0 )]
ID.NatCommunDE.genes50 = unique( unlist( sapply( coreGenesNCDE, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) ) ) )

### expression
samples = c(colnames(gene_1pass.SE.organoids)[colnames(gene_1pass.SE.organoids)%in% high_pur_samples & !str_detect(colnames(gene_1pass.SE.organoids),"SINET")],
            colnames(gene_1pass.SE.GigaSc),colnames(gene_1pass.SE.CarcGeo))

expr_genes_LNEN = DESeqDataSetFromMatrix(assay(cbind(gene_1pass.SE.organoids[,colnames(gene_1pass.SE.organoids)%in% high_pur_samples & !str_detect(colnames(gene_1pass.SE.organoids),"SINET")],
                                                     gene_1pass.SE.GigaSc,gene_1pass.SE.CarcGeo),"counts"), 
                                         design = ~1 , colData = tibble(Sample=samples,Experiment="Test") )
expr_genes_LNEN.vst = varianceStabilizingTransformation(expr_genes_LNEN,blind = T)

write_tsv( bind_cols(gene= rownames(expr_genes_LNEN.vst), as_tibble(assay(expr_genes_LNEN.vst)) ) , "/data/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/expr_matrices/gene_vst_LNET_LCNEC.tsv" )
vvst = apply(assay(expr_genes_LNEN.vst),1,var)
write_tsv( bind_cols(gene= rownames(expr_genes_LNEN.vst[rank(vvst,ties.method = "min")[1:5000],]), as_tibble(assay(expr_genes_LNEN.vst[rank(vvst,ties.method = "min")[1:5000],])) ) , "/data/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/expr_matrices/gene_vst_LNET_LCNEC_5000.tsv" )

write_tsv( as_tibble(rowData(gene_1pass.SE.organoids)) , "/data/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/expr_matrices/rowData.tsv" )


expr_genes_LNEN.vst.mat = t(assay(expr_genes_LNEN.vst))[,ID.NatCommunDE.genes50]

conf = umap.defaults
conf$n_neighbors = nrow(expr_genes_LNEN.vst.mat)
umap_full = as_tibble(umap(expr_genes_LNEN.vst.mat, config = conf )$layout) #,n_neighbors=100
umap_full$Sample = colnames(expr_genes_LNEN.vst)
umap_full$Type = umap_full$`Molecular clusters` = NA
umap_full$Type[sapply(Attributes$Sample_ID,function(x) which(umap_full$Sample==x))] = Attributes$Histopathology_simplified
umap_full$`Molecular clusters`[sapply(Attributes$Sample_ID,function(x) which(umap_full$Sample==x))] = Attributes$Molecular_clusters
umap_full$`Molecular clusters`[1:25] = str_remove(umap_full$Sample[1:25],"[TM][p0-9.]*$")

umap_full = umap_full %>% mutate( `Molecular clusters` = case_when(`Molecular clusters`=="LC1"~"Carcinoid-A1",
                                                       `Molecular clusters`=="LC2"~"Carcinoid-B",
                                                       `Molecular clusters`=="LC3"~"Carcinoid-A2",
                                                       `Molecular clusters`%in%c("LCNEC/TypeI","LCNEC/TypeII","LCNEC/NA","SCLC/LCNEC-like")~"LCNEC-like",
                                                       `Molecular clusters`%in%c("LCNEC/SCLC-like","SCLC/SCLC-like")~"SCLC-like",
                                                       TRUE~`Molecular clusters`) )

colors_org_LP = colors_org[as.character(unique(umap_full[1:25,]$`Molecular clusters`))]

## plot
gg_UMAP_full <- ggplot(umap_full[-(1:25),] %>% mutate(`Molecular clusters`=fct_drop(`Molecular clusters`)),aes(x=V1,y=V2,col=`Molecular clusters`) ) +
  geom_point(shape=16) +#theme(legend.position = "none",legend.key.size = unit(10,"pt")) + 
  scale_color_manual( values=alpha(c(colors_clusters,"black"),0.5))+ theme_minimal()  +
  new_scale_color() +
  geom_path(data=umap_full[1:25,],aes(x=V1,y=V2,col=`Molecular clusters`),inherit.aes = F,size=1.5)+
  geom_point(data=umap_full %>% filter(str_detect(Sample,"^[A-Z]+[0-9]+[MT]$")), aes(col=`Molecular clusters`),size=3) +
  geom_point(data=umap_full %>% filter(str_detect(Sample,"^[A-Z]+[0-9]+[MT]p")), aes(col=`Molecular clusters`),fill="white",shape=21,size=1.5,stroke=1.5) + 
  scale_color_manual( values=colors_org_LP, limits = names(colors_org_LP) ) + 
  xlab("UMAP dimension 1") + ylab("UMAP dimension 2") + 
  xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
  theme(legend.key.size = unit(10,"pt"), legend.title = element_text(size = 10,face = "bold"), 
        legend.text = element_text(size = 8) , legend.spacing = unit(0,'pt'), 
        legend.spacing.y = unit(1,"pt") ) + theme_classic() +
  geom_mark_hull()
gg_UMAP_full


## panel D: SINET umap

### read metadata 
Attributes_SINET_SRA = read_csv("/data/lungNENomics/work/organoids/metadata/SINET_SraRunTable_SI.txt",col_names = F)
Attributes_SINET_EGA = read_xlsx("/data/lungNENomics/work/organoids/metadata/SINET2_7samples.xlsx")

samples = c(colnames(gene_1pass.SE.organoids)[colnames(gene_1pass.SE.organoids)%in% high_pur_samples & str_detect(colnames(gene_1pass.SE.organoids),"SINET")],
            colnames(gene_1pass.SE.SINET),colnames(gene_1pass.SE.SINET2))

### read gene data
require(readxl)
SINETgenes = lapply(1:4, function(i) read_xlsx("/data/lungNENomics/work/organoids/metadata/41588_2018_138_MOESM6_ESM.xlsx",sheet = i)$MR )
SINETgenes = unique(unlist(SINETgenes))
SINETgenes[SINETgenes=="MYCL1"] = "MYCL"
SINETgenes[SINETgenes=="TCEB2"] = "ELOB"
SINETgenes[SINETgenes=="MLL3"] = "KMT2C"
SINETgenes[SINETgenes=="JUB"] = "AJUBA"
SINETgenes[SINETgenes=="GPR128"] = "ADGRG7"
SINETgenes[SINETgenes=="GPR126"] = "ADGRG6"
SINETgenes[SINETgenes=="ZCCHC6"] = "TUT7"
SINETgenes[SINETgenes=="RDBP"] = "NELFE"
SINETgenes[SINETgenes=="SUV420H2"] = "KMT5C"
SINETgenes[SINETgenes=="LASS5"] = "CERS5"
SINETgenes[SINETgenes=="CCDC101"] = "SGF29"
SINETgenes[SINETgenes=="PTRF"] = "CAVIN1"
SINETgenes[SINETgenes=="MS4A8B"] = "MS4A8"
SINETgenes[SINETgenes=="ELTD1"] = "ADGRL4"
SINETgenes[SINETgenes=="WISP2"] = "CCN5"
SINETgenes[SINETgenes=="GPR124"] = "ADGRA2"
SINETgenes[SINETgenes=="DARC"] = "ACKR1"
SINETgenes[SINETgenes=="ERN2"] = "AC008870.1"
SINETgenes[SINETgenes=="GAS7"] = "AC005747.1"
SINETgenes[SINETgenes=="CSDA"] = "YBX3"
SINETgenes = SINETgenes[which(SINETgenes %in% rowData(gene_1pass.SE.organoids)$gene_name )]

ID.SINETtypes.genes = sapply( SINETgenes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )

### expression
expr_genes_SINET = DESeqDataSetFromMatrix(assay(cbind(gene_1pass.SE.organoids[,colnames(gene_1pass.SE.organoids)%in% high_pur_samples & str_detect(colnames(gene_1pass.SE.organoids),"SINET")],
                                                     gene_1pass.SE.SINET,gene_1pass.SE.SINET2),"counts"), 
                                         design = ~1 , colData = tibble(Sample=samples,Experiment="Test") )
expr_genes_SINET.vst = varianceStabilizingTransformation(expr_genes_SINET,blind = T)

write_tsv( bind_cols(gene= rownames(expr_genes_SINET.vst), as_tibble(assay(expr_genes_SINET.vst)) ) , "/data/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/expr_matrices/gene_vst_SINET.tsv" )
vvst = apply(assay(expr_genes_SINET.vst),1,var)
write_tsv( bind_cols(gene= rownames(expr_genes_SINET.vst[rank(vvst,ties.method = "min")[1:5000],]), as_tibble(assay(expr_genes_SINET.vst[rank(vvst,ties.method = "min")[1:5000],])) ) , "/data/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/expr_matrices/gene_vst_SINET_5000.tsv" )

expr_genes_SINET.vst.mat = t(assay(expr_genes_SINET.vst))
expr_genes_SINET.vst.mat = expr_genes_SINET.vst.mat[,ID.SINETtypes.genes]

conf = umap.defaults
conf$n_neighbors = nrow(expr_genes_SINET.vst.mat)
umap_fullSI = as_tibble(umap(expr_genes_SINET.vst.mat, config = conf )$layout) #,n_neighbors= nrow(expr_genes_SINET.vst.mat))#100)
umap_fullSI$Sample = colnames(expr_genes_SINET.vst)
umap_fullSI$Type = NA
umap_fullSI$Type[sapply(Attributes_SINET_SRA$X1,function(x) which(umap_fullSI$Sample==x))] = Attributes_SINET_SRA$X28
umap_fullSI$Type[1:9] = str_remove(umap_fullSI$Sample[1:9],"[TM][p0-9.]*$")
umap_fullSI$Type[91:97] = Attributes_SINET_EGA$`Tumour Metastasis Site`
umap_fullSI$Type[umap_fullSI$Type=="LN"] = "lymph node metastasis"
umap_fullSI$Type[umap_fullSI$Type=="Hepatic"] = "liver metastasis"

colors_org_SI = colors_org[as.character(unique(umap_fullSI[1:11,]$Type))]

## plot data
gg_umap_fullSI <- ggplot(umap_fullSI[-(1:9),] %>% mutate(Type=fct_drop(Type)),aes(x=V1,y=V2,col=Type) ) +
  geom_point(shape=16) +#theme(legend.position = "none",legend.key.size = unit(10,"pt")) + 
  scale_color_manual( values=alpha(c(colors_SIgroups,"black"),0.5))+ theme_minimal()  +
  new_scale_color() +
  geom_path(data=umap_fullSI[1:9,],aes(x=V1,y=V2,col=Type),inherit.aes = F,size=1.5)+
  geom_point(data=umap_fullSI %>% filter(str_detect(Sample,"^[A-Z]+[0-9]+[MT]$")), aes(col=Type),size=3) +
  geom_point(data=umap_fullSI %>% filter(str_detect(Sample,"^[A-Z]+[0-9]+[MT]p")), aes(col=Type),fill="white",shape=21,size=1.5,stroke=1.5) + 
  scale_color_manual( values=colors_org_SI, limits = names(colors_org_SI) ) + 
  xlab("UMAP dimension 1") + ylab("UMAP dimension 2") + 
  xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
  theme(legend.key.size = unit(10,"pt"), legend.title = element_text(size = 10,face = "bold"), 
        legend.text = element_text(size = 8) , legend.spacing = unit(0,'pt'), 
        legend.spacing.y = unit(1,"pt") ) + theme_classic() +
  geom_mark_hull()

## save
ggsave("/data/lungNENomics/work/organoids/figures/Fig3CD_raw.svg", gg_umap_fullSI + gg_UMAP_full,height = 2.2*1.2,width = 4*2*1.05)

write_tsv(gg_UMAP_full ,path = "Fig3C_draft_data_SINET.tsv")
write_tsv(gg_umap_fullSI ,path = "Fig3D_draft_data_SINET.tsv")
