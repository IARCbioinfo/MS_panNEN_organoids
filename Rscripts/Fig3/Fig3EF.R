library(DESeq2)
library(umap)
library(ggnewscale)
library(gghighlight)
library(ggalt)
library(RColorBrewer)

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
                    
## panel B: NET/NEC umap
### genes
require(readxl)
NatCommunDE.genes = read_xlsx("/data/gcs/lungNENomics/work/organoids/metadata/41467_2019_11276_MOESM15_ESM.xlsx",skip = 40)
coreGenesNCDE = NatCommunDE.genes$Gene[which( NatCommunDE.genes$core.genes.A+NatCommunDE.genes$core.genes.B+NatCommunDE.genes$core.genes.LCNEC>0 )]
ID.NatCommunDE.genes50 = unique( unlist( sapply( coreGenesNCDE, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) ) ) )

### expression
samples = c(colnames(gene_1pass.SE.organoids)[colnames(gene_1pass.SE.organoids)%in% high_pur_samples & !str_detect(colnames(gene_1pass.SE.organoids),"SINET")],
            colnames(gene_1pass.SE.GigaSc),colnames(gene_1pass.SE.CarcGeo))

expr_genes_LNEN = DESeqDataSetFromMatrix(assay(cbind(gene_1pass.SE.organoids[,colnames(gene_1pass.SE.organoids)%in% high_pur_samples & !str_detect(colnames(gene_1pass.SE.organoids),"SINET")],
                                                     gene_1pass.SE.GigaSc,gene_1pass.SE.CarcGeo),"counts"), 
                                         design = ~1 , colData = tibble(Sample=samples,Experiment="Test") )
expr_genes_LNEN.vst = varianceStabilizingTransformation(expr_genes_LNEN,blind = T)

write_tsv( bind_cols(gene= rownames(expr_genes_LNEN.vst), as_tibble(assay(expr_genes_LNEN.vst)) ) , "/data/gcs/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/expr_matrices/gene_vst_LNET_LCNEC.tsv" )
vvst = apply(assay(expr_genes_LNEN.vst),1,var)
write_tsv( bind_cols(gene= rownames(expr_genes_LNEN.vst[rank(vvst,ties.method = "min")[1:5000],]), as_tibble(assay(expr_genes_LNEN.vst[rank(vvst,ties.method = "min")[1:5000],])) ) , "/data/gcs/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/expr_matrices/gene_vst_LNET_LCNEC_5000.tsv" )

write_tsv( as_tibble(rowData(gene_1pass.SE.organoids)) , "/data/gcs/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/expr_matrices/rowData.tsv" )


expr_genes_LNEN.vst.mat = t(assay(expr_genes_LNEN.vst))[,ID.NatCommunDE.genes50]
#vv = apply(expr_genes_LNEN.vst.mat,2,var)
#expr_genes_LNEN.vst.mat = expr_genes_LNEN.vst.mat[,order(vv,decreasing = T)[1:5000]]

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
  #geom_encircle(data = subset(umap_full,`Molecular clusters`=="Carcinoid-A1/LC1") , colour=distinctive_cols[1], fill=distinctive_cols[1], s_shape=0.5, expand=0.02)+ 
  #geom_encircle(data = subset(umap_full,`Molecular clusters`=="Carcinoid-A2/LC3") , colour=distinctive_cols[2], fill=distinctive_cols[2], s_shape=0.5, expand=0.02)+ 
  #geom_encircle(data = subset(umap_full,`Molecular clusters`=="Carcinoid-B/LC2") , colour=distinctive_cols[3], fill=distinctive_cols[3], s_shape=0.5, expand=0.02)+ 
  #geom_encircle(data = subset(umap_full,`Molecular clusters`=="LCNEC-like (LCNEC & SCLC)") , colour=distinctive_cols[5], fill=distinctive_cols[5], s_shape=0.5, expand=0.02)+ 
  #geom_encircle(data = subset(umap_full,`Molecular clusters`=="SCLC-like (LCNEC & SCLC)") , colour=distinctive_cols[6], fill=distinctive_cols[6], s_shape=0.5, expand=0.02)
gg_UMAP_full


## panel F: SINET umap

### read metadata 
Attributes_SINET_SRA = read_csv("/data/gcs/lungNENomics/work/organoids/metadata/SINET_SraRunTable_SI.txt",col_names = F)
Attributes_SINET_EGA = read_xlsx("/data/gcs/lungNENomics/work/organoids/metadata/SINET2_7samples.xlsx")

samples = c(colnames(gene_1pass.SE.organoids)[colnames(gene_1pass.SE.organoids)%in% high_pur_samples & str_detect(colnames(gene_1pass.SE.organoids),"SINET")],
            colnames(gene_1pass.SE.SINET),colnames(gene_1pass.SE.SINET2))

### read gene data
require(readxl)
SINETgenes = lapply(1:4, function(i) read_xlsx("/data/gcs/lungNENomics/work/organoids/metadata/41588_2018_138_MOESM6_ESM.xlsx",sheet = i)$MR )
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


write_tsv( bind_cols(gene= rownames(expr_genes_SINET.vst), as_tibble(assay(expr_genes_SINET.vst)) ) , "/data/gcs/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/expr_matrices/gene_vst_SINET.tsv" )
vvst = apply(assay(expr_genes_SINET.vst),1,var)
write_tsv( bind_cols(gene= rownames(expr_genes_SINET.vst[rank(vvst,ties.method = "min")[1:5000],]), as_tibble(assay(expr_genes_SINET.vst[rank(vvst,ties.method = "min")[1:5000],])) ) , "/data/gcs/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/expr_matrices/gene_vst_SINET_5000.tsv" )

expr_genes_SINET.vst.mat = t(assay(expr_genes_SINET.vst))
#vv = apply(expr_genes_SINET.vst.mat,2,var)
expr_genes_SINET.vst.mat = expr_genes_SINET.vst.mat[,ID.SINETtypes.genes]#[,order(vv,decreasing = T)[1:5000]]

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
#geom_encircle(data = subset(umap_fullSI,`Molecular clusters`=="Carcinoid-A1/LC1") , colour=distinctive_cols[1], fill=distinctive_cols[1], s_shape=0.5, expand=0.02)+ 
#geom_encircle(data = subset(umap_fullSI,`Molecular clusters`=="Carcinoid-A2/LC3") , colour=distinctive_cols[2], fill=distinctive_cols[2], s_shape=0.5, expand=0.02)+ 
#geom_encircle(data = subset(umap_fullSI,`Molecular clusters`=="Carcinoid-B/LC2") , colour=distinctive_cols[3], fill=distinctive_cols[3], s_shape=0.5, expand=0.02)+ 
#geom_encircle(data = subset(umap_fullSI,`Molecular clusters`=="LCNEC-like (LCNEC & SCLC)") , colour=distinctive_cols[5], fill=distinctive_cols[5], s_shape=0.5, expand=0.02)+ 
#geom_encircle(data = subset(umap_fullSI,`Molecular clusters`=="SCLC-like (LCNEC & SCLC)") , colour=distinctive_cols[6], fill=distinctive_cols[6], s_shape=0.5, expand=0.02)


ggsave("/data/gcs/lungNENomics/work/organoids/figures/Fig3BC_raw.svg", gg_umap_fullSI + gg_UMAP_full,height = 2.2*1.2,width = 4*2*1.05)


## old
require(umap)
meta_SINET = read_csv("SINET_data/SraRunTable_SI.txt",col_names = F) #,col_names = c("ID","Technique","read_length","ID2","ID3","ID4","Repo","access","format",""))

Attributes.SINETtypes = Attributes.VST[colnames(expr_genes_SINETtypes),]
Attributes.SINETtypes[meta_SINET$X1,]$Molecular_clusters = meta_SINET$X28
Attributes.SINETtypes$Molecular_clusters[Attributes.SINETtypes$Molecular_clusters=="SINET"] = "primary"
Attributes.SINETtypes$Molecular_clusters = factor(gsub("Experiment ","",Attributes.SINETtypes$Molecular_clusters), levels=c("SINET7","SINET8","SINET12","SINET21","SINET22","primary","liver metastasis","lymph node metastasis","mesenteric metastasis"))

umap_SINETgenes = as_tibble(umap(t(expr_genes_SINETtypes),n_neighbors=ncol(expr_genes_SINETtypes) )$layout) 
colnames(umap_SINETgenes) = c("UMAP_dimension1","UMAP_dimension2")
umap_SINETgenes$`Molecular clusters` = Attributes.SINETtypes$Molecular_clusters
umap_SINETgenes = umap_SINETgenes %>% mutate(`Molecular clusters`= fct_drop(fct_relabel(`Molecular clusters`,~gsub("Experiment ","",levels(`Molecular clusters`))) ) )

colors_org_SI = colors_org[as.character(unique(umap_SINETgenes[1:11,]$`Molecular clusters`))]

gg_UMAP_SINETgenes <- ggplot(umap_SINETgenes[-(1:11),] %>% mutate(`Molecular clusters`=fct_drop(`Molecular clusters`)),
                             aes(x=UMAP_dimension1,y=UMAP_dimension2,col=`Molecular clusters`) ) +
  geom_point() +theme(legend.position = "none",legend.key.size = unit(10,"pt")) + 
  scale_color_manual( values=alpha(c(distinctive_cols[9:12],"black"),0.8))+ theme_classic()  +
  geom_encircle(data = subset(umap_SINETgenes,`Molecular clusters`=="primary") , colour=alpha(distinctive_cols[9],0.5), fill=alpha(distinctive_cols[9],0.5), s_shape=0.5, expand=0.02)+ 
  geom_encircle(data = subset(umap_SINETgenes,`Molecular clusters`=="liver metastasis") , colour=alpha(distinctive_cols[10],0.5), fill=alpha(distinctive_cols[10],0.5), s_shape=0.5, expand=0.02)+
  geom_encircle(data = subset(umap_SINETgenes,`Molecular clusters`=="lymph node metastasis") , colour=alpha(distinctive_cols[11],0.5), fill=alpha(distinctive_cols[11],0.5), s_shape=0.5, expand=0.02)+ 
  new_scale_color() +
  geom_path(data=umap_SINETgenes[c(1,2),],aes(x=UMAP_dimension1,y=UMAP_dimension2,col=`Molecular clusters`),inherit.aes = F,size=1.5)+
  geom_path(data=umap_SINETgenes[-2,],aes(x=UMAP_dimension1,y=UMAP_dimension2,col=`Molecular clusters`),inherit.aes = F,size=1.5)+
  geom_point(data=umap_SINETgenes[which(!str_detect(samples_SI,"p")),],aes(x=UMAP_dimension1,y=UMAP_dimension2,col=`Molecular clusters`,fill=`Molecular clusters`),size=2,stroke=1,shape=21,inherit.aes = F) +
  geom_point(data=umap_SINETgenes[which(str_detect(samples_SI,"p")),],aes(x=UMAP_dimension1,y=UMAP_dimension2,col=`Molecular clusters`),fill="white",size=2,stroke=1,shape=21,inherit.aes = F) +
  scale_color_manual( values=colors_org_SI, limits = names(colors_org_SI) ) + scale_fill_manual( values=colors_org_SI, limits = names(colors_org_SI) ) +
  xlab("UMAP dimension 1") + ylab("UMAP dimension 2") + 
  xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
  theme(legend.key.size = unit(10,"pt"), legend.title = element_text(size = 10,face = "bold"), 
        legend.text = element_text(size = 8) , legend.spacing = unit(0,'pt'), 
        legend.spacing.y = unit(1,"pt") )  #+

gg_UMAP_SINETgenes

ggsave("Fig_3D_UMAP_SINET_raw_16082021.png",gg_UMAP_SINETgenes,height = 3,width = 5)
ggsave("Fig_3D_UMAP_SINET_raw_16082021.svg",gg_UMAP_SINETgenes,height = 3,width = 5)
#ggsave("Fig_3D_UMAP_raw.pdf",gg_UMAP,height = 4,width = 6)

write_tsv(umap_SINETgenes,path = "Fig3D_draft_data_SINET.tsv")
