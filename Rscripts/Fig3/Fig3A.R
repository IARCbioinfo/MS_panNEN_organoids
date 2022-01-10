library(tidyverse)
library(reshape2)
library(ggbeeswarm)
library(readxl)
library(ggforce)

#define colors
colors_org = c(LNET2="#aade87ff",LNET6="#5fd38dff",LNET13="#16502dff",LNET14="#6f917cff",
               LNET5="#e6a73cff",LNET10="#ff9955ff",LNET15="#ffd42aff", LNET16 = "#ff6600ff", LNET18= "#d0742fff", 
               LNET19="#2aff80ff", 
               LNET20 = "#f6e62bff", 
               LCNEC3="#ff8080ff",LCNEC4="#d35f5fff", LCNEC23 = "#ff5555ff", 
               LCNEC11="#ff5599ff",PANEC1="#8d5fd3ff",
               SINET7="#2ad4ffff",SINET8="#80b3ffff",SINET9="#5f8dd3ff",SINET12="#5fbcd3ff", SINET21="#0066ffff", SINET22="#2c5aa0ff")

colors_types = c(SI= "#5bc2f0ff", Lung = "#9b9972ff", Pancreas = "#8d5fd3ff")

colors_grades = c(G1= "#58b873f9", G2 = "#ff9955ff", G3 = "#f0677dff", "G1/G2" = "#58c1f0ff")

# load expression data 
#organoids
load("/data/gcs/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/Robjects/gene_1pass.SE.rda")
gene_1pass.SE.organoids = gene_1pass.SE

#reference lung NEN data
load("/data/gcs/lungNENomics/work/voegelec/RNASEQ/ANALYSES-2/Carcinoids-GEO-30/out_RNAseq-transcript-nf-2.2/Robjects/gene_1pass.SE.rda")
gene_1pass.SE.CarcGeo = gene_1pass.SE
load("/data/gcs/lungNENomics/work/voegelec/RNASEQ/ANALYSES-2/Carcinoids-LCNEC-SCLC-158+52/out_RNAseq-transcript-nf-2.2/Robjects/gene_1pass.SE.rda")
gene_1pass.SE.GigaSc  = gene_1pass.SE
colnames(gene_1pass.SE.GigaSc)[colnames(gene_1pass.SE.GigaSc)=="S01103_T2"] = "S01103"
colnames(gene_1pass.SE.GigaSc) = str_replace(colnames(gene_1pass.SE.GigaSc),"AB","")
load("/data/gcs/lungNENomics/work/organoids/RNAseq/processing/SINET/quantification_SINET/Robjects/gene_1pass.SE.rda")
gene_1pass.SE.SINET  = gene_1pass.SE
load("/data/gcs/lungNENomics/work/organoids/RNAseq/processing/SINET/quantification_SINET2_23042021/Robjects/gene_1pass.SE.rda")
gene_1pass.SE.SINET2  = gene_1pass.SE

rm(gene_1pass.SE)
gc()

#metadata
Attributes = read.table(unzip('/data/gcs/lungNENomics/work/organoids/metadata/DRMetrics/data/Attributes.txt.zip')[1], sep = '\t', header = T)[,1:9] # Gigascience
Attributes_organoids = read_xlsx("/data/gcs/lungNENomics/work/organoids/metadata/Data_Sharing_organoids_11062021.xlsx",sheet = 1,skip = 1)

# getting gene data 
NEdiff.genes = c("CHGA", "NCAM1", "SYP")  # , "HES6", "DDC", "UCHL1", "CALCA" ,"GRP",) #  "NEUROG3"
ID.NEdiff.genes = sapply( NEdiff.genes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )
NElin.genes = c("ASCL1","INSM1","NEUROD1") #c("GAP43", "FEZ1", "CHRNA3", "DCX")
ID.NElin.genes = sapply( NElin.genes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )

org_LNETg1 = c("LNET13T","LNET13Tp1","LNET14T","LNET14Tp1","LNET19T","LNET19Tp2","LNET6T","LNET6Tp1")
org_LNETg2 = c("LNET5T","LNET5Tp2.2","LNET5Tp4","LNET5Tp7","LNET10T","LNET10Tp11","LNET10Tp4","LNET15M","LNET15Mp2","LNET16M","LNET16Mp1","LNET16T","LNET16Tp2",
               "LNET18Tp2","LNET20M","LNET20Mp2")
high_pur_samples = c("LNET13T","LNET13Tp1","LNET14T","LNET14Tp1","LNET6T","LNET6Tp1","LNET10T","LNET10Tp11","LNET10Tp4","LNET15M","LNET15Mp2","LNET16M","LNET16Mp1")

expr_genes = assay(cbind(gene_1pass.SE.organoids,gene_1pass.SE.GigaSc,gene_1pass.SE.CarcGeo,gene_1pass.SE.SINET,gene_1pass.SE.SINET2),
                        "abundance_TPM")[c(ID.NEdiff.genes,ID.NElin.genes),]
rownames(expr_genes) = c(NEdiff.genes,NElin.genes)
expr_genes.tib = bind_cols(Sample= colnames(expr_genes), as_tibble(t(expr_genes)) ) %>% pivot_longer(-Sample,values_to = "Expression") %>% 
  mutate(Type="SINET",Experiment = "Reference")
for(x in names(colors_org)) expr_genes.tib$Experiment[str_detect(expr_genes.tib$Sample,x)] = x

expr_genes.tib$Experiment = factor(expr_genes.tib$Experiment,levels = c(names(colors_org),"Reference"))
expr_genes.tib$Expression[expr_genes.tib$Expression<=0.01] = 0.01
expr_genes.tib$Type[sapply(Attributes$Sample_ID,function(x) which(expr_genes.tib$Sample==x))] = Attributes$Histopathology_simplified
expr_genes.tib$Grade = "G1/G2"
expr_genes.tib$Grade[expr_genes.tib$Type=="Carcinoid"] = NA
expr_genes.tib$Grade[expr_genes.tib$Type%in%c("Atypical","Supra_carcinoid")] = "G2"
expr_genes.tib$Grade[expr_genes.tib$Type%in%c("Typical")] = "G1"
expr_genes.tib$Grade[expr_genes.tib$Type%in%c("LCNEC","SCLC")] = "G3"

expr_genes.tib$Grade[str_detect(expr_genes.tib$Experiment,"LCNEC")] = "G3"
expr_genes.tib$Grade[str_detect(expr_genes.tib$Experiment,"PANEC")] = "G3"
expr_genes.tib$Grade[expr_genes.tib$Sample%in%org_LNETg1] = "G1"
expr_genes.tib$Grade[expr_genes.tib$Sample%in%org_LNETg2] = "G2"

expr_genes.tib$Type[expr_genes.tib$Type %in% c("Atypical","Typical","Carcinoid","Supra_carcinoid")] = "LNET"
expr_genes.tib$Type[str_detect(expr_genes.tib$Experiment,"LNET")] = "LNET"
expr_genes.tib$Type[str_detect(expr_genes.tib$Experiment,"LCNEC|PANEC")] = "LCNEC"

expr_genes.tib$Grade = factor(expr_genes.tib$Grade,levels=c("G1/G2","G1","G2","G3"))
expr_genes.tib$Type = factor(expr_genes.tib$Type,levels=c("LNET","LCNEC","SCLC","SINET"))
expr_genes.tib$TypeGrade = interaction(expr_genes.tib$Type,expr_genes.tib$Grade,drop = T)

expr_genes.tib$TypeGrade = factor(expr_genes.tib$TypeGrade,levels = c("LNET.G1", "LNET.G2", "LCNEC.G3", "SCLC.G3","SINET.G1/G2"))

expr_genes.tib$Gene_group = "NE differentiation"
expr_genes.tib$Gene_group[expr_genes.tib$name %in% NElin.genes] = "NE lineage TF"

## function for plotting trajectory 
traject_RNA <- function(dat=expr_genes.tib,samples,col=colors_org[2],xpos=1-0.30,xoff=-0.001,fills=c(col,"white"),lty=1){
  res = list(trajRNA = geom_curve(data = left_join(dat[dat$Sample==samples[1],],dat[dat$Sample==samples[2],],
                                                   by=c("Type","Experiment","name"),suffix=c("",".y")) %>% mutate(xpos=xpos+rep(0:2,2)), 
                                  aes(x =  xpos, y = Expression , xend = xpos+xoff,yend = Expression.y), size=0.5,colour = col,inherit.aes = T ,curvature = -0.2,lineend = "round") , 
             posRNA =  geom_ellipse(data = dat[sapply(samples[1], function(x) which(dat$Sample==x)),]%>% mutate(Sample=factor(Sample,levels=samples)) %>% 
                                      mutate(xpos=xpos+rep(0:2,2)), 
                                    aes(x0 = xpos, y0 = Expression,a=0.04,b=0.2,angle=0), inherit.aes = F, lty=lty, lwd= 0.5, 
                                    colour = col,fill=fills[1]) ,
             posRNA2 =  geom_ellipse(data = dat[sapply(samples[2], function(x) which(dat$Sample==x)),]%>% mutate(Sample=factor(Sample,levels=samples)) %>% 
                                       mutate(xpos=xpos+rep(0:2,2)+xoff), 
                                     aes(x0 = xpos, y0 = Expression,a=0.04,b=0.2,angle=0), inherit.aes = F, lty=lty, lwd= 0.5, 
                                     colour = col,fill=fills[2])
  )
  
  return(res)
}

# plot
data_summary <- function(x) {
  res <- quantile(x,c(0.25,0.5,0.75))
  names(res) = c("ymin","y","ymax")
  return(res)
}

#traj5.1.4 = traject_RNA(dat=expr_genes.tib,samples=c("LNET5T","LNET5Tp4"),col=colors_org["LNET5"],xpos = 1,xoff = 0.03)
#traj5.1.7 = traject_RNA(dat=expr_genes.tib,samples=c("LNET5Tp4","LNET5Tp7"),col=colors_org["LNET5"],xpos = 1+0.03,xoff = 0.03,fills = c("white","white"))
#traj5.2   = traject_RNA(dat=expr_genes.tib,samples=c("LNET5T","LNET5Tp2.2"),col=colors_org["LNET5"],xpos = 1,xoff = -0.03)
traj6 = traject_RNA(dat=expr_genes.tib,samples=c("LNET6T","LNET6Tp1"),col=colors_org["LNET6"],xpos = 1-0.2,xoff=0.01)
#traj19 = traject_RNA(dat=expr_genes.tib,samples=c("LNET19T","LNET19Tp2"),col=colors_org["LNET19"],xpos = 1-0.15,xoff=0.01)
traj13 = traject_RNA(dat=expr_genes.tib,samples=c("LNET13T","LNET13Tp1"),col=colors_org["LNET13"],xpos = 1+0,xoff=0.01)
traj14 = traject_RNA(dat=expr_genes.tib,samples=c("LNET14T","LNET14Tp1"),col=colors_org["LNET14"],xpos = 1+0.2,xoff=0.01)

traj15 = traject_RNA(dat=expr_genes.tib,samples=c("LNET15M","LNET15Mp2"),col=colors_org["LNET15"],xpos = 1-0.2,xoff=0.01)
traj10.4 = traject_RNA(dat=expr_genes.tib,samples=c("LNET10T","LNET10Tp4"),col=colors_org["LNET10"],xpos = 1-0.025,xoff=0.05)
traj10.11 = traject_RNA(dat=expr_genes.tib,samples=c("LNET10Tp4","LNET10Tp11"),col=colors_org["LNET10"],xpos = 1+0.025,xoff=0.05,fills = c("white","white"))
#traj20 = traject_RNA(dat=expr_genes.tib,samples=c("LNET20M","LNET20Mp2"),col=colors_org["LNET20"],xpos = 1,xoff=0.01)
traj16M = traject_RNA(dat=expr_genes.tib,samples=c("LNET16M","LNET16Mp1"),col=colors_org["LNET16"],xpos = 1+0.2,xoff=0.01,lty="12")
#traj16T = traject_RNA(dat=expr_genes.tib,samples=c("LNET16T","LNET16Tp2"),col=colors_org["LNET16"],xpos = 1+0.3,xoff=0.01)

# LCNEC
traj1.4  = traject_RNA(dat=expr_genes.tib,samples=c("PANEC1T","PANEC1Tp4"),col=colors_org["PANEC1"],xpos = 1-0.15,xoff=+0.01)
traj1.14 = traject_RNA(dat=expr_genes.tib,samples=c("PANEC1Tp4","PANEC1Tp14"),col=colors_org["PANEC1"],xpos = 1-0.15+0.01,xoff=+0.03,fills = c("white","white"))
traj3.17 = traject_RNA(dat=expr_genes.tib,samples=c("LCNEC3T","LCNEC3Tp17.2"),col=colors_org["LCNEC3"],xpos = 1,xoff=0.03)
traj3.24 = traject_RNA(dat=expr_genes.tib,samples=c("LCNEC3Tp17.2","LCNEC3Tp24"),col=colors_org["LCNEC3"],xpos = 1+0.03,xoff=0.03,fills = c("white","white"))
traj11   = traject_RNA(dat=expr_genes.tib,samples=c("LCNEC11M","LCNEC11Mp3"),col=colors_org["LCNEC11"],xpos = 1+0.15,xoff=-0.05,lty="12")
traj4.7  = traject_RNA(dat=expr_genes.tib,samples=c("LCNEC4T","LCNEC4Tp7"),col=colors_org["LCNEC4"],xpos = 1+0.25,xoff=0.03)
traj4.24 = traject_RNA(dat=expr_genes.tib,samples=c("LCNEC4Tp7","LCNEC4Tp24"),col=colors_org["LCNEC4"],xpos = 1+0.28,xoff=0.03,fills = c("white","white"))

# SINET
traj7  = traject_RNA(dat=expr_genes.tib,samples=c("SINET7M","SINET7Mp2"),col=colors_org["SINET7"],xpos = 1-0.3,xoff=0.01,lty="12")
traj8  = traject_RNA(dat=expr_genes.tib,samples=c("SINET8M","SINET8Mp2"),col=colors_org["SINET8"],xpos = 1-0.15,xoff=0.01,lty="12")
traj21 = traject_RNA(dat=expr_genes.tib,samples=c("SINET21M","SINET21Mp2"),col=colors_org["SINET21"],xpos = 1,xoff=0.01,lty="12")
traj12.1 = traject_RNA(dat=expr_genes.tib,samples=c("SINET12M","SINET12Mp1.1"),col=colors_org["SINET12"],xpos = 1.15,xoff=-0.04,lty="12")
traj12.3 = traject_RNA(dat=expr_genes.tib,samples=c("SINET12M","SINET12Mp1.3"),col=colors_org["SINET12"],xpos = 1.15,xoff=0.04,lty="12")
traj22 = traject_RNA(dat=expr_genes.tib,samples=c("SINET22M","SINET22Mp2"),col=colors_org["SINET22"],xpos = 1+0.3,xoff=-0.02,lty="12")

## plot
Fig3A_raw <- ggplot( expr_genes.tib %>% filter(Experiment=="Reference",Type!="SCLC",!is.na(Grade)) , aes(x=name,y=Expression,color=Experiment) ) + 
  geom_violin(aes(x=name,y=Expression,color=NA),draw_quantiles = c(0.25,0.5,0.75), color=NA,scale = "width",fill=rgb(0.85,0.85,0.85),width=1,position=position_dodge(width=1)) + 
  #geom_quasirandom(method = "quasirandom",size=0.5,cex = 1.5,priority = "density") + 
  geom_hline(yintercept =1,linetype="dashed")+
  #stat_summary(fun.data=data_summary,geom = "crossbar", width=0.5,col="white",cex=0.5) +
  traj6 + traj13 +traj14 + 
  #traj19 +  #traj5.2 + traj5.1.4 + traj5.1.7 +#traj16T + traj20 + # mixed tumors
  traj10.4 + traj10.11 + traj15 + traj16M + 
  traj3.17 + traj3.24 + traj4.7 + traj4.24 + traj11 +
  traj1.4 + traj1.14+
  traj7 + traj8 + traj12.1 + traj12.3 +traj21 +traj22 +
  geom_ellipse(data = expr_genes.tib %>% filter(Sample=="LCNEC23Mp3") %>% mutate(xpos=xpos+rep(0:2,2)+xoff), aes(x0 = xpos, y0 = Expression,a=0.04,b=0.2,angle=0), 
               inherit.aes = F, lty=lty, lwd=0.5,colour = colors_org["LCNEC23"],fill="white") +
  theme_minimal() + theme(legend.title=element_blank() ,axis.text.x = element_text(face="italic",angle = 90, vjust = 0.5, hjust=1), 
                          strip.text.x =element_text(face = "bold",size = 10),axis.title.x = element_blank() , panel.grid.major.x = element_blank() ) + 
  facet_grid(TypeGrade~Gene_group,scales = "free_x", space = "free",drop = T) + 
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  labs(x="Markers",y="Gene expression (TPM)" ) + 
  scale_y_log10(breaks=c(0.01,1,10**2,10**4),limits=c(0.01,200000),labels=c("≤0.01",1,100,10000)) + 
  scale_color_manual(values=c(colors_org, Reference=rgb(0.75,0.75,0.75))[levels(droplevels(expr_genes.tib$Experiment))],
                     limits = levels(droplevels(expr_genes.tib$Experiment)), 
                     labels = levels(droplevels(expr_genes.tib$Experiment)) ) #+ coord_fixed()

Fig3A_raw

ggsave("/data/gcs/lungNENomics/work/organoids/figures/Fig3A_raw_19112021.svg",Fig3A_raw,width = 3.8,height = 3.8*842/783)
ggsave("/data/gcs/lungNENomics/work/organoids/figures/Fig3A_raw_19112021.png",Fig3A_raw,width = 3.8,height = 3.9*842/783)
