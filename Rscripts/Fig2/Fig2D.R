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

# getting MKI67 data 
ID.MKI67 = which(rowData(gene_1pass.SE.organoids)$gene_name=="MKI67")

org_LNETg1 = c("LNET13T","LNET13Tp1","LNET14T","LNET14Tp1","LNET19T","LNET19Tp2","LNET6T","LNET6Tp1")
org_LNETg2 = c("LNET5T","LNET5Tp2.2","LNET5Tp4","LNET5Tp7","LNET10T","LNET10Tp11","LNET10Tp4","LNET15M","LNET15Mp2","LNET16M","LNET16Mp1","LNET16T","LNET16Tp2",
               "LNET18Tp2","LNET20M","LNET20Mp2")

expr_genes_KI67 = assay(cbind(gene_1pass.SE.organoids,gene_1pass.SE.GigaSc,gene_1pass.SE.CarcGeo,gene_1pass.SE.SINET,gene_1pass.SE.SINET2),
                         "abundance_TPM")[ID.MKI67,]
expr_genes_KI67.tib = tibble(Sample= names(expr_genes_KI67), Expression=expr_genes_KI67, Type="SINET",Experiment = "Reference")
for(x in names(colors_org)) expr_genes_KI67.tib$Experiment[str_detect(expr_genes_KI67.tib$Sample,x)] = x

expr_genes_KI67.tib$Experiment = factor(expr_genes_KI67.tib$Experiment,levels = c(names(colors_org),"Reference"))
expr_genes_KI67.tib$Expression[expr_genes_KI67.tib$Expression<=0.01] = 0.01
expr_genes_KI67.tib$Type[sapply(Attributes$Sample_ID,function(x) which(expr_genes_KI67.tib$Sample==x))] = Attributes$Histopathology_simplified
expr_genes_KI67.tib$Grade = "G1/G2"
expr_genes_KI67.tib$Grade[expr_genes_KI67.tib$Type=="Carcinoid"] = NA
expr_genes_KI67.tib$Grade[expr_genes_KI67.tib$Type%in%c("Atypical","Supra_carcinoid")] = "G2"
expr_genes_KI67.tib$Grade[expr_genes_KI67.tib$Type%in%c("Typical")] = "G1"
expr_genes_KI67.tib$Grade[expr_genes_KI67.tib$Type%in%c("LCNEC","SCLC")] = "G3"

#expr_genes_KI67.tib$Grade = factor(expr_genes_KI67.tib$Type,levels=c("G1","G2","G3","SCLC","SINET"))
expr_genes_KI67.tib$Grade[str_detect(expr_genes_KI67.tib$Experiment,"LCNEC")] = "G3"
expr_genes_KI67.tib$Grade[str_detect(expr_genes_KI67.tib$Experiment,"PANEC")] = "G3"
expr_genes_KI67.tib$Grade[expr_genes_KI67.tib$Sample%in%org_LNETg1] = "G1"
expr_genes_KI67.tib$Grade[expr_genes_KI67.tib$Sample%in%org_LNETg2] = "G2"

expr_genes_KI67.tib$Type[expr_genes_KI67.tib$Type %in% c("Atypical","Typical","Carcinoid","Supra_carcinoid")] = "LNET"
expr_genes_KI67.tib$Type[str_detect(expr_genes_KI67.tib$Experiment,"LNET")] = "LNET"
expr_genes_KI67.tib$Type[str_detect(expr_genes_KI67.tib$Experiment,"LCNEC|PANEC")] = "LCNEC"

expr_genes_KI67.tib$Grade = factor(expr_genes_KI67.tib$Grade,levels=c("G1/G2","G1","G2","G3"))

expr_genes_KI67.tib$Type = factor(expr_genes_KI67.tib$Type,levels=c("LNET","LCNEC","SCLC","SINET"))

## alternative: comparisons of distributions 
library(ggridges)
Fig2D_alt <- ggplot( expr_genes_KI67.tib %>% filter(Experiment=="Reference",Type!="SCLC",!is.na(Grade)) , aes(y=Grade,x=Expression,fill=Grade) ) + 
  geom_density_ridges2(scale = 1.5,col=NA) + 
  geom_point(data = expr_genes_KI67.tib %>% filter(Experiment!="Reference",str_detect(Sample,"p[0-9.]*$")) , size=4, pch=16, col="white",fill="black") +
  geom_point(data = expr_genes_KI67.tib %>% filter(Experiment!="Reference",str_detect(Sample,"p[0-9.]*$")) , size=2.7, pch=16, col="black",fill="black") + 
  theme_classic()  + labs(y="Histopathological type",x=expression(italic(MKI67)~" Expression (TPM)") ) + 
  geom_vline(xintercept = 1,linetype="dashed") + 
  scale_x_log10(breaks=c(0.01,0.1,1,10,10**2),limits=c(0.01,200),labels=c("≤0.01",0.1,1,10,100)) + scale_fill_manual(values=alpha(colors_grades,0.5))
  
#  theme_joy()

ggsave("/data/gcs/lungNENomics/work/organoids/figures/Fig2D_raw_17112021.svg",Fig2D_alt,width = 3.5,height = 3.5)


## function for plotting trajectory 
traject_RNA_1gene <- function(dat=expr_genes_NET,samples,col=colors_org[3],xpos=1-0.30,xoff=-0.001,fills=c(col,"white"),lty=1){
  res = list(trajRNA = geom_curve(data = left_join(dat[dat$Sample==samples[1],],dat[dat$Sample==samples[2],],
                                                   by=c("Type","Experiment"),suffix=c("",".y"))%>% arrange(Sample) %>% mutate(xpos=xpos), 
                                  aes(x =  xpos, y = Expression , xend = xpos+xoff,yend = Expression.y), size=0.5,colour = col,inherit.aes = T ,curvature = -0.2,lineend = "round") , 
               posRNA =  geom_ellipse(data = dat[sapply(samples[1], function(x) which(dat$Sample==x)),]%>% mutate(Sample=factor(Sample,levels=samples)) %>% 
                                       arrange(Sample) %>% 
                                       mutate(xpos=xpos), 
                                     aes(x0 = xpos, y0 = Expression,a=0.05,b=0.08,angle=0), inherit.aes = F, lty=lty, lwd= 1, 
                                     colour = col,fill=fills[1]) ,
             posRNA2 =  geom_ellipse(data = dat[sapply(samples[2], function(x) which(dat$Sample==x)),]%>% mutate(Sample=factor(Sample,levels=samples)) %>% 
                             arrange(Sample) %>% 
                             mutate(xpos=xpos+xoff), 
                           aes(x0 = xpos, y0 = Expression,a=0.05,b=0.08,angle=0), inherit.aes = F, lty=lty, lwd= 1, 
                           colour = col,fill=fills[2])
             )
               
               #geom_point(data = dat[sapply(samples, function(x) which(dat$Sample==x)),]%>% mutate(Sample=factor(Sample,levels=samples)) %>% 
              #                     arrange(Sample) %>% 
              #                     mutate(xpos=c(xpos,xpos+xoff)), 
              #                   aes(x = xpos, y = Expression), size=2,stroke=1,colour = col,fill=rep(fills,each=length(xpos)),shape=21,inherit.aes = T) )
  return(res)
}

# plot
data_summary <- function(x) {
  res <- quantile(x,c(0.25,0.5,0.75))
  names(res) = c("ymin","y","ymax")
  return(res)
}

traj5.1.4 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET5T","LNET5Tp4"),col=colors_org["LNET5"],xpos = 1,xoff = 0.03)
traj5.1.7 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET5Tp4","LNET5Tp7"),col=colors_org["LNET5"],xpos = 1+0.03,xoff = 0.03,fills = c("white","white"))
traj5.2   = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET5T","LNET5Tp2.2"),col=colors_org["LNET5"],xpos = 1,xoff = -0.03)

traj6 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET6T","LNET6Tp1"),col=colors_org["LNET6"],xpos = 1-0.3,xoff=0.01)
traj19 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET19T","LNET19Tp2"),col=colors_org["LNET19"],xpos = 1-0.15,xoff=0.01)
traj14 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET14T","LNET14Tp1"),col=colors_org["LNET14"],xpos = 1+0.15,xoff=0.01)
traj13 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET13T","LNET13Tp1"),col=colors_org["LNET13"],xpos = 1+0.3,xoff=0.01)

traj15 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET15M","LNET15Mp2"),col=colors_org["LNET15"],xpos = 2-0.2,xoff=0.01)
traj10.4 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET10T","LNET10Tp4"),col=colors_org["LNET10"],xpos = 2-0.1,xoff=0.03)
traj10.11 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET10Tp4","LNET10Tp11"),col=colors_org["LNET10"],xpos = 2-0.07,xoff=0.03,fills = c("white","white"))
traj20 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET20M","LNET20Mp2"),col=colors_org["LNET20"],xpos = 2,xoff=0.01)
traj16M = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET16M","LNET16Mp1"),col=colors_org["LNET16"],xpos = 2+0.15,xoff=0.01,lty="12")
traj16T = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LNET16T","LNET16Tp2"),col=colors_org["LNET16"],xpos = 2+0.3,xoff=0.01)

# LCNEC
traj3.17 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LCNEC3T","LCNEC3Tp17.2"),col=colors_org["LCNEC3"],xpos = 1-0.25,xoff=0.03)
traj3.24 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LCNEC3Tp17.2","LCNEC3Tp24"),col=colors_org["LCNEC3"],xpos = 1-0.22,xoff=0.03,fills = c("white","white"))
traj1.4  = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("PANEC1T","PANEC1Tp4"),col=colors_org["PANEC1"],xpos = 1,xoff=-0.01)
traj1.14 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("PANEC1Tp4","PANEC1Tp14"),col=colors_org["PANEC1"],xpos = 1-0.01,xoff=-0.03,fills = c("white","white"))
traj11   = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LCNEC11M","LCNEC11Mp3"),col=colors_org["LCNEC11"],xpos = 1+0.15,xoff=-0.05,lty="12")
traj4.7  = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LCNEC4T","LCNEC4Tp7"),col=colors_org["LCNEC4"],xpos = 1+0.25,xoff=0.03)
traj4.24 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("LCNEC4Tp7","LCNEC4Tp24"),col=colors_org["LCNEC4"],xpos = 1+0.28,xoff=0.03,fills = c("white","white"))

# SINET
traj7  = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("SINET7M","SINET7Mp2"),col=colors_org["SINET7"],xpos = 1-0.3,xoff=0.01,lty="12")
traj8  = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("SINET8M","SINET8Mp2"),col=colors_org["SINET8"],xpos = 1-0.15,xoff=0.01,lty="12")
traj12.1 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("SINET12M","SINET12Mp1.1"),col=colors_org["SINET12"],xpos = 1,xoff=-0.04,lty="12")
traj12.3 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("SINET12M","SINET12Mp1.3"),col=colors_org["SINET12"],xpos = 1,xoff=0.04,lty="12")
traj22 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("SINET22M","SINET22Mp2"),col=colors_org["SINET22"],xpos = 1+0.15,xoff=-0.02,lty="12")
traj21 = traject_RNA_1gene(dat=expr_genes_KI67.tib,samples=c("SINET21M","SINET21Mp2"),col=colors_org["SINET21"],xpos = 1+0.3,xoff=0.01,lty="12")

## plot
Fig2D_raw <- ggplot( expr_genes_KI67.tib %>% filter(Experiment=="Reference",Type!="SCLC",!is.na(Grade)) , aes(x=Grade,y=Expression,color=Experiment) ) + 
  geom_violin(aes(x=Grade,y=Expression,color=NA),draw_quantiles = c(0.25,0.5,0.75), color=NA,scale = "width",fill=rgb(0.85,0.85,0.85),width=1,position=position_dodge(width=1)) + 
  geom_quasirandom(method = "quasirandom",size=0.5,cex = 1.5,priority = "density") + 
  geom_hline(yintercept =1,linetype="dashed")+
  stat_summary(fun.data=data_summary,geom = "crossbar", width=0.5,col="white",cex=0.5) +
  traj6 + traj13 +traj14 + 
  traj19 + 
  traj10.4 + traj10.11 + traj15 + traj16M + 
  traj5.2 + traj5.1.4 + traj5.1.7 +
  traj16T + traj20 +
  traj3.17 + traj3.24 + traj4.7 + traj4.24 + traj11 +
  traj1.4 + traj1.14+
  traj7 + traj8 + traj12.1 + traj12.3 +traj21 +traj22 +
  geom_ellipse(data = dat %>% filter(Sample=="LNET18Tp2") %>% mutate(xpos=xpos+xoff), aes(x0 = 2-0.3, y0 = Expression,a=0.05,b=0.08,angle=0), inherit.aes = F, lty=lty, lwd= 1, 
                colour = colors_org["LNET18"],fill="white") +
  geom_ellipse(data = dat %>% filter(Sample=="LCNEC23Mp3") %>% mutate(xpos=xpos+xoff), aes(x0 = 1-0.40, y0 = Expression,a=0.05,b=0.08,angle=0), inherit.aes = F, lty=lty, lwd= 1, 
               colour = colors_org["LCNEC23"],fill="white") +
  theme_minimal() + theme(legend.title=element_blank() ,axis.text.x = element_text(face="italic",angle = 90, vjust = 0.5, hjust=1), 
                          strip.text.x =element_text(face = "bold",size = 10),axis.title.x = element_blank() , panel.grid.major.x = element_blank() ) + 
  facet_grid(~Type,scales = "free_x", space = "free",drop = T) + 
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  labs(x="Markers",y=expression(italic(MKI67)~" Expression (TPM)") ) + 
  scale_y_log10(breaks=c(0.01,0.1,1,10,10**2),limits=c(0.01,200),labels=c("≤0.01",0.1,1,10,100)) + 
  scale_color_manual(values=c(colors_org, Reference=rgb(0.75,0.75,0.75))[levels(droplevels(expr_genes_KI67.tib$Experiment))],
                     limits = levels(droplevels(expr_genes_KI67.tib$Experiment)), 
                     labels = levels(droplevels(expr_genes_KI67.tib$Experiment)) ) #+ coord_fixed()

Fig2D_raw

ggsave("/data/gcs/lungNENomics/work/organoids/figures/Fig2D_raw_02112021.svg",Fig2D_raw,width = 6/0.94,height = 3)
ggsave("/data/gcs/lungNENomics/work/organoids/figures/Fig2D_raw_02112021.png",Fig2D_raw,width = 6/0.94,height = 3)

# compute a few stats
expr_genes_KI67.tib %>% group_by(Experiment=="Reference",Type,Grade) %>% summarize(mean(Expression))

expr_genes_KI67.tib %>% filter(str_detect(Sample,"p[0-9]")) %>% group_by(Experiment=="Reference",Type,Grade) %>% summarize(mean(Expression))

write_tsv(Fig2D_raw$data,"/data/gcs/lungNENomics/work/organoids/figures/Fig2D_raw_02112021.tsv")


## Figure 3
Fig3Bg2_raw <- ggplot( expr_genes_NETg1[expr_genes_NETg1$Experiment=="Reference",] , aes(x=ID,y=value,color=Experiment) ) + 
  geom_violin(aes(x=ID,y=value,color=NA),draw_quantiles = c(0.25,0.5,0.75), color=NA,scale = "width",fill=rgb(0.85,0.85,0.85)) + 
  geom_quasirandom(method = "quasirandom",size=0.5,cex = 1.5,priority = "density") + 
  geom_hline(yintercept =1,linetype="dashed")+
  stat_summary(fun.data=data_summary,geom = "crossbar", width=0.25,col="white",cex=0.5) +
  traj10.4[[1]]+traj10.4[[2]]+traj10.11[[1]]+traj10.11[[2]]+ 
  traj15[[1]]+traj15[[2]]+ traj16T[[1]]+traj16T[[2]]+ traj16M[[1]]+traj16M[[2]]+ traj20[[1]]+ traj20[[2]]+
  geom_point(data = expr_genes_NETg2 %>% filter(Sample=="LNET18Tp2") %>% arrange(Set,ID) %>% mutate(xpos=c(1:3,1:3,1:4)+0.18), 
             aes(x = xpos, y = value), size=2,stroke=1,colour = colors_org[9],fill="white",shape=21,inherit.aes = T) + 
  theme_minimal() + theme(legend.title=element_blank() ,axis.text.x = element_text(face="italic",angle = 90, vjust = 0.5, hjust=1), 
                          strip.text.x =element_text(face = "bold",size = 10),axis.title.x = element_blank() , panel.grid.major.x = element_blank() ) + 
  facet_grid(~Set,scales = "free", space = "free",drop = T) + 
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  labs(x="Markers",y="Expression (TPM)") + 
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,10**2,10**3,10**4,10**5),limits=c(0.001,160000),labels=c("≤0.001",0.01,0.1,1,10,100,"1,000","10,000","100,000")) + 
  scale_color_manual(values=c(colors_org[levels(expr_genes_NETg2$Experiment)[-6]],rgb(0.7,0.7,0.7)),limits = levels(expr_genes_NETg2$Experiment), labels = levels(expr_genes_NETg2$Experiment) ) 

Fig3Bg2_raw

ggsave("Fig3Bg2_raw_21072021.svg",Fig3Bg2_raw,width = 3*3,height = 3)
ggsave("Fig3Bg2_raw_21072021.png",Fig3Bg2_raw,width = 3*3,height = 3)

#Fig. 3C SINET
expr_genes_SINET = assay(cbind(gene_1pass.SE.organoids[,str_detect(colnames(gene_1pass.SE.organoids),"SINET")],gene_1pass.SE.SINET,gene_1pass.SE.SINET2),
                         "abundance_TPM")[c(ID.NEdiff.genes,ID.NElin.genes,ID.clin.genes,ID.LNET.genes),]
rownames(expr_genes_SINET) = c(NEdiff.genes,NElin.genes,clinical.genes,LNET.genes )
expr_genes_SINET = melt(expr_genes_SINET,varnames = c("ID","Sample"))
expr_genes_SINET$Set = "NE differentiation"
expr_genes_SINET$Set[expr_genes_SINET$ID%in%NElin.genes ]  = "NE lineage TF"
expr_genes_SINET$Set[expr_genes_SINET$ID%in%clinical.genes ]  = "Clinically relevant"
expr_genes_SINET$Set[expr_genes_SINET$ID%in%c(LNET.genes) ]  = "Clinically relevant"
expr_genes_SINET$Set = factor(expr_genes_SINET$Set,levels = unique(expr_genes_SINET$Set))

#expr_genes_SINET$ID = factor(expr_genes_SINET$ID, levels=names( sort(geom_strat_mean,decreasing = T) ) )#sort(c( PNEN.genes,Neuro.genes,EGF.genes,TRK.genes,clinical.genes,LNET.genes,Other.genes,Request.genes )) )
expr_genes_SINET$Experiment = "Reference"
for(x in names(colors_org)) expr_genes_SINET$Experiment[str_detect(expr_genes_SINET$Sample,x)] = x
expr_genes_SINET$Experiment = factor(expr_genes_SINET$Experiment,levels = c("SINET7","SINET8","SINET12","SINET21","SINET22","Reference"))
expr_genes_SINET$value[expr_genes_SINET$value==0] = 0.001


traj7 = traject_RNA(dat=expr_genes_SINET,samples=c("SINET7M","SINET7Mp2"),col=colors_org["SINET7"],xpos = c(1:3,1:3,1:4)-0.3,geneset = levels(expr_genes_SINET$Set))
traj8 = traject_RNA(dat=expr_genes_SINET,samples=c("SINET8M","SINET8Mp2"),col=colors_org["SINET8"],xpos = c(1:3,1:3,1:4)-0.15,geneset = levels(expr_genes_SINET$Set))
traj21 = traject_RNA(dat=expr_genes_SINET,samples=c("SINET21M","SINET21Mp2"),col=colors_org["SINET21"],xpos = c(1:3,1:3,1:4)+0.15,geneset = levels(expr_genes_SINET$Set))
traj22 = traject_RNA(dat=expr_genes_SINET,samples=c("SINET22M","SINET22Mp2"),col=colors_org["SINET22"],xpos = c(1:3,1:3,1:4)+0.3,geneset = levels(expr_genes_SINET$Set))

traj12.1 = traject_RNA(dat=expr_genes_SINET,samples=c("SINET12M","SINET12Mp1.1"),col=colors_org["SINET12"],xpos = c(1:3,1:3,1:4)-0,xoff=-0.05,geneset = levels(expr_genes_NEN$Set))
traj12.2 =  traject_RNA(dat=expr_genes_SINET,samples=c("SINET12M","SINET12Mp1.3"),col=colors_org["SINET12"],xpos = c(1:3,1:3,1:4)-0,xoff=0.05,geneset = levels(expr_genes_NEN$Set),fills=c("white","white"))

Fig3D_raw <- ggplot( expr_genes_SINET[expr_genes_SINET$Experiment=="Reference",] , aes(x=ID,y=value,color=Experiment) ) + 
  geom_violin(aes(x=ID,y=value,color=NA),draw_quantiles = c(0.25,0.5,0.75), color=NA,scale = "width",fill=rgb(0.85,0.85,0.85)) + 
  geom_quasirandom(method = "quasirandom",size=0.5,cex = 1.5,priority = "density") + 
  geom_hline(yintercept =1,linetype="dashed")+
  stat_summary(fun.data=data_summary,geom = "crossbar", width=0.25,col="white",cex=0.5) +
  traj7[[1]] + traj7[[2]]+traj8[[1]] + traj8[[2]]+traj12.1[[1]] +traj12.1[[2]] +traj12.2[[1]] + traj12.2[[2]]+
  traj21[[1]] + traj21[[2]]+traj22[[1]] + traj22[[2]]+
  theme_minimal() + theme(legend.title=element_blank() ,axis.text.x = element_text(face="italic",angle = 90, vjust = 0.5, hjust=1), 
                          strip.text.x =element_text(face = "bold",size = 10),axis.title.x = element_blank() , panel.grid.major.x = element_blank() ) + 
  facet_grid(~Set,scales = "free", space = "free",drop = T) + 
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  labs(x="Markers",y="Expression (TPM)") + 
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,10**2,10**3,10**4,10**5),limits=c(0.001,160000),labels=c("≤0.001",0.01,0.1,1,10,100,"1,000","10,000","100,000")) + 
  scale_color_manual(values=c(colors_org[c(17:18,20:22)],rgb(0.7,0.7,0.7)),limits = levels(expr_genes_SINET$Experiment), labels = levels(expr_genes_SINET$Experiment) ) 

Fig3D_raw

ggsave("Fig3Cg1_raw_21072021.svg",Fig3B_raw,width = 3*3,height = 3)
ggsave("Fig3Cg1_raw_21072021.png",Fig3B_raw,width = 3*3,height = 3)

# Fig. 3B :  Expression of canonical NE marker genes for LCNEC organoids and ref panel
NEdiff.genes = c("CHGA", "NCAM1", "SYP")  # , "HES6", "DDC", "UCHL1", "CALCA" ,"GRP",) #  "NEUROG3"
ID.NEdiff.genes = sapply( NEdiff.genes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )
NElin.genes = c("ASCL1","INSM1","NEUROD1") #c("GAP43", "FEZ1", "CHRNA3", "DCX")
ID.NElin.genes = sapply( NElin.genes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )
#Other.genes = c("FLT1","KDR", "FLT4","CD274","CTLA4", "POU2F3","YAP1")
#ID.Other.genes = sapply( Other.genes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )
#EGF.genes = c("EGFR","ERBB2", "ERBB3", "ERBB4", "EGF", "NRG1", "TGFA", "HBEGF", "AREG", "EPGN", "NRG2", "NRG3", "NRG4", "TNF") #TNF=TNFA?
#ID.EGF.genes = sapply( EGF.genes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )
#TRK.genes = c("NTRK1", "NTRK2", "NTRK3", "NGF", "NTF3", "BDNF", "NTF4")
#ID.TRK.genes = sapply( TRK.genes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )
clinical.genes = c("MKI67","SSTR2") #c("MGMT", "ASIC5", "CFC1", "HOXB13","IL13RA2","INSL5") #ASIC5 = ACCN5?
ID.clin.genes = sapply( clinical.genes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )
#SINET.genes = "SSTR2"
#ID.SINET.genes = sapply( SINET.genes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )
LNET.genes = c("DLL3","OTP")
ID.LNET.genes = sapply( LNET.genes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )

#Request.genes = c("CD200")
#ID.request.genes = sapply( Request.genes, function(x) which(rowData(gene_1pass.SE.organoids)$gene_name==x) )

expr_genes_NEN = assay(cbind(gene_1pass.SE.organoids[,str_detect(colnames(gene_1pass.SE.organoids),"NEC")],
                             gene_1pass.SE.GigaSc[,colnames(gene_1pass.SE.GigaSc)%in%Attributes$Sample_ID[str_detect(Attributes$Histopathology,"LCNEC")]]),
                       "abundance_TPM")[c(ID.NEdiff.genes,ID.NElin.genes,ID.clin.genes,ID.LNET.genes),] #c(ID.PNEN.genes,ID.Neuro.genes,ID.EGF.genes,ID.TRK.genes,ID.clin.genes,ID.LNET.genes,ID.Other.genes,ID.request.genes),]
rownames(expr_genes_NEN) = c(NEdiff.genes,NElin.genes,clinical.genes,LNET.genes )#PNEN.genes,Neuro.genes,EGF.genes,TRK.genes,clinical.genes,LNET.genes,Other.genes,Request.genes)
expr_genes_NEN = melt(expr_genes_NEN,varnames = c("ID","Sample"))
expr_genes_NEN$Set = "NE differentiation"
expr_genes_NEN$Set[expr_genes_NEN$ID%in%NElin.genes ]  = "NE lineage TF"
#expr_genes_NEN$Set[expr_genes_NEN$ID%in%EGF.genes ]  = "EGF signalling"
#expr_genes_NEN$Set[expr_genes_NEN$ID%in%TRK.genes ]  = "TRK signalling"
expr_genes_NEN$Set[expr_genes_NEN$ID%in%clinical.genes ]  = "Clinically relevant"
expr_genes_NEN$Set[expr_genes_NEN$ID%in%c(LNET.genes) ]  = "Clinically relevant"
#expr_genes_NEN$Set[expr_genes_NEN$ID%in%Request.genes ] = "Special request"
#expr_genes_NEN$Set[expr_genes_NEN$ID%in%Other.genes ] = "Other"
expr_genes_NEN$Set = factor(expr_genes_NEN$Set,levels = unique(expr_genes_NEN$Set))

# order
#expr_genes_all = assay(cbind(gene_1pass.SE.GigaSc[,colnames(gene_1pass.SE.GigaSc)%in%Attributes$Sample_ID],gene_1pass.SE.SINET),
#      "abundance_TPM")[c(ID.PNEN.genes,ID.Neuro.genes,ID.EGF.genes,ID.TRK.genes,ID.clin.genes,ID.LNET.genes,ID.Other.genes,ID.request.genes),]
#rownames(expr_genes_all) = c(PNEN.genes,Neuro.genes,EGF.genes,TRK.genes,clinical.genes,LNET.genes,Other.genes,Request.genes)
#meanexpr = cbind( LCNEC = rowMeans(expr_genes_all[,colnames(expr_genes_all)%in%Attributes$Sample_ID[str_detect(Attributes$Histopathology,"LCNEC")]]),
#       LNET  = rowMeans(expr_genes_all[,colnames(expr_genes_all)%in%Attributes$Sample_ID[str_detect(Attributes$Histopathology,"Atypical|Carcinoid|Typical")]]),
#       SINET = rowMeans(expr_genes_all[,colnames(expr_genes_all)%in%colnames(gene_1pass.SE.SINET)]) )
#apply(meanexpr,1,function(x) 3/sum(1/x) )
#geom_strat_mean = apply(meanexpr,1,function(x) prod(x)**(1/3) ) # geometric mean - to manage different orders of magnitude

#expr_genes_NEN$ID  = factor(expr_genes_NEN$ID, levels=names( sort(geom_strat_mean,decreasing = T) ) )
#sort(c( PNEN.genes,Neuro.genes,EGF.genes,TRK.genes,clinical.genes,LNET.genes,Other.genes )) )

expr_genes_NEN$Experiment = "Reference"
for(x in names(colors_org)) expr_genes_NEN$Experiment[str_detect(expr_genes_NEN$Sample,x)] = x
expr_genes_NEN$Experiment = factor(expr_genes_NEN$Experiment,levels = c("LCNEC3","LCNEC4","LCNEC11","LCNEC23","PANEC1","Reference"))
expr_genes_NEN$value[expr_genes_NEN$value==0] = 0.001

data_summary <- function(x) {
  res <- quantile(x,c(0.25,0.5,0.75))
  names(res) = c("ymin","y","ymax")
  return(res)
}

traj3 = traject_RNA(dat=expr_genes_NEN,samples=c("LCNEC3T","LCNEC3Tp17.2"),col=colors_org[12],xpos = c(1:3,1:3,1:4)-0.25,geneset = levels(expr_genes_NEN$Set))
traj4.1 = traject_RNA(dat=expr_genes_NEN,samples=c("LCNEC4T","LCNEC4Tp7"),col=colors_org[13],xpos = c(1:3,1:3,1:4)-0.15,xoff=0.05,geneset = levels(expr_genes_NEN$Set))
traj4.2 =  traject_RNA(dat=expr_genes_NEN,samples=c("LCNEC4Tp7","LCNEC4Tp24"),col=colors_org[13],xpos = c(1:3,1:3,1:4)-0.1,xoff=0.05,geneset = levels(expr_genes_NEN$Set),fills=c("white","white"))
traj11 = traject_RNA(dat=expr_genes_NEN,samples=c("LCNEC11M","LCNEC11Mp3"),col=colors_org[15],xpos = c(1:3,1:3,1:4)+0.05,xoff=0.05,geneset = levels(expr_genes_NEN$Set))
traj1.1 = traject_RNA(dat=expr_genes_NEN,samples=c("PANEC1T","PANEC1Tp4"),col=colors_org[16],xpos = c(1:3,1:3,1:4)+0.3,xoff=+0.05,geneset = levels(expr_genes_NEN$Set))
traj1.2 = traject_RNA(dat=expr_genes_NEN,samples=c("PANEC1Tp4","PANEC1Tp14"),col=colors_org[16],xpos = c(1:3,1:3,1:4)+0.35,xoff=+0.05,geneset = levels(expr_genes_NEN$Set),fills=c("white","white"))

Fig3Cg1_raw <- ggplot( expr_genes_NEN[expr_genes_NEN$Experiment=="Reference",] , aes(x=ID,y=value,color=Experiment) ) + 
  geom_violin(aes(x=ID,y=value,color=NA),draw_quantiles = c(0.25,0.5,0.75), color=NA,scale = "width",fill=rgb(0.85,0.85,0.85)) + 
  geom_quasirandom(method = "quasirandom",size=0.5,cex = 1.5,priority = "density") + 
  geom_hline(yintercept =1,linetype="dashed")+
  stat_summary(fun.data=data_summary,geom = "crossbar", width=0.25,col="white",cex=0.5) +
  traj3[[1]] + traj3[[2]]+traj4.1[[1]] + traj4.1[[2]]+traj4.2[[1]] + traj4.2[[2]]+
  traj11[[1]] + traj11[[2]]+traj1.1[[1]] + traj1.1[[2]]+traj1.2[[1]] + traj1.2[[2]]+ #traj23[[2]] + 
  geom_point(data = expr_genes_NEN %>% filter(Sample=="LCNEC23Mp3") %>% arrange(Set,ID) %>% mutate(xpos=c(1:3,1:3,1:4)+0.2), 
             aes(x = xpos, y = value), size=2,stroke=1,colour = colors_org[14],fill="white",shape=21,inherit.aes = T) + 
  theme_minimal() + theme(legend.title=element_blank() ,axis.text.x = element_text(face="italic",angle = 90, vjust = 0.5, hjust=1), 
                          strip.text.x =element_text(face = "bold",size = 10),axis.title.x = element_blank() , panel.grid.major.x = element_blank() ) + 
  facet_grid(~Set,scales = "free", space = "free",drop = T) + 
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  labs(x="Markers",y="Expression (TPM)") + 
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,10**2,10**3,10**4,10**5),limits=c(0.001,160000),labels=c("≤0.001",0.01,0.1,1,10,100,"1,000","10,000","100,000")) + 
  scale_color_manual(values=c(colors_org[12:16],rgb(0.7,0.7,0.7)),limits = levels(expr_genes_NEN$Experiment), labels = levels(expr_genes_NEN$Experiment) ) 

Fig3Cg1_raw

ggsave("Fig3Cg1_raw_21072021.svg",Fig3B_raw,width = 3*3,height = 3)
ggsave("Fig3Cg1_raw_21072021.png",Fig3B_raw,width = 3*3,height = 3)


# combine panels 
library(patchwork)

ggsave("Fig3ABCD_raw_16082021.svg",(Fig3A_raw+theme(axis.text.x = element_blank(),axis.title.y = element_text(size=12),axis.text.y = element_text(size=11),
                                                    strip.text.x = element_text(size=13)) )/
         (Fig3Bg2_raw+theme(axis.text.x = element_blank(),strip.text.x = element_blank(),axis.title.y = element_text(size=12),axis.text.y = element_text(size=11)) )/
         (Fig3B_raw+theme(axis.text.x = element_blank(),strip.text.x = element_blank(),axis.title.y = element_text(size=12),axis.text.y = element_text(size=11)) )/
         (Fig3D_raw+theme(axis.text.x =element_text(size=15),strip.text.x = element_blank(),axis.title.y = element_text(size=12),
                          axis.title.x = element_blank(),
                          axis.text.y = element_text(size=11)) )+ plot_annotation(tag_levels = 'A')&theme(plot.tag = element_text(size=20,face = "bold"),
                                                                                                          plot.margin = unit(c(0,0,0,0),"pt")),
       width = 2.5*3.5,height = 2*4)

write_tsv(Fig3A_raw$data,path = "Fig3A_draft_data_21072021.tsv")
write_tsv(Fig3B_raw$data,path = "Fig3B_draft_data_21072021.tsv")
write_tsv(Fig3Bg2_raw$data,path = "Fig3C_draft_data_21072021.tsv")
write_tsv(Fig3D_raw$data,path = "Fig3D_draft_data_21072021.tsv")


