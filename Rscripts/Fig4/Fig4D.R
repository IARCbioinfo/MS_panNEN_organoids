# load librairies
library(tidyverse)
library(ggpubr)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)

# theme 
theme_organoids <- theme_classic()+ grids()

colors_org = c(LNET2="#aade87ff",LNET6="#5fd38dff",LNET13="#16502dff",LNET14="#6f917cff",
               LNET5="#e6a73cff",LNET10="#ff9955ff",LNET15="#ffd42aff", LNET16 = "#ff6600ff", LNET18= "#d0742fff", 
               LNET19="#2aff80ff", 
               LNET20 = "#f6e62bff", 
               LCNEC3="#ff8080ff",LCNEC4="#d35f5fff", LCNEC23 = "#ff5555ff", 
               LCNEC11="#ff5599ff",PANEC1="#8d5fd3ff",
               SINET7="#2ad4ffff",SINET8="#80b3ffff",SINET9="#5f8dd3ff",SINET12="#5fbcd3ff", SINET21="#0066ffff", SINET22="#2c5aa0ff")

# metadata
require(readxl)
require(tidyverse)
mutated_gene_list1 = read_xlsx("/data/lungNENomics/work/organoids/metadata/Lung NEN and SI NEN mutated genes_16Oct2020.xlsx",sheet = 1)
mutated_gene_list1$Type = "LNEN"
mutated_gene_list2 = read_xlsx("/data/lungNENomics/work/organoids/metadata/Lung NEN and SI NEN mutated genes_16Oct2020.xlsx",sheet = 2)
mutated_gene_list2$Type = "SI"

mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="DICER"] = "DICER1"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="KMT2C/MLL3"] = "KMT2C"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="KMT2D/MLL2"] = "KMT2D"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="EML4-ALK"] = "ALK"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="NKX2-1/TTF1"] = "NKX2-1"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="GAS7"] = "AC005747.1"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="PIKSCA"] = "PIK3CA"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="ERCC9L"] = "ERCC6L"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="KAT9B"] = "ELP3"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="AR1D2"] = "ARID2"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="KRA3"] = "KRAS"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="TOP3"] = "TOP3A"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="HSP90"] = "HSP90A" # add others?
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="PDGFR"] = "PDGFRA" # add others
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="MLL4"] = "KMT2B"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="Hist1H2AC"] = "HIST1H2AC" 
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="DCTN"] = "DCTN1"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="FRMDP4"] = "FRMPD4"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="KIAA1598"] = "SHTN1"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="CNOT"] = "CNOT1"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="PIK3CA2"] = "PIK3C2A"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="TMEM49"] = "VMP1"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="MYH5"] = "MYH6"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="GTBP4"] = "GTPBP4"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="TCG7KL2"] = "TCF7L2"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="URB5"] = "UBR5"

mutated_gene_list  = bind_rows(mutated_gene_list1,mutated_gene_list2)
mutated_gene_list  = mutated_gene_list[!duplicated(mutated_gene_list$`gene name`),]
mutated_gene_list  = mutated_gene_list[!str_detect(mutated_gene_list$`gene name`,"hromosome"),]

mutated_gene_list = mutated_gene_list %>% mutate(Type=case_when((`gene name`%in% mutated_gene_list1$`gene name`) & (`gene name`%in% mutated_gene_list2$`gene name`) ~ "both",
                                                                `gene name`%in% mutated_gene_list1$`gene name` ~ "LNEN",
                                                                `gene name`%in% mutated_gene_list2$`gene name` ~ "SINET") )

mutated_gene_listB = mutated_gene_list %>% dplyr::select(`gene name`,Type) %>% mutate(WhichType=2*(Type=="both")+1*(Type=="LNEN")) %>% 
  filter(!`gene name`%in%c("ALK"))

mutated_gene_listB$`gene name`[mutated_gene_listB$`gene name`=="HSP90A"] = "HSP90AA1"
mutated_gene_listB$`gene name`[mutated_gene_listB$`gene name`=="HIST1H2AC"] = "H2AC6"

# Fig. 4D : CNV + SV profiles
CNVsummary = read_tsv("/data/lungNENomics/work/organoids/WGS/CNV_calling/release2_PURPLE_multisample_05072021/purple_summary.txt" )
CNVsummary = CNVsummary %>% mutate(experiment =str_remove(tumor_id,"[TMN][p0-9.]*$"))

## sup figure with purities and ploidies
FigS4C_A <- ggplot(CNVsummary , aes(x=ploidy,y=purity,col=experiment) ) + geom_point() + 
  geom_point(data = CNVsummary %>% filter(str_detect(tumor_id,"p")) , pch=16 , col="white",size=0.6 ) + scale_color_manual(values = colors_org[CNVsummary$experiment], drop = TRUE) + 
  geom_hline(yintercept = 0.5, linetype="dashed")  + theme_organoids + xlab("Ploidy") + ylab("Purity")
ggsave("/data/lungNENomics/work/organoids/figures/FigS4C_raw_26112021.svg",FigS4C_A, height = 3,width = 4)


CNVfiles = list.files("/data/lungNENomics/work/organoids/WGS/CNV_calling/release2_PURPLE_multisample_05072021/PURPLE/", recursive = T , 
                      pattern="purple.cnv.somatic.tsv",full.names = T)
CNVnames = gsub("[A-Z]+[0-9]+[NMT][p0-9.]*_PURPLE/|.purple.cnv.somatic.tsv","", list.files("/data/lungNENomics/work/organoids/WGS/CNV_calling/release2_PURPLE_multisample_05072021/PURPLE/", recursive = T , 
                                                         pattern="purple.cnv.somatic.tsv" ) )

CNVl = lapply(CNVfiles, read_tsv )
for(i in 1:length(CNVl)) {
  print(i)
  CNVl[[i]] = CNVl[[i]] %>% 
    mutate(Segment_Mean=log(copyNumber,2),Chromosome=chromosome,Sample=CNVnames[i] ) %>%
    filter( Chromosome %in% paste0("chr",c(1:23,"X","Y") ) , !(is.na(minorAlleleCopyNumber)|is.na(majorAlleleCopyNumber)) )
  CNVl[[i]]$Chromosome[CNVl[[i]]$Chromosome=="chr23"] = "chrX"
}

#seqinf = seqinfo(CNV_GR)
CNV_GR = lapply(CNVl, function(x) GRanges(seqnames = x$chromosome, ranges = IRanges(x$start,x$end), cn=x$majorAlleleCopyNumber + x$minorAlleleCopyNumber, 
                                          loh= x$minorAlleleCopyNumber==0, #clonal_fraction=x$cf.em, 
                                          minor_CN=x$minorAlleleCopyNumber) ) # ,seqinfo = seqinf) )
names(CNV_GR) = CNVnames

for(i in 1:length(CNVl)) names(CNV_GR[[i]]) = paste0(CNVl[[i]]$Chromosome,":",CNVl[[i]]$start,"-",CNVl[[i]]$end)

sampleOrderAll = c("LNET2T","LNET2Tp12","LNET2Np12",
  "LNET5T",	"LNET5Tp4",
  "LNET6T",	"LNET6Tp1",
  "LNET10T",	"LNET10Tp4",
  "LCNEC3T","LCNEC3Tp17",
  "LCNEC4T",	"LCNEC4Tp7", "LCNEC4Tp24",	
  "PANEC1T","PANEC1Tp4","PANEC1Tp14",
  "SINET7M",	"SINET7Mp2",
  "SINET8M",	"SINET8Mp2",
  "SINET9M",	"SINET9Mp1")


Data_Fig4D = bind_rows( lapply(1:length(sampleOrderAll),function(i){ 
  res=as.data.frame(CNV_GR[sampleOrderAll[i]])
  colnames(res)=c("seqnames","start","end","width","strand","cn","loh","minor_CN")
  res$Sample=sampleOrderAll[i];return(res)} ) )

write_tsv(Data_Fig4D,path = "Fig4D_draft_data_9122021.tsv")

# get gene info
gene_spans = read_tsv("/data/lungNENomics/work/organoids/references/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans",
                      col_names = c("ENSG","chr","start","end","strand","symbol","type") )
load("/data/lungNENomics/work/organoids/RNAseq/quantification/release2_all_03072021/Robjects/transcript_1pass.SE.rda")

## cross tx with CNV data
CNV_GR.S = CNV_GR
for(i in 1:length(CNV_GR.S)) CNV_GR.S[[i]] = CNV_GR[[i]][CNV_GR[[i]]$cn > 2.2 | CNV_GR[[i]]$cn < 1.8 | CNV_GR[[i]]$minor_CN < 0.8,]

mutated_gene_list_tx = rowRanges(transcript_1pass.SE)[which(rowData(transcript_1pass.SE)$gene_name %in% mutated_gene_listB$`gene name`)] # cdstmp.g[gene_spans[gene_spans$symbol %in% mutated_gene_listB$`gene name`,]$ENSG]
mutated_gene_list_tx.names = rowData(transcript_1pass.SE)$gene_name[which(rowData(transcript_1pass.SE)$gene_name %in% mutated_gene_listB$`gene name`)] # cdstmp.g[gene_spans[gene_spans$symbol %in% mutated_gene_listB$`gene name`,]$ENSG]


CNVmaf=c()
for(i in 1:length(CNVnames)){
  OV = findOverlaps(CNV_GR.S[CNVnames[i]][[1]],mutated_gene_list_tx)
  if(length(OV)>0){
    CNVmaf = rbind(CNVmaf, 
                   data.frame( gene = mutated_gene_list_tx.names[subjectHits(OV)] , 
                               Sample = CNVnames[i] , Variant_Classification = as_tibble( mcols(CNV_GR.S[CNVnames[i]][[1]][queryHits(OV),])) %>% 
                                 mutate(Variant_Classification=case_when( cn>2~"Amp", cn<2~"Del",loh~"LOH"  ) ) %>% pull(Variant_Classification) ) )
    CNVmaf = CNVmaf[!duplicated(CNVmaf),]
  }
}


##
bind_rows(CNVl) %>% dplyr::select( Sample, chromosome, start, end, bafCount , copyNumber )  %>% write_tsv(path=paste0("CNV_tsv_all.tsv"))

load("/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/SURVIVOR/SV_uniq_recovered.Rdata")


# create maf file
maf_subset = organoids_somatic.subset@data
vcf.SVl.filtered$type = sapply(vcf.SVl.filtered$msg, function(x) str_split(x," from| between")[[1]][1] )

maf_subset.wSV = rbind(maf_subset[,c(1:18,101:105)], data.table(Chromosome=vcf.SVl.filtered$seqnames,Start_Position=vcf.SVl.filtered$start,End_Position=vcf.SVl.filtered$altpos,
                                                                Reference_Allele=vcf.SVl.filtered$ref,Tumor_Seq_Allele2=vcf.SVl.filtered$alt,Tumor_Sample_Barcode=vcf.SVl.filtered$Sample,
                                                                Hugo_Symbol=vcf.SVl.filtered$geneSymbol,Variant_Classification="Structural_Variant",tx=NA,exon=vcf.SVl.filtered$exonGene,
                                                                txChange=NA,aaChange=NA,Variant_Type="SV",sample_id="multisample_annovar",Func.refGene="exonic",
                                                                Gene.refGene=vcf.SVl.filtered$geneSymbol,GeneDetail.refGene=".",ExonicFunc.refGene="Structural_Variant",
                                                                AD=vcf.SVl.filtered$AD,RD=vcf.SVl.filtered$DP,VAF=NA,Otherinfo16=NA,dam_pred=TRUE) )

write_tsv(maf_subset.wSV,"/data_alcalan/Organoids/organoids_coding_SV.maf")

#read CNVs
CNVs_org = bind_rows(CNVl) %>% select( Sample, Chromosome, Start, End, Num_Probes , Segment_Mean,tcn.em,lcn.em ) %>% 
  mutate(majorCN = tcn.em-lcn.em,minorCN = lcn.em) #read_tsv("CNV_tsv_all.tsv")

## plot
library(circlize)
library(RColorBrewer)

chr_lengths = tibble( chr=seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)@seqnames, 
                      length=seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)@seqlengths) %>% 
  filter( chr %in% c(paste0("chr",c(1:23,"X","Y")))  )

## useful metadata
colCNVs = c(rgb(0.8,0.2,0.1),rgb(1,1,1),rgb(0,0.4,0.7,0.25),rgb(0,0.4,0.7,0.5),rgb(0,0.4,0.7,0.75),rgb(0,0.4,0.7,1))
colCNVs = c( colorRampPalette(c(rgb(0.8,0.2,0.1), rgb(1,1,1)))(100) , colorRampPalette(c(rgb(1,1,1), rgb(0,0.4,0.7,1) ))(401) )

samplesTime = as.data.frame(matrix(c("LNET2T","LNET2Tp12",NA,
                                     "LNET5T",	"LNET5Tp4",NA,
                                     "LNET6T",	"LNET6Tp1",NA,
                                     "LNET10T",	"LNET10Tp4",NA,
                                     "LCNEC3T","LCNEC3Tp17",NA,
                                     "LCNEC4T",	"LCNEC4Tp7", "LCNEC4Tp24",	
                                     "PANEC1T","PANEC1Tp4","PANEC1Tp14",
                                     "SINET7M",	"SINET7Mp2",NA,
                                     "SINET8M",	"SINET8Mp2",NA,
                                     "SINET9M",	"SINET9Mp1",NA) , ncol=3,byrow=T,dimnames = list(c("LNET2" , "LNET5", "LNET6", "LNET10","LCNEC3","LCNEC4","PANEC1","SINET7","SINET8","SINET9"),
                                                                                                 c("Parental","PDTO1","PDTO2"))) )

## check if all driver genes in annotation
all(mutated_gene_listB$`gene name` %in% ref.33.cds$geneSymbol)

hg38_bed = read.table("/data/lungNENomics/work/organoids/references/hs38DH/hs38DH_main.bed",row.names = 1,header = F)
rownames(hg38_bed) = c("start","end")
## plot function
plot_crc1 <- function(genes_to_add=driver_uniq_recovered, experiment){
  vcfs = SV_uniq_recovered[SV_uniq_recovered$Experiment==experiment,]
  
  par(mfrow=c(1,3),mar = c(3, 3, 3, 3)*0.5)
  for(i in 1:sum(!is.na(samplesTime[experiment,]))){  
    sample.tmp = samplesTime[experiment,i]
    circos.clear()
    col_text <- rgb(0.2,0.2,0.2)#"grey40"
    circos.par("track.height" = 0.20, start.degree = 90, circle.margin = 0.35 , #canvas.xlim = c(-1, 1), canvas.ylim = c(-1, 1)
               gap.degree = 3,
               cell.padding = c(0, 0, 0, 0))
    circos.initialize(sectors = rownames(hg38_bed),xlim = hg38_bed )
    
    # add gene labels
    if(any(!(unlist((genes_to_add[as.data.frame(genes_to_add)[,12+i]!="NO",] %>% filter(Experiment==experiment))$Gene) %in% ref.33.cds$geneSymbol ) ) ) warning("Some driver genes not in annotation")
    genes_to_add.tmp = ref.33.cds[ref.33.cds$geneSymbol %in% unlist((genes_to_add[as.data.frame(genes_to_add)[,12+i]!="NO",] %>% filter(Experiment==experiment))$Gene),]
    genes_to_add.tmp = genes_to_add.tmp[!duplicated(genes_to_add.tmp$geneSymbol),] %>% mutate(chr=paste0("chr",seqnames),gene=as.character(geneSymbol))
    if(nrow(genes_to_add.tmp)>0 ){
      circos.genomicLabels( genes_to_add.tmp %>% dplyr::select(chr,start,end,gene) %>% as.data.frame(),
                            labels.column=4, cex=0.8, col="black",line_lwd=1,line_col="black", #"grey80",
                            side="outside",connection_height=0.07,labels_height=0.04)
    }
    
    # genomes
    circos.track(ylim=c(0,1),panel.fun=function(x,y) {
      chr=str_remove(CELL_META$sector.index,"chr")
      xlim=CELL_META$xlim
      ylim=CELL_META$ylim
      circos.text(mean(xlim),mean(ylim),chr,cex=0.8,col=col_text,
                  facing="bending.inside",niceFacing=TRUE)
    },bg.col=rep(c(rgb(0.85,0.85,0.85),rgb(0.7,0.7,0.7)),12),#c(brewer.pal(12,"Paired"),brewer.pal(12,"Set3")),
    bg.border=F,track.height=0.06)#,track.index=3)
    
    # add copy number
    CNVs_org.tmp = CNVl[[which(CNVnames==sample.tmp)]] %>% #%>% filter(majorCN!=1 | minorCN!=1 ) %>% 
      mutate(chr = chromosome,value=1,value1=majorAlleleCopyNumber,value2=minorAlleleCopyNumber)
    circos.genomicTrack( CNVs_org.tmp %>% dplyr::select(chr,start,end,value,value1,value2) , ylim=c(0,4), bg.border=NA, #track.index=4 , 
                         panel.fun = function(region, value, ...) {
                           region1 = region[value$value1!=1,]
                           region2 = region[value$value2!=1,]
                           if(nrow(region1)>0) circos.genomicRect(region1,value[value$value1!=1,1],
                                                                  col=colCNVs[round(value$value1[value$value1!=1]*100+1)],
                                                                  border = NA, ytop = 2,ybottom = 0 )
                           if(nrow(region2)>0) circos.genomicRect(region2,value[value$value2!=1,1],
                                                                  col=colCNVs[round(value$value2[value$value2!=1]*100+1)],
                                                                  border = NA, ytop = 4,ybottom = 2 )
                         })
    
    
    # add rearrangements
    Link1 = vcfs[as.data.frame(vcfs)[,12+i]!="NO",] %>% mutate(chrom=CHROM,start=start.B1,end=end.B1) %>% dplyr::select(chrom,start,end) %>% as.data.frame()
    Link2 = vcfs[as.data.frame(vcfs)[,12+i]!="NO",] %>% mutate(chrom=CHROM2,start=start.B2,end=end.B2) %>% dplyr::select(chrom,start,end) %>% as.data.frame()
    circos.genomicLink(Link1,Link2,col = c(rgb(0.8,0.4,0.4),rgb(0,0.7,0.8))[as.numeric(vcfs[as.data.frame(vcfs)[,12+i]!="NO",]$type=="intra")+1]  ) #,col=rcols,border=NA)
    title( sample.tmp )
  }
}

for(exp.tmp in rownames(samplesTime)){
  pdf(paste0("/data/lungNENomics/work/organoids/figures/Fig4D_raw_15122021_",exp.tmp,".pdf"),height = 4*1,width = 4*2 )
  plot_crc1(experiment=exp.tmp)
  dev.off()
}

### Table 2 SV proportions
SVburdens = SV_uniq_recovered %>% group_by(Experiment) %>% summarize(tot = , sumP = sum(Parental!="NO"), sumPDTO1 = sum(PDTO1!="NO") , sumPDTO2 = sum(PDTO2!="NO",na.rm=T) )
SV_Venn = SV_uniq_recovered %>% group_by(Experiment) %>% summarize(Tot = n(), 
                                                         #Common     = sum(Parental!="NO" & PDTO1!="NO" & (PDTO2!="NO"| is.na(PDTO2)) ), 
                                                         "P&O1&O2"  = sum(Parental!="NO" & PDTO1!="NO" & (!is.na(PDTO2) & (PDTO2!="NO") )),
                                                         "P&O1"  = sum(Parental!="NO" & PDTO1!="NO" & (is.na(PDTO2) | (PDTO2=="NO") )),
                                                         "P&O2"  = sum(Parental!="NO" & PDTO1=="NO" & ((!is.na(PDTO2)) & (PDTO2!="NO") )),
                                                         P     = sum(Parental!="NO" & PDTO1=="NO" & (is.na(PDTO2) | (PDTO2=="NO") )),
                                                         O1 = sum(Parental=="NO" & PDTO1!="NO" & (is.na(PDTO2) | (PDTO2=="NO") )),
                                                         O2 = sum(Parental=="NO" & PDTO1=="NO" ) )

write_tsv( SV_Venn %>% mutate(Prop_P_conserved = 1-P/(`P&O1&O2`+ `P&O1`+ `P&O2`+P) ) %>% dplyr::select(Experiment,Prop_P_conserved) , 
           "/data/lungNENomics/work/organoids/figures/Table2_SVprop.tsv")

library(eulerr)
error_plot()

# Fig. 4D
# Fig4D
require(readxl)
require(tidyverse)
mutated_gene_list1 = read_xlsx("Lung NEN and SI NEN mutated genes_16Oct2020.xlsx",sheet = 1)
mutated_gene_list1$Type = "LNEN"
mutated_gene_list2 = read_xlsx("Lung NEN and SI NEN mutated genes_16Oct2020.xlsx",sheet = 2)
mutated_gene_list2$Type = "SI"

mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="DICER"] = "DICER1"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="KMT2C/MLL3"] = "KMT2C"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="KMT2D/MLL2"] = "KMT2D"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="EML4-ALK"] = "ALK"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="NKX2-1/TTF1"] = "NKX2-1"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="GAS7"] = "AC005747.1"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="PIKSCA"] = "PIK3CA"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="ERCC9L"] = "ERCC6L"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="KAT9B"] = "ELP3"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="AR1D2"] = "ARID2"
mutated_gene_list1$`gene name`[mutated_gene_list1$`gene name`=="KRA3"] = "KRAS"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="TOP3"] = "TOP3A"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="HSP90"] = "HSP90A" # add others?
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="PDGFR"] = "PDGFRA" # add others
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="MLL4"] = "KMT2B"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="Hist1H2AC"] = "HIST1H2AC" 
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="DCTN"] = "DCTN1"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="FRMDP4"] = "FRMPD4"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="KIAA1598"] = "SHTN1"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="CNOT"] = "CNOT1"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="PIK3CA2"] = "PIK3C2A"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="TMEM49"] = "VMP1"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="MYH5"] = "MYH6"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="GTBP4"] = "GTPBP4"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="TCG7KL2"] = "TCF7L2"
mutated_gene_list2$`gene name`[mutated_gene_list2$`gene name`=="URB5"] = "UBR5"

mutated_gene_list  = bind_rows(mutated_gene_list1,mutated_gene_list2)
mutated_gene_list  = mutated_gene_list[!duplicated(mutated_gene_list$`gene name`),]
mutated_gene_list  = mutated_gene_list[!str_detect(mutated_gene_list$`gene name`,"hromosome"),]

mutated_gene_list = mutated_gene_list %>% mutate(Type=case_when((`gene name`%in% mutated_gene_list1$`gene name`) & (`gene name`%in% mutated_gene_list2$`gene name`) ~ "both",
                                                                `gene name`%in% mutated_gene_list1$`gene name` ~ "LNEN",
                                                                `gene name`%in% mutated_gene_list2$`gene name` ~ "SINET") )

mutated_gene_listB = mutated_gene_list %>% select(`gene name`,Type) %>% mutate(WhichType=2*(Type=="both")+1*(Type=="LNEN")) %>% 
  filter(!`gene name`%in%c("ALK"))


# Fig 4D
sampleOrder = c(#"LNET2T","LNET2Tp12",#"LNET2Np12",
  #"LNET5T",	"LNET5Tp4",
  "LNET6T",	"LNET6Tp1",
  "LNET10T",	"LNET10Tp4",
  "LCNEC3T","LCNEC3Tp17",
  "LCNEC4T",	"LCNEC4Tp7", "LCNEC4Tp24",	
  "PANEC1T","PANEC1Tp4","PANEC1Tp14",
  "SINET7M",	"SINET7Mp2",
  "SINET8M",	"SINET8Mp2",
  "SINET9M",	"SINET9Mp1")

organoids_somatic.subset = subsetMaf(organoids_somatic,tsb = sampleOrder,query = "dam_pred==TRUE" )#,query = "dam_pred==TRUE")

# add SVs 
maf_subset.wSV = read.maf("/data_alcalan/Organoids/organoids_coding_SV.maf",
                          vc_nonSyn = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", 
                                        "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","Structural_Variant"))
# 

maf_subset.wSV@data = maf_subset.wSV@data[!duplicated(cbind(maf_subset.wSV@data$Hugo_Symbol,maf_subset.wSV@data$Variant_Classification,maf_subset.wSV@data$Tumor_Sample_Barcode)),]

tabgenes = maf_subset.wSV@data %>% filter(Hugo_Symbol%in%mutated_gene_listB$`gene name` ) %>% 
  mutate(Exp=str_remove(Tumor_Sample_Barcode,"[MT]p*[0-9]*$")) %>% 
  group_by(Hugo_Symbol) %>% summarise(n=n(),nexp=length(unique(Exp)),Exp ) %>% left_join(mutated_gene_listB,by=c("Hugo_Symbol"="gene name")) %>% 
  mutate(WhichType=(WhichType==2)*(-1.5)+WhichType,Type1=WhichType==1,Type2=WhichType==2 ) %>% 
  mutate(Exp=factor(Exp,levels =  levels_nmut)) %>% 
  arrange(Exp,desc(nexp),desc(n))
levels_normal = c("SINET7","SINET8","SINET9","LNET6","LNET10","LCNEC3","LCNEC4","PANEC1")
levels_nmut = c("LCNEC4", "PANEC1", "LCNEC3", "LNET10", "SINET8", "SINET7", "SINET9", "LNET6")
#tabgenes$WhichType[tabgenes$WhichType==2] = 0.5

sampleOrder_nmut = sampleOrder[c(7:12,5:6,3:4,15:16,13:14,17:18,1:2)]
sampleOrder_onco = c(sampleOrder[13:18],sampleOrder[1:12])
vc_col = maftools:::get_vcColors(websafe = FALSE)
vc_col["Structural_Variant"] = "#a02c89ff"
vc_col["Missense_Mutation"] = "#1f78b4ff"
vc_col["Frame_Shift_Del"] = "#782121ff"
vc_col["Frame_Shift_Ins"] = "#2ca089ff"

pdf("Fig4D_raw_30042021c.pdf",h=8.5,w=5.5)
oncoplot(maf = maf_subset.wSV, showTumorSampleBarcodes = T, barcode_mar = 6, drawRowBar = F,drawColBar = F,  showTitle = F,
         #additionalFeature = c("dam_pred075",TRUE), 
         leftBarData = tabgenes[!duplicated(tabgenes$Hugo_Symbol),c(1,6)], 
         sampleOrder = sampleOrder_nmut,#sampleOrder_onco,
         genes = unique(tabgenes$Hugo_Symbol), 
         keepGeneOrder = TRUE,
         removeNonMutated = FALSE,
         colors = vc_col,
         borderCol= "white",
         bgCol = "#fef8e4ff")
dev.off()

# write data
write_tsv( maf_subset.wSV@data ,path = "Fig4D_draft_data_16042021.tsv" )
