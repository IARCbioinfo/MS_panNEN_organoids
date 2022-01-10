# load libraries
library(maftools)
library(VariantAnnotation)
library(tidyverse)

# create sample and passage data frames
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

samplesTimeRNA = as.data.frame(matrix(c("LNET2T","LNET2Tp12",NA,NA,
                                     "LNET5T",	"LNET5Tp4","LNET5Tp7","LNET5Tp2.2",
                                     "LNET6T",	"LNET6Tp1",NA,NA,
                                     "LNET10T",	"LNET10Tp4","LNET10Tp11",NA,
                                     "LNET13T",	"LNET13Tp1",NA,NA,
                                     "LNET14T",	"LNET14Tp1",NA,NA,
                                     "LNET15M",	"LNET15Mp2",NA,NA,
                                     "LNET16T",	"LNET16Tp2","LNET16M",	"LNET16Mp1",
                                     "LNET18T",	"LNET18Tp2",NA,NA,
                                     "LNET19T",	"LNET19Tp2",NA,NA,
                                     "LNET20M",	"LNET20Mp2",NA,NA,
                                     "LCNEC3T","LCNEC3Tp17.2","LCNEC3Tp24",NA,
                                     "LCNEC4T",	"LCNEC4Tp7", "LCNEC4Tp24",	NA,
                                     "LCNEC11M","LCNEC11Mp3",NA,NA,
                                     NA,"LCNEC23Mp3",NA,NA,
                                     "PANEC1T","PANEC1Tp4","PANEC1Tp14",NA,
                                     "SINET7M",	"SINET7Mp2",NA,NA,
                                     "SINET8M",	"SINET8Mp2",NA,NA,
                                     "SINET12M",	"SINET12Mp1.1","SINET12Mp1.3",NA,
                                     "SINET21M",	"SINET21Mp2",NA,NA,
                                     "SINET22M",	"SINET22Mp2",NA,NA
                                     ) , ncol=4,byrow=T,dimnames = list(c("LNET2" , "LNET5", "LNET6", "LNET10","LNET13","LNET14","LNET15","LNET16","LNET18","LNET19","LNET20",
                                                                          "LCNEC3","LCNEC4","LCNEC11","LCNEC23","PANEC1","SINET7","SINET8","SINET12","SINET21","SINET22"),
                                                                                                 c("Parental","PDTO1","PDTO2","PDTO3"))) )

sampleOrder = c("LNET2T","LNET2Tp12",
  "LNET6T",	"LNET6Tp1",
  "LNET5T",	"LNET5Tp4",
  "LNET10T",	"LNET10Tp4",
  "LCNEC3T","LCNEC3Tp17",
  "LCNEC4T",	"LCNEC4Tp7", "LCNEC4Tp24",	
  "PANEC1T","PANEC1Tp4","PANEC1Tp14",
  "SINET7M",	"SINET7Mp2",
  "SINET8M",	"SINET8Mp2",
  "SINET9M",	"SINET9Mp1")

sampleOrderRNA = c("SINET21M",	"SINET21Mp2",
                   "SINET22M",	"SINET22Mp2",
                   "SINET7M",	"SINET7Mp2",
                   "SINET8M",	"SINET8Mp2",
                   "SINET12M",	"SINET12Mp1.1","SINET12Mp1.3",
                   "LNET6T",	"LNET6Tp1",
                   "LNET13T",	"LNET13Tp1",
                   "LNET14T",	"LNET14Tp1",
                   "LNET19T",	"LNET19Tp2",
                "LNET5T",	"LNET5Tp4","LNET5Tp7","LNET5Tp2.2",
                "LNET10T",	"LNET10Tp4","LNET10Tp11",
                "LNET15M",	"LNET15Mp2",
                "LNET16T",	"LNET16Tp2","LNET16M",	"LNET16Mp1",
                "LNET18Tp2",
                "LNET20M",	"LNET20Mp2",
                "LCNEC3T","LCNEC3Tp17.2","LCNEC3Tp24",
                "LCNEC4T",	"LCNEC4Tp7", "LCNEC4Tp24",	
                "LCNEC11M","LCNEC11Mp3",
                "LCNEC23Mp3",
                "PANEC1T","PANEC1Tp4","PANEC1Tp14"
                )

## convert annovar to maf
annovar_org = list.files("/data/lungNENomics/work/organoids/WGS/variant_calling/release3_25102021/",pattern=".vcf.gz$",full.names = T)
annovar_org.names = str_remove( list.files("/data/lungNENomics/work/organoids/WGS/variant_calling/release3_25102021/",pattern=".vcf.gz$"), 
                                "_final.vcf.gz")

# create MAF file
vcf.tib.long = c()
alt_summary = c()
for(i in 1:length(annovar_org)){
vcff <- VcfFile(annovar_org[i])
open(vcff)
param <- ScanVcfParam(info=c("CONTQ","DP","Func.ensGene","Gene.ensGene","GeneDetail.ensGene","ExonicFunc.ensGene","AAChange.ensGene",
                                                 "REVEL","cosmic92_coding","cosmic92_noncoding"))
vcf = readVcf(vcff,"hg38",param = param)

vcf.tib.tmp = tibble( Chr=as.character(seqnames(rowRanges(vcf))) , Start=start(rowRanges(vcf)), End=end(rowRanges(vcf)), Ref=as.character(rowRanges(vcf)$REF) , 
        Alt=as.character(unlist(rowRanges(vcf)$ALT)) , Func.ensGene=as.character(info(vcf)$Func.ensGene) , 
        Gene.ensGene = str_replace(as.character(info(vcf)$Gene.ensGene),"\\\\x3b",";") , 
        GeneDetail.ensGene = str_replace_all(str_replace_all(as.character(info(vcf)$GeneDetail.ensGene),"\\\\x3d",":"),"\\\\x3b",";")  , 
        ExonicFunc.ensGene = as.character(info(vcf)$ExonicFunc.ensGene),    
        AAChange.ensGene = sapply( info(vcf)$AAChange.ensGene, function(x) paste(x,collapse = ";")) , 
        REVEL = as.numeric(as.character(info(vcf)$REVEL)) )

# find columns
ID.P = which(colnames(geno(vcf)$AD)==samplesTime[annovar_org.names[i],]$Parental)
ID.O1 = which(colnames(geno(vcf)$AD)==samplesTime[annovar_org.names[i],]$PDTO1)
ID.O2 = which(colnames(geno(vcf)$AD)==samplesTime[annovar_org.names[i],]$PDTO2)
vcf.tib.tmp$DP.P = geno(vcf)$DP[,ID.P]
vcf.tib.tmp$DP.O1 = geno(vcf)$DP[,ID.O1]
vcf.tib.tmp$AD.P = sapply( geno(vcf)$AD[,ID.P], function(x) x[2])
vcf.tib.tmp$AD.O1 = sapply( geno(vcf)$AD[,ID.O1], function(x) x[2])
if(length(ID.O2)>0){
  vcf.tib.tmp$AD.O2 = sapply( geno(vcf)$AD[,ID.O2], function(x) x[2])
  vcf.tib.tmp$DP.O2 = geno(vcf)$DP[,ID.O2]
}else{
  vcf.tib.tmp$AD.O2 = NA
  vcf.tib.tmp$DP.O2 = NA
}
vcf.tib.tmp$VAF.P  = vcf.tib.tmp$AD.P/vcf.tib.tmp$DP.P
vcf.tib.tmp$VAF.O1 = vcf.tib.tmp$AD.O1/vcf.tib.tmp$DP.O1
vcf.tib.tmp$VAF.O2 = vcf.tib.tmp$AD.O2/vcf.tib.tmp$DP.O2

# find number of shared and private with DP >=30
alt_summary = bind_rows(alt_summary, vcf.tib.tmp %>% filter(DP.P>=50 & DP.O1>=50 & (is.na(DP.O2) | DP.O2>=50)  ) %>% summarize( "P&O1&O2"= sum(AD.P>0 & AD.O1>0 & !is.na(AD.O2) & AD.O2>0) , 
                                                                                           "P&O1" = sum(AD.P>0 & AD.O1>0 & (is.na(AD.O2) | AD.O2==0)) , 
                                                                                           "P&O2" = sum(AD.P>0 & AD.O1==0 & !is.na(AD.O2) & AD.O2>0) ,
                                                                                           "O1&O2"= sum(AD.P==0 & AD.O1>0 & !is.na(AD.O2) & AD.O2>0) , 
                                                                                           P = sum(AD.P>0 & AD.O1==0 & (is.na(AD.O2) | AD.O2==0)) ,
                                                                                           O1 = sum(AD.P==0 & AD.O1>0 & (is.na(AD.O2) | AD.O2==0)) , 
                                                                                           O2 = sum(AD.P==0 & AD.O1==0)
                                                                                           ) %>% mutate(Experiment=annovar_org.names[i]) )

vcf.tib.tmp.long = vcf.tib.tmp %>% pivot_longer(cols = AD.P:AD.O2,names_to = c("geno","Sample"), names_sep = "\\.",values_to = "AD") %>% filter(!is.na(AD) & AD>0)
vcf.tib.tmp.long$DP = vcf.tib.tmp.long$DP.O2
vcf.tib.tmp.long$DP[vcf.tib.tmp.long$Sample=="P"] = vcf.tib.tmp.long$DP.P[vcf.tib.tmp.long$Sample=="P"]
vcf.tib.tmp.long$DP[vcf.tib.tmp.long$Sample=="O1"] = vcf.tib.tmp.long$DP.O1[vcf.tib.tmp.long$Sample=="O1"]

vcf.tib.tmp.long$VAF = vcf.tib.tmp.long$VAF.O2
vcf.tib.tmp.long$VAF[vcf.tib.tmp.long$Sample=="P"] = vcf.tib.tmp.long$VAF.P[vcf.tib.tmp.long$Sample=="P"]
vcf.tib.tmp.long$VAF[vcf.tib.tmp.long$Sample=="O1"] = vcf.tib.tmp.long$VAF.O1[vcf.tib.tmp.long$Sample=="O1"]
vcf.tib.tmp.long = vcf.tib.tmp.long %>% dplyr::select(-(DP.P:VAF.O2),-geno ) %>% mutate(Sample = case_when(Sample=="P"~samplesTime[annovar_org.names[i],]$Parental,
                                               Sample=="O1"~samplesTime[annovar_org.names[i],]$PDTO1,
                                               TRUE~samplesTime[annovar_org.names[i],]$PDTO2) )
vcf.tib.long = bind_rows(vcf.tib.long , vcf.tib.tmp.long)
}

# Switch to recognized annovar values
vcf.tib.long$ExonicFunc.ensGene = str_replace(vcf.tib.long$ExonicFunc.ensGene, "_", " ")

annovar_values = c(exonic = "RNA", splicing = "Splice_Site", 
                   UTR5 = "5'UTR", UTR3 = "3'UTR", intronic = "Intron", 
                   upstream = "5'Flank", downstream = "3'Flank", intergenic = "IGR", 
                   `frameshift insertion` = "Frame_Shift_Ins", `frameshift deletion` = "Frame_Shift_Del", 
                   `frameshift block substitution` = "Frameshift_INDEL", 
                   `frameshift substitution` = "Frameshift_INDEL", stopgain = "Nonsense_Mutation", 
                   stoploss = "Nonstop_Mutation", startloss = "Translation_Start_Site", 
                   startgain = "Unknown", `nonframeshift insertion` = "In_Frame_Ins", 
                   `nonframeshift deletion` = "In_Frame_Del", `nonframeshift block substitution` = "Inframe_INDEL", 
                   `nonframeshift substitution` = "Inframe_INDEL", `nonsynonymous SNV` = "Missense_Mutation", 
                   `synonymous SNV` = "Silent", unknown = "Unknown", ncRNA_exonic = "RNA", 
                   ncRNA_intronic = "RNA", ncRNA_UTR3 = "RNA", ncRNA_UTR5 = "RNA", 
                   ncRNA = "RNA", ncRNA_splicing = "RNA")

all(vcf.tib.long$ExonicFunc.ensGene[vcf.tib.long$ExonicFunc.ensGene!="."] %in% names(annovar_values))

write_tsv(alt_summary , file = "/data/lungNENomics/work/organoids/WGS/variant_calling/release3_25102021/mutations_venn.tsv")
write_tsv(vcf.tib.long , file = "/data/lungNENomics/work/organoids/WGS/variant_calling/release3_25102021/multisample_annovar.maf")

cbind(alt_summary$Experiment , 1 - alt_summary$P / (rowSums(alt_summary[,c(1:3)]) + alt_summary$P ))

# Figure 4B
# load MAF file
organoids_somatic = annovarToMaf(annovar = "/data/lungNENomics/work/organoids/WGS/variant_calling/release3_25102021/multisample_annovar.maf", 
                                 Center = 'Utrecht', refBuild = 'hg38', tsbCol = 'Sample', table = 'ensGene',MAFobj = TRUE,ens2hugo = FALSE,)

## test if all in exonic
table(organoids_somatic@data$ExonicFunc.refGene, organoids_somatic@data$Variant_Classification )
organoids_somatic@data %>% filter( ExonicFunc.refGene=="." )
table(organoids_somatic@maf.silent$ExonicFunc.refGene, organoids_somatic@maf.silent$Variant_Classification)
organoids_somatic@maf.silent %>% filter( !ExonicFunc.refGene %in% c(".","synonymous SNV") )
### OK

## add SVs
maf.data = organoids_somatic@data

SV.maf = driver_uniq_recovered %>% filter(ID!="") %>% pivot_longer(cols = Parental:PDTO2,names_to = "Sample") %>% filter(!is.na(value) & value!="NO" )
SV.maf$Sample = apply(SV.maf ,1 , function(x) samplesTime[x$Experiment,x$Sample])

SV.maf = SV.maf %>% mutate(Chromosome= CHROM, Start_Position=start.B1 , End_Position=end.B1 , Reference_Allele = NA , Tumor_Seq_Allele2 = NA , Tumor_Sample_Barcode = Sample , 
                           Hugo_Symbol = Gene, Variant_Classification =  SVTYPE, tx=NA , exon=NA , txChange=NA , 
                           aaChange=NA , Variant_Type = "SV" , sample_id = "multisample_SV", Func.refGene ="exonic", Gene.refGene=Gene , 
                           GeneDetail.refGene =ID , ExonicFunc.refGene = SVTYPE , AAChange.refGene =NA , REVEL =NA , AD =NA,
                           DP=NA , VAF=NA) %>% dplyr::select(Chromosome:VAF)
SV.maf$Hugo_Symbol = sapply(SV.maf$Gene.refGene,function(x){ res=NA; if(length(x)>0){res=as.character(x)} ; return(res) })
SV.maf$Gene.refGene = SV.maf$Hugo_Symbol

SV.maf = SV.maf[which(!duplicated(SV.maf[,c(1,6:8,13:14)])),]
maf.data = rbind(maf.data,SV.maf)

# check if damaging
maf.data = maf.data %>% mutate(damaging = case_when(as.numeric(REVEL)>=0.5 | Variant_Classification !="Missense_Mutation" ~"YES",
                                         as.numeric(REVEL)<0.5 | Variant_Classification !="Missense_Mutation" ~"MAYBE",
                                         TRUE~"NO") )

write_tsv( maf.data ,file = "/data/lungNENomics/work/organoids/WGS/variant_calling/release3_25102021/multisample_annovar_wSV.maf"  )

organoids_somatic2 = read.maf("/data/lungNENomics/work/organoids/WGS/variant_calling/release3_25102021/multisample_annovar_wSV.maf" ,
                              vc_nonSyn = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", 
                                            "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","INV","TRA","DEL","DUP"))
                                 
## plot 
#sampleOrder_nmut = sampleOrder[c(7:12,5:6,3:4,15:16,13:14,17:18,1:2)]
sampleOrder_onco = c(sampleOrder[17:22],sampleOrder[1:16])
levels_normal = c("SINET7","SINET8","SINET9","LNET2","LNET6","LNET5","LNET10","LCNEC3","LCNEC4","PANEC1")
#levels_nmut = c("LCNEC4", "PANEC1", "LCNEC3", "LNET10", "SINET8", "SINET7", "SINET9", "LNET6")

organoids_somatic2.dam = subsetMaf(organoids_somatic2,query = "damaging=='YES'")
tabgenes = organoids_somatic2.dam@data %>% filter(Hugo_Symbol%in%mutated_gene_listB$`gene name` ) %>% 
  mutate(Exp=str_remove(Tumor_Sample_Barcode,"[MT]p*[0-9]*$")) %>% 
  group_by(Hugo_Symbol) %>% summarise(n=n(),nexp=length(unique(Exp)),Exp ) %>% left_join(mutated_gene_listB,by=c("Hugo_Symbol"="gene name")) %>% 
  mutate(WhichType=(WhichType==2)*(-1.5)+WhichType,Type1=WhichType==1,Type2=WhichType==2 ) %>% 
  mutate(Exp=factor(Exp,levels =  levels_normal)) %>% 
  arrange(Exp,desc(nexp),desc(n))
#tabgenes$Exp = 
#tabgenes$WhichType[tabgenes$WhichType==2] = 0.5

vc_col = maftools:::get_vcColors(websafe = FALSE)
vc_col["INV"] = "#a02c89ff"
vc_col["TRA"] = "#42b6f5ff"
vc_col["DEL"] = "#d647b9ff"
vc_col["DUP"] = "#5c47d6ff"
vc_col["Missense_Mutation"] = "#1f78b4ff"
vc_col["Frame_Shift_Del"] = "#782121ff"
vc_col["Frame_Shift_Ins"] = "#2ca089ff"

pdf("/data/lungNENomics/work/organoids/figures/Fig4B_raw_17122021.pdf",h=10.5,w=6.5)
oncoplot(maf = organoids_somatic2.dam, showTumorSampleBarcodes = T, barcode_mar = 6, drawRowBar = F,drawColBar = F,  showTitle = F,
         #additionalFeature = c("dam_pred075",TRUE), 
         leftBarData = tabgenes[!duplicated(tabgenes$Hugo_Symbol),c(1,6)], 
         sampleOrder = sampleOrder_onco,#sampleOrder_onco,
         genes = unique(tabgenes$Hugo_Symbol), 
         keepGeneOrder = TRUE,
         removeNonMutated = FALSE,
         colors = vc_col,
         borderCol= "white",
         bgCol = "#fef8e4ff")
dev.off()

# write data
write_tsv( organoids_somatic2@data ,path = "Fig4B_draft_data_16042021.tsv" )

## C: RNAseq oncoplot panel
annovar_org_RNA = list.files("/data/lungNENomics/work/organoids/integration/Tonly_RNAseq_vcf_drivers_normalized_annotated/",pattern=".vcf.gz$",full.names = T)
annovar_org_RNA.names = str_remove( list.files("/data/lungNENomics/work/organoids/integration/Tonly_RNAseq_vcf_drivers_normalized_annotated/",pattern=".vcf.gz$"), 
                                "_calls_driverRegions_norm.vcf.hg38_multianno.vcf.gz")

# create MAF file
vcfRNA.tib.long = c()
alt_summary_RNA = c()
for(i in 1:length(annovar_org_RNA)){
  vcff <- VcfFile(annovar_org_RNA[i])
  open(vcff)
  param <- ScanVcfParam(info=c("CONTQ","DP","Func.ensGene","Gene.ensGene","GeneDetail.ensGene","ExonicFunc.ensGene","AAChange.ensGene",
                               "REVEL","cosmic92_coding","cosmic92_noncoding","non_cancer_AF_popmax","ExAC_nontcga_ALL","1000g2015aug_all"))
  
  vcf = readVcf(vcff,"hg38",param = param)
  
  vcfRNA.tib.tmp = tibble( Chr=as.character(seqnames(rowRanges(vcf))) , Start=start(rowRanges(vcf)), End=end(rowRanges(vcf)), Ref=as.character(rowRanges(vcf)$REF) , 
                        Alt=as.character(unlist(rowRanges(vcf)$ALT)) , Func.ensGene=as.character(info(vcf)$Func.ensGene) , 
                        Gene.ensGene = str_replace(as.character(info(vcf)$Gene.ensGene),"\\\\x3b",";") , 
                        GeneDetail.ensGene = str_replace_all(str_replace_all(as.character(info(vcf)$GeneDetail.ensGene),"\\\\x3d",":"),"\\\\x3b",";")  , 
                        ExonicFunc.ensGene = as.character(info(vcf)$ExonicFunc.ensGene),    
                        AAChange.ensGene = sapply( info(vcf)$AAChange.ensGene, function(x) paste(x,collapse = ";")) , 
                        REVEL = as.numeric(as.character(info(vcf)$REVEL)),cosmic92_coding = sapply(info(vcf)$cosmic92_coding , paste, collapse=";"),
                        cosmic92_noncoding = sapply(info(vcf)$cosmic92_noncoding , paste, collapse=";"), 
                        non_cancer_AF_popmax = as.numeric(as.character(info(vcf)$non_cancer_AF_popmax)),
                        ExAC_nontcga_ALL = as.numeric(as.character(info(vcf)$ExAC_nontcga_ALL)),
                        `1000g2015aug_all` = as.numeric(as.character(info(vcf)$`1000g2015aug_all`))
                        )
  
  # find columns
  ID.P = which(colnames(geno(vcf)$AD)==samplesTimeRNA[annovar_org_RNA.names[i],]$Parental)
  ID.O1 = which(colnames(geno(vcf)$AD)==samplesTimeRNA[annovar_org_RNA.names[i],]$PDTO1)
  ID.O2 = which(colnames(geno(vcf)$AD)==samplesTimeRNA[annovar_org_RNA.names[i],]$PDTO2)
  ID.O3 = which(colnames(geno(vcf)$AD)==samplesTimeRNA[annovar_org_RNA.names[i],]$PDTO3)
  
  vcfRNA.tib.tmp$DP.O1 = geno(vcf)$DP[,ID.O1]
  if(length(ID.P)>0){
    vcfRNA.tib.tmp$DP.P = geno(vcf)$DP[,ID.P]
    vcfRNA.tib.tmp$AD.P = sapply( geno(vcf)$AD[,ID.P], function(x) x[2])
  }else{
    vcfRNA.tib.tmp$DP.P = NA
    vcfRNA.tib.tmp$AD.P = NA
  }
  vcfRNA.tib.tmp$AD.O1 = sapply( geno(vcf)$AD[,ID.O1], function(x) x[2])
  if(length(ID.O2)>0){
    vcfRNA.tib.tmp$AD.O2 = sapply( geno(vcf)$AD[,ID.O2], function(x) x[2])
    vcfRNA.tib.tmp$DP.O2 = geno(vcf)$DP[,ID.O2]
  }else{
    vcfRNA.tib.tmp$AD.O2 = NA
    vcfRNA.tib.tmp$DP.O2 = NA
  }
  if(length(ID.O3)>0){
    vcfRNA.tib.tmp$AD.O3 = sapply( geno(vcf)$AD[,ID.O3], function(x) x[2])
    vcfRNA.tib.tmp$DP.O3 = geno(vcf)$DP[,ID.O3]
  }else{
    vcfRNA.tib.tmp$AD.O3 = NA
    vcfRNA.tib.tmp$DP.O3 = NA
  }
  vcfRNA.tib.tmp$VAF.P  = vcfRNA.tib.tmp$AD.P/vcfRNA.tib.tmp$DP.P
  vcfRNA.tib.tmp$VAF.O1 = vcfRNA.tib.tmp$AD.O1/vcfRNA.tib.tmp$DP.O1
  vcfRNA.tib.tmp$VAF.O2 = vcfRNA.tib.tmp$AD.O2/vcfRNA.tib.tmp$DP.O2
  vcfRNA.tib.tmp$VAF.O3 = vcfRNA.tib.tmp$AD.O3/vcfRNA.tib.tmp$DP.O3
  
  # find number of shared and private with DP >=30
  alt_summary_RNA = bind_rows(alt_summary_RNA, vcfRNA.tib.tmp %>% filter( (is.na(DP.P) | DP.P>=50) & DP.O1>=50 & (is.na(DP.O2) | DP.O2>=50)  ) %>% 
                                summarize( "P&O1&O2&O3"= sum(AD.P>0 & AD.O1>0 & !is.na(AD.O2) & AD.O2>0& !is.na(AD.O3) & AD.O3>0) , 
                                           "P&O1&O2"= sum(AD.P>0 & AD.O1>0 & !is.na(AD.O2) & AD.O2>0 & (is.na(AD.O3 | AD.O3==0)) ) ,
                                           "P&O1&O3"= sum(AD.P>0 & AD.O1>0 & !is.na(AD.O3) & AD.O3>0 & (is.na(AD.O2 | AD.O2==0)) ) , 
                                           "P&O2&O3"= sum(AD.P>0 & AD.O1==0 & !is.na(AD.O2) & AD.O2>0 & !is.na(AD.O3) & AD.O3>0 ) ,
                                           "O1&O2&O3"= sum(AD.P==0 & AD.O1>0 & !is.na(AD.O2) & AD.O2>0 & !is.na(AD.O3) & AD.O3>0 ) , 
                                           "P&O1" = sum(!is.na(AD.P) & AD.P>0 & AD.O1>0 & (is.na(AD.O2) | AD.O2==0) & (is.na(AD.O3) | AD.O3==0) ) , 
                                           "P&O2" = sum(AD.P>0 & AD.O1==0 & !is.na(AD.O2) & AD.O2>0 & (is.na(AD.O3) | AD.O3==0) ) ,
                                           "P&O3" = sum(AD.P>0 & AD.O1==0 & !is.na(AD.O3) & AD.O3>0 & (is.na(AD.O2 | AD.O2==0)) ) ,
                                           "O1&O2"= sum(AD.P==0 & AD.O1>0 & !is.na(AD.O2) & AD.O2>0 & (is.na(AD.O3) | AD.O3==0) ) , 
                                           "O1&O3"= sum(AD.P==0 & AD.O1>0 & !is.na(AD.O3) & AD.O3>0 & (is.na(AD.O2) | AD.O2==0) ) , 
                                           "O2&O3"= sum(AD.P==0 & AD.O1==0 & !is.na(AD.O2) & AD.O2>0 & !is.na(AD.O3) & AD.O3>0 ) , 
                                           P = sum(AD.P>0 & AD.O1==0 & (is.na(AD.O2) | AD.O2==0)) ,
                                           O1 = sum( (AD.P==0 | is.na(AD.P)) & AD.O1>0 & (is.na(AD.O2) | AD.O2==0)) , 
                                           O2 = sum( (AD.P==0 | is.na(AD.P)) & AD.O1==0 & (is.na(AD.O3) | AD.O3==0) ),
                                           O3 = sum( (AD.P==0 | is.na(AD.P)) & AD.O1==0 & (is.na(AD.O2) | AD.O2==0) )) %>% mutate(Experiment=annovar_org_RNA.names[i]) )
  
  vcfRNA.tib.tmp.long = vcfRNA.tib.tmp %>% pivot_longer(cols = AD.P:AD.O3,names_to = c("geno","Sample"), names_sep = "\\.",values_to = "AD") %>% filter(!is.na(AD)) # & AD>0 & )
  vcfRNA.tib.tmp.long$DP = vcfRNA.tib.tmp.long$DP.O3
  vcfRNA.tib.tmp.long$DP[vcfRNA.tib.tmp.long$Sample=="P"] = vcfRNA.tib.tmp.long$DP.P[vcfRNA.tib.tmp.long$Sample=="P"]
  vcfRNA.tib.tmp.long$DP[vcfRNA.tib.tmp.long$Sample=="O1"] = vcfRNA.tib.tmp.long$DP.O1[vcfRNA.tib.tmp.long$Sample=="O1"]
  vcfRNA.tib.tmp.long$DP[vcfRNA.tib.tmp.long$Sample=="O2"] = vcfRNA.tib.tmp.long$DP.O1[vcfRNA.tib.tmp.long$Sample=="O2"]
  # remove not found and low coverage
  vcfRNA.tib.tmp.long$AD[vcfRNA.tib.tmp.long$AD==0 & vcfRNA.tib.tmp.long$DP<10] = -1
  
  vcfRNA.tib.tmp.long$VAF = vcfRNA.tib.tmp.long$VAF.O3
  vcfRNA.tib.tmp.long$VAF[vcfRNA.tib.tmp.long$Sample=="P"]  = vcfRNA.tib.tmp.long$VAF.P[vcfRNA.tib.tmp.long$Sample=="P"]
  vcfRNA.tib.tmp.long$VAF[vcfRNA.tib.tmp.long$Sample=="O1"] = vcfRNA.tib.tmp.long$VAF.O1[vcfRNA.tib.tmp.long$Sample=="O1"]
  vcfRNA.tib.tmp.long$VAF[vcfRNA.tib.tmp.long$Sample=="O2"] = vcfRNA.tib.tmp.long$VAF.O2[vcfRNA.tib.tmp.long$Sample=="O2"]
  vcfRNA.tib.tmp.long = vcfRNA.tib.tmp.long %>% dplyr::select(-(DP.P:VAF.O3),-geno ) %>% mutate(Sample = case_when(Sample=="P"~samplesTimeRNA[annovar_org_RNA.names[i],]$Parental,
                                                                                                             Sample=="O1"~samplesTimeRNA[annovar_org_RNA.names[i],]$PDTO1,
                                                                                                             Sample=="O2"~samplesTimeRNA[annovar_org_RNA.names[i],]$PDTO2,
                                                                                                             TRUE~samplesTimeRNA[annovar_org_RNA.names[i],]$PDTO3) )
  vcfRNA.tib.long = bind_rows(vcfRNA.tib.long , vcfRNA.tib.tmp.long)
}

# Switch to recognized annovar values
vcfRNA.tib.long$ExonicFunc.ensGene = str_replace(vcfRNA.tib.long$ExonicFunc.ensGene, "_", " ")
vcfRNA.tib.long$non_cancer_AF_popmax[is.na(vcfRNA.tib.long$non_cancer_AF_popmax)] = 0
vcfRNA.tib.long$ExAC_nontcga_ALL[is.na(vcfRNA.tib.long$ExAC_nontcga_ALL)] = 0
vcfRNA.tib.long$`1000g2015aug_all`[is.na(vcfRNA.tib.long$`1000g2015aug_all`)] = 0

all(vcfRNA.tib.long$ExonicFunc.ensGene[vcfRNA.tib.long$ExonicFunc.ensGene!="."] %in% names(annovar_values))
## example 
vcfRNA.tib.long%>% filter(Gene.ensGene=="TP53", !ExonicFunc.ensGene %in% c("synonymous SNV","."), cosmic92_coding!=".", non_cancer_AF_popmax<=0.01,
                          ExAC_nontcga_ALL<=0.01 , `1000g2015aug_all`<=0.01 , REVEL>=0.5 |ExonicFunc.ensGene!="nonsynonymous SNV" )  %>% 
  dplyr::select(Start,Ref,Alt,ExonicFunc.ensGene,REVEL,cosmic92_coding,Sample,AD,DP,VAF)

vcfRNA.tib.long%>% filter(Gene.ensGene=="RB1", !ExonicFunc.ensGene %in% c("synonymous SNV","."), #cosmic92_coding!=".", 
                          non_cancer_AF_popmax<=0.01,
                          ExAC_nontcga_ALL<=0.01 , `1000g2015aug_all`<=0.01 , REVEL>=0.5 |ExonicFunc.ensGene!="nonsynonymous SNV" )  %>% 
  dplyr::select(Start,Ref,Alt,ExonicFunc.ensGene,REVEL,cosmic92_coding,Sample,AD,DP,VAF)

vcfRNA.tib.long%>% filter(Gene.ensGene=="STK11", !ExonicFunc.ensGene %in% c("synonymous SNV","."), #cosmic92_coding!=".", 
                          non_cancer_AF_popmax<=0.01,
                          ExAC_nontcga_ALL<=0.01 , `1000g2015aug_all`<=0.01 , REVEL>=0.5 |ExonicFunc.ensGene!="nonsynonymous SNV" )  %>% 
  dplyr::select(Start,Ref,Alt,ExonicFunc.ensGene,REVEL,cosmic92_coding,Sample,AD,DP,VAF)

write_tsv(alt_summary_RNA , file = "/data/lungNENomics/work/organoids/integration/mutations_venn_RNA.tsv")
write_tsv(vcfRNA.tib.long , file = "/data/lungNENomics/work/organoids/integration/multisample_annovar_RNA.maf")

cbind(alt_summary_RNA$Experiment , 1 - alt_summary_RNA$P / (rowSums(alt_summary_RNA[,c(1:4,6:8)]) + alt_summary_RNA$P ))

# load file
organoids_somaticRNA = annovarToMaf(annovar = "/data/lungNENomics/work/organoids/integration/multisample_annovar_RNA.maf", 
                                 Center = 'Utrecht', refBuild = 'hg38', tsbCol = 'Sample', table = 'ensGene',MAFobj = TRUE,ens2hugo = FALSE,)

## test if all in exonic
table(organoids_somaticRNA@data$ExonicFunc.refGene, organoids_somaticRNA@data$Variant_Classification )
organoids_somaticRNA@data %>% filter( ExonicFunc.refGene=="." )
table(organoids_somaticRNA@maf.silent$ExonicFunc.refGene, organoids_somaticRNA@maf.silent$Variant_Classification)
organoids_somaticRNA@maf.silent %>% filter( !ExonicFunc.refGene %in% c(".","synonymous SNV") )
### OK
plotmafSummary(organoids_somaticRNA)

# write to file
mafRNA = organoids_somaticRNA@data %>% filter(Variant_Classification!="Splice_Site" )
mafRNA$Variant_Classification = as.character(mafRNA$Variant_Classification)
mafRNA$Variant_Classification[mafRNA$AD==-1] = "No_Coverage"
write_tsv(  mafRNA ,file = "/data/lungNENomics/work/organoids/integration/multisample_annovar_RNA2.maf"  )


## add fusions # to do
maf.data = organoids_somaticRNA@data

SV.maf = driver_uniq_recovered %>% filter(ID!="") %>% pivot_longer(cols = Parental:PDTO2,names_to = "Sample") %>% filter(!is.na(value) & value!="NO" )
SV.maf$Sample = apply(SV.maf ,1 , function(x) samplesTime[x$Experiment,x$Sample])

SV.maf = SV.maf %>% mutate(Chromosome= CHROM, Start_Position=start.B1 , End_Position=end.B1 , Reference_Allele = NA , Tumor_Seq_Allele2 = NA , Tumor_Sample_Barcode = Sample , 
                           Hugo_Symbol = Gene, Variant_Classification =  SVTYPE, tx=NA , exon=NA , txChange=NA , 
                           aaChange=NA , Variant_Type = "SV" , sample_id = "multisample_SV", Func.refGene ="exonic", Gene.refGene=Gene , 
                           GeneDetail.refGene =ID , ExonicFunc.refGene = SVTYPE , AAChange.refGene =NA , REVEL =NA , AD =NA,
                           DP=NA , VAF=NA) %>% dplyr::select(Chromosome:VAF)
SV.maf$Hugo_Symbol = sapply(SV.maf$Gene.refGene,function(x){ res=NA; if(length(x)>0){res=as.character(x)} ; return(res) })
SV.maf$Gene.refGene = SV.maf$Hugo_Symbol

SV.maf = SV.maf[which(!duplicated(SV.maf[,c(1,6:8,13:14)])),]
maf.data = rbind(maf.data,SV.maf)

# check if damaging
maf.data = maf.data %>% mutate(damaging = case_when(as.numeric(REVEL)>=0.5 | Variant_Classification !="Missense_Mutation" ~"YES",
                                                    as.numeric(REVEL)<0.5 | Variant_Classification !="Missense_Mutation" ~"MAYBE",
                                                    TRUE~"NO") )

write_tsv( maf.data ,file = "/data/lungNENomics/work/organoids/WGS/variant_calling/release3_25102021/multisample_annovar_wSV.maf"  )

# read maf
organoids_somaticRNA2 = read.maf("/data/lungNENomics/work/organoids/integration/multisample_annovar_RNA2.maf" ,
                              vc_nonSyn = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", 
                                            "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","No_Coverage"))

## plot 

#sampleOrder_nmut = sampleOrder[c(7:12,5:6,3:4,15:16,13:14,17:18,1:2)]
sampleOrder_onco = c(sampleOrder[17:22],sampleOrder[1:16])
levels_normalRNA = c("SINET21","SINET22","SINET7","SINET8","SINET12","LNET6","LNET13","LNET14","LNET19",
                  "LNET5","LNET10","LNET15","LNET16","LNET18","LNET20",
                  "LCNEC3","LCNEC4","LCNEC11","LCNEC23","PANEC1")
#levels_nmut = c("LCNEC4", "PANEC1", "LCNEC3", "LNET10", "SINET8", "SINET7", "SINET9", "LNET6")

organoids_somaticRNA2sub = subsetMaf(organoids_somaticRNA2,
                                     query = "(Gene.refGene=='MEN1' & (aaChange=='p.D33V' |aaChange=='p.V436M')) | (non_cancer_AF_popmax<=0.01 & ExAC_nontcga_ALL<=0.01 & `1000g2015aug_all`<=0.01 & cosmic92_coding !='.' & (REVEL>=0.7 | is.na(REVEL)) )")
tabgenesRNA = organoids_somaticRNA2sub@data %>% filter(Hugo_Symbol%in%mutated_gene_listB$`gene name` ) %>% 
  mutate(Exp=str_remove(Tumor_Sample_Barcode,"[MT]p*[0-9.]*$")) %>% 
  group_by(Hugo_Symbol) %>% summarise(n=n(),nexp=length(unique(Exp)),Exp ) %>% left_join(mutated_gene_listB,by=c("Hugo_Symbol"="gene name")) %>% 
  mutate(WhichType=(WhichType==2)*(-1.5)+WhichType,Type1=WhichType==1,Type2=WhichType==2 ) %>% 
  mutate(Exp=factor(Exp,levels =  levels_normalRNA)) %>% 
  arrange(Exp,desc(nexp),desc(n))
#tabgenes$Exp = 
#tabgenes$WhichType[tabgenes$WhichType==2] = 0.5

vc_col = maftools:::get_vcColors(websafe = FALSE)
#vc_col["INV"] = "#a02c89ff"
#vc_col["TRA"] = "#42b6f5ff"
#vc_col["DEL"] = "#d647b9ff"
#vc_col["DUP"] = "#5c47d6ff"
vc_col["No_Coverage"] = rgb(0.5,0.5,0.5)
vc_col["Missense_Mutation"] = "#1f78b4ff"
vc_col["Frame_Shift_Del"] = "#782121ff"
vc_col["Frame_Shift_Ins"] = "#2ca089ff"

pdf("/data/lungNENomics/work/organoids/figures/FigS4C_raw_18122021.pdf",h=12,w=8)
oncoplot(maf = organoids_somaticRNA2sub, showTumorSampleBarcodes = T, barcode_mar = 6, drawRowBar = F,drawColBar = F,  showTitle = F,
         #additionalFeature = c("dam_pred075",TRUE), 
         leftBarData = tabgenesRNA[!duplicated(tabgenesRNA$Hugo_Symbol),c(1,6)], 
         sampleOrder = sampleOrderRNA,#sampleOrder_onco,
         genes = unique(tabgenesRNA$Hugo_Symbol), 
         keepGeneOrder = TRUE,
         removeNonMutated = FALSE,
         colors = vc_col,
         borderCol= "white",
         bgCol = "#fef8e4ff")
dev.off()

pdf("/data/lungNENomics/work/organoids/figures/Fig4C_raw_18122021.pdf",h=12,w=8)
oncoplot(maf = organoids_somaticRNA2sub, showTumorSampleBarcodes = T, barcode_mar = 6, drawRowBar = F,drawColBar = F,  showTitle = F,
         #additionalFeature = c("dam_pred075",TRUE), 
         leftBarData = tabgenesRNA[!duplicated(tabgenesRNA$Hugo_Symbol),c(1,6)], 
         sampleOrder = sampleOrderRNA[!str_detect(sampleOrderRNA,"SINET7|SINET8|LNET6|LNET5|LNET10|LCNEC3|LCNEC4|PANEC1|LNET16T|LNET18|LNET19|LNET20")],#sampleOrder_onco,
         genes = unique(tabgenesRNA$Hugo_Symbol[!str_detect(tabgenesRNA$Exp,"SINET7|SINET8|LNET6|LNET5|LNET10|LCNEC3|LCNEC4|PANEC1|LNET18|LNET19|LNET20")]), 
         keepGeneOrder = TRUE,
         removeNonMutated = TRUE,
         colors = vc_col,
         borderCol= "white",
         bgCol = "#fef8e4ff")
dev.off()


# write data
write_tsv( organoids_somatic2.dam@data ,path = "/data/lungNENomics/work/organoids/figures/Fig4C_draft_data_18122021.tsv" )
write_tsv( organoids_somaticRNA2sub@data ,path = "/data/lungNENomics/work/organoids/figures/FigS4C_draft_data_18122021.tsv" )


#######  S figs
# corrected TMB
predcod_exonic.Allgenes.longmat <- predcod_exonic2 %>% 
  group_by(Sample,QUERYID) %>% count(var_predict,variant_class) %>% dplyr::slice(which(n==max(n)) ) %>% 
  dplyr::slice(which.min(as.numeric(var_predict)) ) 

TMBnew = predcod_exonic.Allgenes.longmat %>% filter(var_predict!="Silent") %>% group_by(Sample) %>% count()

# load cosmic genes
cosmic_genes = read_xlsx("/run/user/1000/gvfs/smb-share:server=shares.hpc.iarc.lan,share=data/mesomics/work/mangiantel/DataBases/COSMIC_cancer_genes_v90.xlsx")
cosmic_hallmarks = read_xlsx("/run/user/1000/gvfs/smb-share:server=shares.hpc.iarc.lan,share=data/mesomics/work/mangiantel/DataBases/COSMIC_cancer_genes_v90.xlsx",sheet = 2,skip = 17)
cosmic_hallmarks_head = read_xlsx("/run/user/1000/gvfs/smb-share:server=shares.hpc.iarc.lan,share=data/mesomics/work/mangiantel/DataBases/COSMIC_cancer_genes_v90.xlsx",sheet = 2,n_max = 12,skip=3)

cosmic_gene_list = unlist(sapply(cosmic_genes$Synonyms,strsplit,","))

#  mutation load
source("TMB_compare2.R")
MAF.list = list( #subsetMaf( organoids_somatic,tsb = grep("LNET2T",levels(organoids_somatic@data$Tumor_Sample_Barcode),value = T),ranges = ss50.hg38 ),
  subsetMaf( organoids_somatic,tsb = grep("LNET6",levels(organoids_somatic@data$Tumor_Sample_Barcode),value = T),ranges = ss50.hg38 ),
  subsetMaf( organoids_somatic,tsb = grep("LNET10",levels(organoids_somatic@data$Tumor_Sample_Barcode),value = T),ranges = ss50.hg38 ),
  subsetMaf( organoids_somatic,tsb = grep("LCNEC3",levels(organoids_somatic@data$Tumor_Sample_Barcode),value = T),ranges = ss50.hg38 ),
  subsetMaf( organoids_somatic,tsb = grep("LCNEC4",levels(organoids_somatic@data$Tumor_Sample_Barcode),value = T),ranges = ss50.hg38 ),
  subsetMaf( organoids_somatic,tsb = grep("PANEC1",levels(organoids_somatic@data$Tumor_Sample_Barcode),value = T),ranges = ss50.hg38 ),
  subsetMaf( organoids_somatic,tsb = grep("SINET7",levels(organoids_somatic@data$Tumor_Sample_Barcode),value = T),ranges = ss50.hg38 ),
  subsetMaf( organoids_somatic,tsb = grep("SINET8",levels(organoids_somatic@data$Tumor_Sample_Barcode),value = T),ranges = ss50.hg38 ),
  subsetMaf( organoids_somatic,tsb = grep("SINET9",levels(organoids_somatic@data$Tumor_Sample_Barcode),value = T),ranges = ss50.hg38 ))
TMBnew2 = bind_rows(lapply( MAF.list , function(x) x@variants.per.sample))

mut_LNEN = read_xlsx("41467_2019_11276_MOESM7_ESM.xlsx",skip = 49)
mut_SCLC = read_xlsx("41586_2015_BFnature14664_MOESM72_ESM.xlsx",sheet = 3,skip=2)
library(data.table)
other.cohort = mut_LNEN %>% group_by(Sample_ID,Histopathology) %>% summarize(total=n()) %>% mutate(Histopathology=fct_collapse(Histopathology,LNET=c("Atypical","Typical")) ) %>% as.data.table
other.cohort$site = "Lung"
colnames(other.cohort)[1:2] = c("Tumor_Sample_Barcode","cohort")
other.cohort2 = mut_SCLC %>% group_by(PAT_ID) %>% summarize(total=n()) %>% as.data.table
other.cohort2$site = "Lung"
other.cohort2$cohort = "SCLC"
colnames(other.cohort2)[1] = c("Tumor_Sample_Barcode")
tcga.cohort2 = rbind(tcga.cohort %>% filter(cohort %in%c("LUAD","LUSC")),other.cohort,other.cohort2 )#,c("","Average",30))

TMB.dat = tcga.cohort2 %>% filter(cohort%in%c("LCNEC","SCLC","LNET","LUAD","LUSC") )
# add SINET TMB
require(janitor)
TMB.SINETref = read_xlsx("TML Talya.xlsx")
TMB.SINETref = TMB.SINETref %>% clean_names() %>% mutate(tumor_mutational_burden=as.numeric(str_remove(tumor_mutational_burden," variants per mb") ) ) %>% 
  rename( "Tumor_Sample_Barcode"="id_no",  "total"="tumor_mutational_burden"  )

TMB.dat = TMBnew2 %>% mutate(cohort="Organoids",site="depends") %>% rename(total=Variants)  %>% bind_rows(TMB.dat)  # organoids_somatic@variants.per.sample
TMB.dat$Exp = "Reference"
TMB.dat$Exp[TMB.dat$cohort =="Organoids"] = gsub(".$","",gsub("p.+$","",TMB.dat$Tumor_Sample_Barcode[TMB.dat$cohort =="Organoids"]) )
TMB.dat$cohort = as.character(TMB.dat$cohort)
TMB.dat$cohort[TMB.dat$cohort =="Organoids"] = gsub("[0-9]+","",TMB.dat$Exp[TMB.dat$cohort =="Organoids"])
TMB.dat$cohort = as.factor(TMB.dat$cohort)
TMB.dat = rbind(TMB.dat %>% ungroup() , cbind(TMB.SINETref[,1], TMB.SINETref[,3]*35.8,cohort="SINET",site="SI",Exp="Reference") )


TMB.dat = TMB.dat %>% filter(!Exp%in%c("LNET2","LNET5"),Tumor_Sample_Barcode!="LNET2Np12") %>%  
  mutate(Tumor_Sample_Barcode=fct_drop(fct_reorder(Tumor_Sample_Barcode,total)),
         cohort=fct_drop(fct_reorder(cohort,total,median)) )

tp <- data.frame(cohort=levels(TMB.dat$cohort),val=rep(c("A","B"),length.out=length(levels(TMB.dat$cohort)) ) )


traject_RNA <- function(dat=expr_genes_NET,samples,col=colors_org[3],xpos=c(1:11,1:2)-0.20,geneset=c("Neuroendocrine","NET subtype"),xoff=-0.001){
  geom_curve(data = left_join(dat[dat$Sample==samples[1],],dat[dat$Sample==samples[2],],
                              by=c("ID","Experiment","Set"),suffix=c("",".y"))%>% filter(Set%in%geneset) %>% mutate(xpos=xpos), 
             aes(x =  xpos, y = value , xend = xpos+xoff,yend = value.y), size=0.5,colour = col,inherit.aes = T ,curvature = -0.2,lineend = "round")
}
traject_RNA_head <- function(dat=expr_genes_NET,samples,col=colors_org[3],xpos=c(1:11,1:2)-0.20,geneset=c("Neuroendocrine","NET subtype"),xoff=-0.001,
                             fills=c(col,"white")){
  geom_point(data = dat[sapply(samples, function(x) which(dat$Sample==x)),]%>% filter(Set%in%geneset), 
             aes(x =  c(xpos,xpos+xoff), y = value), size=2,stroke=1,colour = col,fill=rep(fills,each=length(xpos)),shape=21,inherit.aes = T)
}

FigS4A <- ggplot( TMB.dat,aes(x=Tumor_Sample_Barcode,y=total/35.8,col=Exp)) +
  geom_rect(data = tp,aes(fill = val),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.5,inherit.aes = F,show.legend = F) +
  scale_fill_manual(values = c("white", "#ffe6d5ff"))+ #new_scale_color() +
  geom_point(size=0.7) + scale_color_manual(values = c(Reference="gray76",colors_org_WGS[-2]), limits = c("Reference",names(colors_org_WGS[-2])))+
  geom_hline(data=TMB.dat%>% filter(Exp=="Reference") %>% group_by(cohort) %>% summarize(median=median(total)/35.8), mapping = aes(yintercept=median), col="red" ) +
  geom_path(data=TMB.dat[1:18,] %>% mutate(Tumor_Sample_Barcode=fct_relevel(fct_drop(Tumor_Sample_Barcode),sampleOrder) )  %>% 
              arrange(Tumor_Sample_Barcode) ,aes(group=Exp),size=1.2,lineend = "round")+ new_scale_fill() +
  geom_point(data=TMB.dat[1:18,][!str_detect(TMB.dat[1:18,]$Tumor_Sample_Barcode,"p"),],aes(fill=Exp),size=2,stroke=1,shape=21) +
  geom_point(data=TMB.dat[1:18,][str_detect(TMB.dat[1:18,]$Tumor_Sample_Barcode,"p"),],size=2,stroke=1,fill="white",shape=21) + 
  scale_fill_manual(values = colors_org)+ 
  scale_y_log10() + 
  facet_wrap(.~cohort,scales = "free_x",ncol = 16,strip.position = "bottom",shrink = F) + 
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "gray76"),
        legend.key = element_blank(), legend.title = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face=c(rep("plain",2),"bold",rep("plain",3)),angle = 90,hjust = 1,size = 12),
        axis.title.x = element_blank()) + labs(y="TMB (nonsynon. mutation / Mb)") + coord_cartesian(clip = "off") + guides(fill=F,color=F)
Fig4A
ggsave("Fig4A_raw_01032021.svg",Fig4A,height = 3,width = 3*2)
ggsave("Fig4A_raw_01032021.png",Fig4A,height = 3,width = 3*2)
ggsave("Fig4A_raw_01032021.pdf",Fig4A,height = 3,width = 3*2)

write_tsv(Fig4A$data, "Fig4A_draft_data_01032021.tsv")

