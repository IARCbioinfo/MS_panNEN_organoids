###################################################
# Script to preprocess mutation data for DPclust. #
# Adapted from DPclust3p github repo              #
###################################################

## Load libraries
library(dpclust3p)
library(DPClust)
library(optparse)
library(tidyverse)

# create useful function
PurpleToBattenberg <- function (samplename, outfile.folder, copynumber.purple.file, samplestatistics.file){
  d = read_xlsx(copynumber.purple.file , sheet=6,skip=2) %>% filter(sample==samplename) %>% 
    mutate(minorAlleleCopyNumber = as.numeric(minorAlleleCopyNumber),majorAlleleCopyNumber = as.numeric(majorAlleleCopyNumber))
  # keep only clonal CN, i.d. with fractional part less than 0.2 or greater than 0.8
  d= d %>% filter(minorAlleleCopyNumber %%1<=0.2 | minorAlleleCopyNumber %%1>=0.8,
               majorAlleleCopyNumber %%1<=0.2 | majorAlleleCopyNumber %%1>=0.8)
  
  subclones = d[, c(2:4)] # chr, start, end
  colnames(subclones) = c("chr", "startpos", "endpos")
  subclones$normal_total = 2
  if(sex=="male") subclones$normal_total[subclones$chr%in%c("chrX","chrY")] = 1
  subclones$nMaj1_A = round(d$majorAlleleCopyNumber)
  subclones$nMin1_A = round(d$minorAlleleCopyNumber)
  subclones$frac1_A = 1
  subclones$nMaj2_A = NA
  subclones$nMin2_A = NA
  subclones$frac2_A = NA
  
  subclones = subclones %>% filter(!is.na(frac1_A))# %>% mutate(chr=str_remove(chr,"chr"))
  
  write_tsv(subclones, paste(outfile.folder, samplename, "_subclones.txt", sep = ""))
  
  samplestats = read_xlsx(copynumber.purple.file , sheet=5,skip=2) %>% filter(tumor_id==samplename)
  purity_ploidy = array(NA, c(3, 5))
  colnames(purity_ploidy) = c("rho", "psi", "ploidy", "distance", "is.best")
  rownames(purity_ploidy) = c("ASCAT", "FRAC_GENOME", "REF_SEG")
  purity_ploidy["FRAC_GENOME", "rho"] = as.numeric(samplestats$purity )
  purity_ploidy["FRAC_GENOME", "psi"] = as.numeric(samplestats$ploidy)
  purity_ploidy["FRAC_GENOME", "ploidy"] = as.numeric(samplestats$ploidy)
  purity_ploidy["FRAC_GENOME", "is.best"] = TRUE
  write.table(purity_ploidy, paste(outfile.folder, samplename, "_rho_and_psi.txt", sep = ""), quote = F, sep = "\t")
}


# get options
if(!exists("opt")){
option_list = list(
  make_option(c("-s", "--samplename"), type="character", default=NULL, help="Samplename of the sample to run", metavar="character"),
  make_option(c("-a", "--alignfolder"), type="character", default=NULL, help="Folder with BAM/CRAM files to extract mutation allele counts from, there must be an index with the same name, but with an extra .bai or .crai extension", metavar="character"),
  make_option(c("-v", "--vcf"), type="character", default=NULL, help="VCF file with mutation data", metavar="character"),
  make_option(c("--samplestats"), type="character", default=NULL, help="Facets samplestatistics output file", metavar="character"),make_option(c("-c", "--copynumber"), type="character", default=NULL, help="Facets output file with copy number data", metavar="character"),
  make_option(c("--sex"), type="character", default=NULL, help="Sex of the sample", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output directory", metavar="character"),
  make_option(c("--fai"), type="character", default=NULL, help="Reference genome index", metavar="character"),
  make_option(c("--ign"), type="character", default=NULL, help="File with a list of contigs to ignore", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
}

samplename  = opt$samplename
align_files = grep(list.files(path=opt$alignfolder,pattern = paste0(samplename)),pattern = "ai$",value = T,invert = T)
samplenames = sapply(align_files, function(x) str_split(x,pattern = "_")[[1]][1] )
vcf_file    = opt$vcf
copynumber.purple.files = opt$copynumber
sex = opt$sex
output_dir = opt$output
fai_file   = opt$fai
ign_file   = opt$ign

# check if files exist
.checkfile = function(infile) {
  if (!file.exists(infile)) {
    stop(paste("File", infile, "does not exist", sep=""))
  }
}

.checkfile(vcf_file)
.checkfile(fai_file)
.checkfile(ign_file)

if (!sex %in% c("male", "female")) {
  stop("Provide male or female as sex")
}

# Define the final output file
dpoutput_files = file.path(output_dir, paste(samplenames, "_allDirichletProcessInfo.txt", sep=""))

# Define various temp files
loci_file = file.path(output_dir, paste(samplename, "_loci.txt", sep=""))
allelecounts_files = file.path(output_dir, paste0(samplenames, "_alleleFrequencies.txt", sep=""))

# Dump loci - this function can take multiple vcf files when multiple samples from same donor
vcf2loci(vcf_file=vcf_file, fai_file=fai_file, ign_file=ign_file, outfile=loci_file)

for(i in 1:length(allelecounts_files)){
  # Fetch allele counts
  alleleCount(locifile=loci_file, bam=paste0(opt$alignfolder,align_files[i]), outfile=allelecounts_files[i], min_baq=20, min_maq=35)
  # Transform ASCAT NGS output into Battenberg output
  PurpleToBattenberg(samplenames[i], output_dir, copynumber.purple.files, samplestatistics.file)
  # Create dpIn file
  runGetDirichletProcessInfo(loci_file=loci_file, allele_frequencies_file=allelecounts_files[i], 
                             cellularity_file=paste(output_dir, samplenames[i], "_rho_and_psi.txt", sep = ""), 
                             subclone_file=paste(output_dir, samplenames[i], "_subclones.txt", sep = ""), 
                             gender=sex, SNP.phase.file="NA", mut.phase.file="NA", output_file=dpoutput_files[i])
  
}
