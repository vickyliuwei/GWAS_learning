setwd("/root/RNA_test/peak/")
dir_list <- dir("/root/RNA_test/peak/")#load peak
library(parallel)
library(GenABEL.data)
library(GenABEL)
library(stringr)
plink_path <- "/w/00/g/g01/user319/data_gwas/ft_test"
gff_path <- "///.gff"
read_regon <- function(i){
    if(file.info(dir_list[i])$size != 0){
      gwas_path <- paste0("path",i,".txt")
      peak_list <- read.table(dir_list[i])
      peak_list <- as.character(peak_list$V1)
      seq <- 10000
    for (i in 1:length(peak_list)){
        a <- str_split_fixed(peak_list[i],":",2)[,1]
        top <- str_split_fixed(peak_list[i],":",2)[,2]
        start <- as.character(as.numeric(top) - seq)
        end <- as.character(as.numeric(top) + seq)
        Locuszoom(LDBlockShow = LDBlockShow,Plink = plink_path,GWAS = gwas_path,GFF = gff_path,Region = paste0(a,":",start,":",end),Output = paste0("/w/00/g/g01/user319/locus/",dir_list[i],a,":",start,":",end))
    }
  } else{
      return(0)
  }
  return(0)
}
Locuszoom <- function(LDBlockShow,Plink,GWAS,GFF,Region,Output){
    cmd_locus = paste(LDBlockShow,"-InPlink",Plink,"-InGWAS",GWAS,"-InGFF",GFF,"-Region",Region,"-Output",Output,"-OutPng","-SeleVar 2","-TopSite",sep = " ")
    cat(cmd_locus,"\n")
    system(cmd_locus) 
}
n <- 1:length(dir_list)
r <- mclapply(n,FUN = read_regon, mc.cores = 10)
