setwd("/root/RNA_test/scan/")#setwd
dir_list <- dir("/root/RNA_test/scan/")#load RData
library(parallel)
library(GenABEL.data)
library(GenABEL)
library(stringr)
myleading <- function(i){  
  load(dir_list[i])
  a <- dir_list[i]
  if (mm@lambda$estimate <= 1.1){
    mm <- mm[which(-log10(mm[,"P1df"])>= thres),]#calculate thres before this 
  }else{
    mm <- mm[which(-log10(mm[,"Pc1df"])>= thres),]
  }
  chr <- as.numeric(unique((mm[,"Chromosome"])))
  snp <- c()
  for ( i in chr){
    chr.id <- as.numeric(((mm[,"Chromosome"])))
    pos_2 <- as.numeric(str_split_fixed(rownames(mm[chr.id %in% i,]), ":",2)[,2])
    snpname <- rownames(mm[chr.id %in% i,])
    window <- seq(pos_2[1],pos_2[length(pos_2)],0.05e06)
    inter.id<- findInterval(pos_2,window)
    id.loop<- as.numeric(unique(inter.id))
    if(length(id.loop)==2 & length(snpname)==2 & inter.id[1] != inter.id[2]){
      snp <- c(snp,snpname)
    }else if (length(id.loop)==1){
      snp <- c(snp,snpname[which.min(mm[chr.id %in% i,]$P1df)])
    }else{
      for(j in 1:length(id.loop)){
        snp <- c(snp,snpname[inter.id %in% id.loop[j]][which.min(mm[chr.id %in% i,]$P1df[inter.id %in% id.loop[j]])])
      }
    }
  }
  write.table(snp,paste0("/root/RNA_test/peak/",a,".txt"),col.names = F,row.names = F,quote = F)#ouput peak information
}
n <- 1:length(dir_list)
r <- mclapply(n,FUN = myleading,mc.cores = 10)
#if use this pipeline in win
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum))
clusterEvalQ(cl,library(GenABEL))
clusterEvalQ(cl,library(stringr))
clusterEvalQ(cl,library(GenABEL.data))
clusterExport(cl,c("dir_list","n"))
r <- parLapply(cl,col_names,fun = GWAS)
stopCluster(cl)
