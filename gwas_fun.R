#output gwas_re
load("/root/RNA_test/TPM_test.RData")#gwaa data_dir
library(parallel)
library(GenABEL.data)
library(GenABEL)
library(stringr)
data <- tpm_a
phe_sort <- phdata(tpm_a)[,c(3:ncol(phdata(tpm_a)))]
col_names <- colnames(phe_sort)
thres <- -log10(0.05/nsnps(data))#use java gec.jar to calculate will be beter 
#calculate effective snps numbers
#java -jar gec.jar -Xmx1g --effect-number --maf 0.03 --no-web --plink-binary bed --genome --out 
thres2 <- 4
kin <- ibs(data,weight = "freq")
export.plink(data,"/w/00/g/g01/user319/data_gwas/tpm")#output plink 
GWAS <- function( i ){
  v_now <- phe_sort[,i]
  cc<- try(polygenic(v_now,kinship.matrix = kin,data = data),silent = T)
  if(!(class(cc)=="try-error")){ 
    cat(i,"pass","\n")
  } else cat(i ,"fail","\n")
  mm <- try(mmscore(cc,data=data),silent = T)
  if(!(class(mm) == "try-error")){
    cat(i,"pass","\n")
  } else cat(i,"fail","\n")
  if(any(-log10(mm[,"Pc1df"]) >= thres2)){ #can change to get different results
    options(bitmapType='cairo')
    png(paste0("/w/00/g/g01/user319/plot/",i,".png"))#a dir to put plotpic in 
    par(mfrow=c(1,3))
    hist(phe_sort[,i],main = i)
    if (mm@lambda$estimate <= 1.1){
      plot(mm)
      abline(h=thres,lty="dashed",col="red")
      abline(h=thres2,lty="dashed",col="red")
      p <- mm@results$P1df
      chr <- as.character(mm@annotation$Chromosome)
      pos <- as.numeric(str_split_fixed(rownames(mm), ":",2)[,2])
      gwas_result <- data.frame(chr = chr,position = pos,Pvalue = p)
      write.table(gwas_result,paste0("/root/RNA_test/",i,".Rdata",".txt",".txt"),row.names = F,quote = F)#dir to put gwasresults in 
    } else{
      plot(mm,df = "Pc1df")
      abline(h=thres,lty="dashed",col="red")
      abline(h=thres2,lty="dashed",col="red")
      p <- mm@results$Pc1df
      chr <- as.character(mm@annotation$Chromosome)
      pos <- as.numeric(str_split_fixed(rownames(mm), ":",2)[,2])
      gwas_result <- data.frame(chr = chr,position = pos,Pvalue = p)
      write.table(gwas_result,paste0("/root/RNA_test/",i,".Rdata",".txt",".txt"),row.names = F,quote = F)#dir to put gwasresults in 
    }
    b <- estlambda(mm[,"P1df"],plot = T)
    garbage <- dev.off()
    save(cc,mm,file=paste0("/root/RNA_test/scan/",i,".Rdata"))#dir to put RData 
  } else{
      return(0)
  }
  return(0) 
}
r <- mclapply(col_names,FUN = GWAS, mc.cores = 10) 
#if use this pipeline in windows
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum))
clusterEvalQ(cl,library(GenABEL))
clusterEvalQ(cl,library(stringr))
clusterEvalQ(cl,library(GenABEL.data))
clusterExport(cl,c("phe_sort","data","kin","thres","thres2","col_names"))
r <- parLapply(cl,col_names,fun = GWAS)
stopCluster(cl)
 
