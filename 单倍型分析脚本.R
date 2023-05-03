geno_sub <- read.table("/w/00/g/g01/user319/test.012")
geno_1 <- geno_sub[,-1]
geno_1 <- geno_1[,seq(1,ncol(geno_1),10)]
snp_1 <- pos[seq(1,nrow(pos),10),2]
hr <- hclust(as.dist(1-cor(t(geno_1),use = "pairwise.complete.obs")), method="complete")
mycl <- cutree(hr, k = 2)
myCol = c("pink1", "violet")
clusterCols <- rainbow(length(unique(mycl))) #sample(myCol,length(unique(mycl)))#
myClusterSideBar <- clusterCols[mycl]
library(gplots)
heatmap.2(geno_1, main="", 
          Rowv=as.dendrogram(hr), Colv = FALSE,scale="row", 
          col=c("slateblue4", "olivedrab"), density.info="none", trace="none", 
          RowSideColors= myClusterSideBar,cexCol=0.8,srtCol = 45,
          key.title = NULL,labCol = snp_1,labRow = ind$V1)
#can also sort indviduals mannually
geno_1 <- geno_1[sort(geno_1[,i]),]

