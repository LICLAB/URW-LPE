
##########################Input##################################
############an extended lncRNA-PCG co-expression network
#############differential gene list with FC that gets from DESeq software output
############all unique lncRNA list from gencode database
#############permutation times
#################################################################


#setwd("H:\\lienmin\\RNAseq¿ÎÌâ7\\paper\\Integrative analyses of transcriptome sequencing-Oncogene -gai\\revised\\URW-LPE")
#####Compute URWScore##############
library(igraph)
FCvalue<-read.table("./DESeqMMNewSample.txt",header=T,stringsAsFactors=FALSE) #differential gene list with logFC
Network<-read.table("./mergecoexpHPRDLNCdis.txt",header=F,sep="\t",quote="",stringsAsFactors=FALSE) # an extended lncRNA-PCG co-expression network.
lncRNA<-read.table("lncgenenew.txt",header=F,sep="\t",quote="",stringsAsFactors=FALSE) # all unique lncRNA list from gencode database.
permutation<-5  # should set permutation= 2000 or other large number.
source("RandomCoexp.r")  #Random Walk method
allNetRe<-GetESCCscoreCutNet(Network,FCvalue)  #Unsupervised Random Walk method
coexpNetScore_final<-allNetRe
allNetRe<-allNetRe[order(allNetRe[,3],decreasing =TRUE),]
lncScorenew<-allNetRe[which(allNetRe[,1]%in%lncRNA[,1]),]
lncScorenew<-lncScorenew[order(lncScorenew[,3],decreasing=TRUE),]

#######permutation  #########
g <- graph.data.frame(Network[,1:2], directed=F, vertices=NULL)
vetextnames<-get.vertex.attribute(g, "name", index=V(g))
coexpNetScore_final<-coexpNetScore_final[match(vetextnames,coexpNetScore_final[,1]),]
index1_final<-which(coexpNetScore_final[,1]%in%lncRNA[,1])
coexpNetScore_lncRNA_final<-coexpNetScore_final[index1_final,]
index2_final<-which(!coexpNetScore_final[,1]%in%lncRNA[,1])
coexpNetScore_gene_final<-coexpNetScore_final[index2_final,]

coexpNetScore<-FCvalue
index1<-which(coexpNetScore[,1]%in%lncRNA[,1])
coexpNetScore_lncRNA<-coexpNetScore[index1,]
index2<-which(!coexpNetScore[,1]%in%lncRNA[,1])
coexpNetScore_gene<-coexpNetScore[index2,]

logfold<-coexpNetScore_lncRNA[,6]
logfold[is.na(logfold)]<-0
logfold<-logfold[!(logfold==Inf)]
logfold<-logfold[!(logfold==(-Inf))]
null_dis_lncRNA<-logfold
lncRNA_per<-c()
#for(i in 1:2000){
for(i in 1:permutation){ 
coexpNetScore1<-coexpNetScore_final
permutation1<-rnorm(length(coexpNetScore_lncRNA_final[,2]),mean(null_dis_lncRNA),sd(null_dis_lncRNA))
permutation2<-sample(coexpNetScore_gene_final[,2],length(coexpNetScore_gene_final[,2]))
coexpNetScore1[index1_final,2]<-permutation1
coexpNetScore1[index2_final,2]<-permutation2
logfold<-coexpNetScore1[,2]
logfold[is.na(logfold)]<-0
logfold[(logfold==Inf)]<-max(logfold[!(logfold==Inf)])
logfold[(logfold==(-Inf))]<-min(logfold[!(logfold==(-Inf))])
resNodeW<-RandomWalk2igraph(g,log10(1+log10(1+abs(logfold))),EdgeWeight=FALSE)
lncRNA_resNodeW<-resNodeW[index1_final]
lncRNA_per<-c(lncRNA_per,lncRNA_resNodeW)
}

#######compute p-value and fdr  #########
lncScore<-lncScorenew
pvalue<-as.double()
for(i in 1:dim(lncScore)[1]){  #compute P-value
total_number<-length(lncRNA_per)
high_number<-length(lncRNA_per[lncRNA_per>lncScore[i,3]])
pvalue[i]<-high_number/total_number
}
fdr.est<-function(p) #fdr funtion
{
    m <- length(ind <- which(!is.na(p)))
    fdr <- rep(NA, length(p))
    stat <- cbind(1:length(p), p, fdr)
    stat[ind, 3] <- unlist(lapply(stat[ind, 2], function(x) {
        c <- length(which(stat[ind, 2] <= x))
        m * x/c
    }))
    stat[ind, ] <- stat[ind, ][order(stat[, 2], decreasing = TRUE), 
        ]
    stat[ind, 3] <- cummin(stat[ind, 3])
    fdr <- stat[order(stat[, 1]), 3]
    fdr[which(fdr > 1)] <- 1
    return(fdr)
}
fdr<-fdr.est(pvalue) #compute fdr
lncScore_with_pvalue<-cbind(lncScore,pvalue,fdr)
fc_rank<-length(lncScore_with_pvalue[,2])+1-rank(abs(lncScore_with_pvalue[,2]),ties.method ="max")
URWScore_rank<-length(lncScore_with_pvalue[,3])+1-rank(lncScore_with_pvalue[,3],ties.method ="max")
lncScore_with_pvalue<-cbind(lncScore,pvalue,fdr,fc_rank,URWScore_rank)
lncScore_with_pvalue[,3]<-lncScore_with_pvalue[,3]*1000
colnames(lncScore_with_pvalue)<-c("LncRNA ID","Logfold","URWScore","pvalue","FDR","Rank(Logfold)","Rank(URWScore)")
write.table(lncScore_with_pvalue,file="./lncScore_with_pvalue.txt",row.names=F,col.names=T,sep="\t",quote=F)

