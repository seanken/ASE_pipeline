library(dplyr)
library(tidyr)
library(MASS)
library(tictoc)
library(lme4)

##
##Takes in a table of SNP counts produced by LoadSNPs.R and tests for AI SNPs 
##For COndition: form=Count~Allele+Sample+Condition*Allele,coef_ret="Alleleref:Conditionyes"
##
TestSNP<-function(dat,minCount=50,minSamp=10,form=Count~Allele+Sample,coef_ret="Allele",numTest=-1)
{
print("Format")
dat<-dat %>% unite(Feature,Gene,SNP,sep="_",remove=F)
tab<-dat %>% group_by(Feature,Sample) %>% summarise(Num=length(unique(Allele)),Count=sum(Count)) %>% group_by(Feature) %>% summarise(Count=sum(Count),NumSamp=sum(Num>1)) %>% as.data.frame()
print(length(unique(dat$Gene)))
feats=tab[tab$NumSamp>minSamp & tab$Count>minCount,"Feature"]
dat=dat[dat$Feature %in% feats,]
print("Number to Test:")
print(length(feats))
print(length(unique(dat$Gene)))



print("Split by SNP")
#return(dat)########
#bySNP=lapply(feats,function(x){dat[dat[,1]==x,]})
#names(bySNP)=feats
bySNP=split(dat,f=dat$Feature)
if(numTest>0){bySNP=bySNP[1:numTest]}
print("Test")
tic()
nams=names(bySNP)
out=lapply(names(bySNP),function(cur_feat){
x=bySNP[[cur_feat]]
#tab=x[,c("Sample","Allele","Count")] %>% spread(Sample,Count,fill=0) %>% gather(Sample,Count,-Allele) %>% unite(Nam,Sample,Allele,sep="_") %>% spread(Nam,Count)
#meta=data.frame(Nam=colnames(tab))
#meta["Sample"]=
#fit=glm.nb(Count~Allele+Sample+Condition*Allele,x)
fit=tryCatch({glm.nb(form,x)},error=function(cond){return(NULL)})
#fit=tryCatch({glm.nb(Count~Allele+Condition*Allele+Sample,x)},error=function(cond){print("Error!");return(NULL)})
if(is.null(fit))
{
return(NULL)
}
theta=fit$theta
warned=F
if(!is.null(fit$th.warn)){warned=T}
#if(!is.null(fit$th.warn)){warned=T;fit=tryCatch({glm(form,x,family=poisson)},error=function(cond){return(NULL)})}
#if(is.null(fit))
#{
#return(NULL)
#}

#if(!is.null(fit$th.warn)){fit=glmer(Count~offset(log(tot))+Allele+(1|sample),x,family=poisson)}
coef=summary(fit)$coefficients
coef=data.frame(coef)
#print(rownames(coef))
coef["Test"]=rownames(coef)
coef=coef[grep(paste("^",coef_ret,sep=""),rownames(coef)),]
coef["SNP"]=cur_feat
coef["NumSamp"]=length(unique(x$Sample))
coef["Theta"]=theta
coef["Count"]=sum(x[,"Count"])
coef["Warn"]=warned
return(coef)
#return(fit)

})
toc()
#names(bySNP)=feats[1:10]
bySNP=out
names(bySNP)=nams
bySNP[sapply(bySNP, is.null)]=NULL

ret=do.call(rbind,bySNP)

ret=data.frame(ret)
#ret["SNP"]=names(bySNP)

colnames(ret)[4]="pval"
#print(head(ret))

ret=ret[order(ret$pval),]
ret["padj"]=p.adjust(ret[,"pval"])

ret["Gene"]=as.character(lapply(ret[,"SNP"],function(x){strsplit(x,split="_")[[1]][1]}))
rownames(ret)=NULL
return(ret)

}
