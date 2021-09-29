library(Seurat)
library(dplyr);
library(tidyr)


##
##Assumes has Seurat object with columns corresponding to the sample (samp), the location of the location of SNPLevelCounts/comb.bed output by the pipeline (snp_col), and the location of the AlleleCounts/counts.txt output by the pipeline (all_col). Assumes the cell names in Seurat are [samp name]_[CBC]
##genes is the list of genes to get the SNP info for
##Returns a list with 3 entries: ref, alt, and meta.
GetSNPs<-function(seur,genes,samp="orig.ident",snp_col="SNP",all_col="Allele",cond="PD")
{
print("Make Nice")
tab=seur@meta.data[,c(samp,snp_col,all_col,cond)]
colnames(tab)=c("Samp","SNP","Allele","Condition")
tab["Name"]=rownames(seur@meta.data)
dat<-tab %>% group_by(Samp,SNP,Allele,Condition) %>% summarise() %>% as.data.frame()

print("Get SNP data")
out=lapply(1:dim(dat)[1],function(i){
print(i)
mat=loadSNPs(dat[i,1],dat[i,2],dat[i,3],tab[tab[,"Samp"]==dat[i,1] & tab[,"SNP"]==dat[i,2] & tab[,"Allele"]==dat[i,3],"Name"],genes)
mat["Condition"]=dat[i,"Condition"]
print(head(mat))
return(mat)
})

dat=do.call(rbind,out)

dat=dat[!is.na(dat[,"Condition"]),]

return(dat)

}


##
##Gets pseudobulk allele specific and SNP specific counts
##Takes in:
##Sample name (samp)
##snp comb.bed file location (snp)
##Allele gene level fount data (allele)
##list of cell names to include (nams)
##list of genes to use (genes)
##
##Returns Psuedobulk
loadSNPs=function(samp,snp,allele,nams,genes,getVals=F)
{
print("Load Allele Data")
cnts=read.table(allele,stringsAsFactors=F)
colnames(cnts)=c("CBC","Gene","Allele","Count")
cnts=cnts[cnts$Allele!="Ambig",]
cnts["Name"]=sub("^",paste(samp,"_",sep=""),cnts[,"CBC"])
print(head(cnts))
print(dim(cnts))
cnts=cnts[cnts$Name %in% nams,]
print(dim(cnts))
print(head(cnts))
print("Make PseudoBulk")
cnts<-cnts %>% group_by(Gene,Allele) %>% summarise(Count=sum(Count)) %>% as.data.frame()
cnts["Sample"]=samp
print(head(cnts))
print("Load SNP info")
bed=read.table(snp,stringsAsFactors=F)
colnames(bed)=c("SNP","Gene","Geno")
bed=bed[bed$Gene %in% genes,]
print(head(bed))
print(length(unique(bed$Gene)))
print(length(unique(cnts$Gene)))
print(length(intersect(unique(bed$Gene),unique(cnts$Gene))))
if(getVals)
{
lst=list()
lst[[1]]=cnts
lst[[2]]=bed
return(lst)
}
print("Combine with Allele data")
cnts=inner_join(cnts,bed)
print(length(unique(cnts$Gene)))
lst=c("alt","ref","ref","alt")
names(lst)=c("All1_1|0","All1_0|1","All2_1|0","All2_0|1")
cnts<-cnts %>% unite(Comb,Allele,Geno,sep="_")
cnts=cnts[cnts$Comb %in% names(lst),]
cnts["Allele"]=lst[cnts[,"Comb"]]
cnts=cnts[,c("Gene","Sample","SNP","Allele","Count")]
print(length(unique(cnts$Gene)))
print("Return!")
return(cnts)
}
