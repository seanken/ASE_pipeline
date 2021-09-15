library(dplyr);library(tidyr)

args = commandArgs(trailingOnly=TRUE)

print("Load data")
tab=read.table(args[2],stringsAsFactors=F)
dat=read.table(args[1],stringsAsFactors=F)

colnames(tab)=c("CBC","Gene","Allele","Count")
colnames(dat)=c("SNP","Gene","Geno")

print("Join Gene and SNP level")
comb=inner_join(tab,dat)

print("Clean up")
comb=comb[comb[,"Allele"]!="Ambig",]
comb<-comb %>% unite(Hap,Allele,Geno)
lst=c("ref","alt","alt","ref")
names(lst)=c("All1_0|1","All1_1|0","All2_0|1","All2_1|0")
comb=comb[comb[,"Hap"] %in% names(lst),]
comb["Allele"]=lst[comb[,"Hap"]]
comb<-comb[,c("CBC","SNP","Gene","Allele","Count")] %>% unite(Feature,Gene,SNP,sep="_")

print("Save")
write.table(comb,args[3],sep="\t",quote=F,row.names=F)

