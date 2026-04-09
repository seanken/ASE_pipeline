library(dplyr);library(tidyr);library(tidytext);library(Matrix)

##
##Changes gene level results into SNP level results, takes 5 arguments
##1) The SNP, gene, Geno info from the vcf
##2) The Gene level allele counts
##3) A file with a list of SNPs
##4) A file with a list of cells
##5) The prefix for the files to save to
##

args = commandArgs(trailingOnly=TRUE)

print("Load data")
tab=read.table(args[2],stringsAsFactors=F)
dat=read.table(args[1],stringsAsFactors=F)
snps=scan(args[3],"")
cbcs=scan(args[4],"")
colnames(tab)=c("CBC","Gene","Allele","Count")
colnames(dat)=c("SNP","Gene","Geno")

print(head(dat))
print(head(snps))

#snps=gsub("chr","",gsub(":","_",snps))

print(head(snps))

print(dim(dat))
print(dim(tab))

tab=tab[tab[,1] %in% cbcs,]
dat=dat[dat[,1] %in% snps,]


print(dim(tab))
print(dim(dat))

print("Join Gene and SNP level")
comb=inner_join(tab,dat)

print(object.size(comb))

print("Clean up")
comb=comb[comb[,"Allele"]!="Ambig",]
comb<-comb %>% unite(Hap,Allele,Geno)
lst=c("ref","alt","alt","ref")
names(lst)=c("All1_0|1","All1_1|0","All2_0|1","All2_1|0")
comb=comb[comb[,"Hap"] %in% names(lst),]
comb["Allele"]=lst[comb[,"Hap"]]
comb<-comb[,c("CBC","SNP","Gene","Allele","Count")] %>% unite(Feature,Gene,SNP,Allele,sep="_")
print("Make Sparse Matrix")
#comb_ref=comb[comb$Allele=="ref",]
#comb_alt=comb[comb$Allele=="alt",]
comb<-comb %>% cast_sparse(Feature,CBC,Count)

print(object.size(comb))

print("Save")
filnam=args[5]
writeMM(comb,paste(filnam,".matrix.mtx",sep=""))
write(rownames(comb),paste(filnam,".snps.txt",sep=""))
write(colnames(comb),paste(filnam,".cbc.txt",sep=""))


