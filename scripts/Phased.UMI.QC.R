library(dplyr);
library(tidyr);
library(Matrix);


getQC=function(count_fil,starsolo_dir,quantType="GeneFull")
{
print("Get Allele Counts")
allele=read.table(count_fil,stringsAsFactors=F)
print(head(allele))
colnames(allele)=c("CBC","Gene","Allele","nUMI")
allele=allele[grep("All",allele[,"Allele"]),]

print("Get Gene Counts")
dat=readMM(paste(starsolo_dir,"/resultsSolo.out/",quantType,"/filtered/matrix.mtx",sep=""))
genes=read.table(paste(starsolo_dir,"/resultsSolo.out/",quantType,"/filtered/features.tsv",sep=""),stringsAsFactors=F,sep="\t")
colnames(genes)=c("Gene_id","Gene","Type")
genes["TotUMI"]=rowSums(dat)
genes<-genes %>% group_by(Gene) %>% summarise(Tot=sum(TotUMI)) %>% as.data.frame()
cells=scan(paste(starsolo_dir,"/resultsSolo.out/",quantType,"/filtered/barcodes.tsv",sep=""),"")

allele=allele[allele$CBC %in% cells,] %>% group_by(Gene,Allele) %>% summarise(nUMI=sum(nUMI)) %>% spread(Allele,nUMI,fill=0) %>% as.data.frame()

print("Make QC!")
dat=inner_join(allele,genes)
write.table(dat,"UMI.counts.by.gene.txt",sep="\t",quote=F,row.names=F)
dat=dat[dat$Tot>0,]

numPhased=sum(dat[,"All1"])+sum(dat[,"All2"])
percentPhased=numPhased/sum(dat[,"Tot"])
numPhasedPerCell=numPhased/length(cells)
numGenes10=sum(dat[,"All1"]+dat[,"All2"]>10)

dat=dat[dat$All1+dat$All2>10,]
dat["Ratio"]=dat[,"All1"]/(dat[,"All1"]+dat[,"All2"])

meanRatio=mean(dat[,"Ratio"])
medRatio=median(dat[,"Ratio"])

QC=data.frame(QC_Name=c("Num Phased UMI","Percent UMI Phased","Number Phased UMI per Cell","Number genes with >10 phased UMIs","Mean AI","Median AI"),QC_value=c(numPhased,percentPhased,numPhasedPerCell,numGenes10,meanRatio,medRatio))

write.table(QC,"Basic.QC.txt",sep="\t",quote=F,row.names=F)

pdf("hist.ratio.pdf")
hist(dat[,"Ratio"],100)
dev.off()

}

if(!interactive())
{
args = commandArgs(trailingOnly=TRUE)
count_fil=args[1]
starsolo_dir=args[2]
quantType=args[3]
getQC(count_fil,starsolo_dir,quantType)
}



