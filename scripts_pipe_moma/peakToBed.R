##
##Turns a peak file from Sierra into a bed file for annotation purposes
##

args=commandArgs(trailingOnly=TRUE)
peaks=args[1]
outfile=args[2]


dat=read.table(peaks,header=T)


out=apply(dat,1,function(x){
chr=x["Chr"]
nam=x["polyA_ID"]
strand=x["Strand"]
if(strand=="-1"){strand="-"}
if(strand!="-"){strand="+"}


start=x["Fit.start"]
end=x["Fit.end"]

#if(strand=="-")
#{
#start=x["Fit.end"]
#end=x["Fit.start"]
#}

if(x["exon.intron"]=="no-junctions")
{
tab=data.frame(Chr=chr,start=start,end=end,name=nam,score=".",strand=strand)
return(tab)
}

val=x["exon.pos"]
lst=strsplit(gsub("(","",val,fixed=T),")",fixed=T)[[1]]

num=length(lst)
start=c(as.character(lapply(lst,function(y){strsplit(y,",")[[1]][1]})))
end=c(as.character(lapply(lst,function(y){strsplit(y,",")[[1]][2]})))



tab=data.frame(Chr=chr,start=start,end=end,name=nam,score=".",strand=strand)

})

bed=do.call(rbind,out)


bed=bed[bed[,2]!=bed[,3],]

bed[2]=gsub(" ","",bed[,2])
bed[3]=gsub(" ","",bed[,3])

#start=bed[,"start"]
#end=bed[,"end"]


#bed[dat[,6]=="-",2]=end[dat[,6]=="-"]
#bed[dat[,6]=="-",3]=start[dat[,6]=="-"]
bed["name"]=gsub("-","_",gsub(":","_",bed[,"name"],fixed=T),fixed=T)

write.table(bed,outfile,sep="\t",quote=F,row.names=F,col.names=F)
