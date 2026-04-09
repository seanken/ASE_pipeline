library(Seurat)
library(tidyr)
library(dplyr)
library(Matrix)
library(glmGamPoi)
library(Azimuth)


##
##dir is the input filtered cell directory, prefix is the prefix used for saving
MakeSeurat<-function(dir,prefix)
{
print("Load data")
dat=ReadInSTARSolo(dir)

print("Make Seurat Object")
minGenes=300
seur<-CreateSeuratObject(dat,"Seurat",min.features=minGenes)
seur<-NormalizeData(seur,normalization.method="LogNormalize",scale.factor=1000000)
seur<-FindVariableFeatures(seur)
seur<-ScaleData(seur,features=seur@assays$RNA@var.features)
seur<-RunPCA(seur,npcs=60)
seur<-RunUMAP(seur,dims=1:20,return.model=T)
saveRDS(seur,paste(prefix,".seur.RDS",sep=""))
write(seur@assays$RNA@var.features,paste(prefix,".var.genes.txt",sep=""))

mn=rowSums(dat)
saveRDS(mn,paste(prefix,".PseudoBulk.RDS",sep=""))

print("Run Azimuth")
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_motorcortex")
meta=RunAzimuth(seur,reference$map)
write.table(meta,paste(prefix,".azimuth.cortex.txt",sep=""),row.names=F,sep="\t",quote=F)
}


##
##Runs Azimuth on a seurat object, seur, using the reference Azimuth object (created with GenerateReference from a Seurat object)
##Assumes the reference of interest is in reference_map@meta.data in the subclass column
##
RunAzimuth<-function(seur,reference_map)
{

print("Apply SCTransform")
seur<-SCTransform(seur,new.assay.name = "refAssay",residual.features = rownames(x = reference_map),reference.SCT.model = reference_map[["refAssay"]]@SCTModel.list$refmodel,method = 'glmGamPoi',ncells = 2000,n_genes = 2000,do.correct.umi = FALSE,do.scale = FALSE,do.center = TRUE)

print("Get anchors")
anchors <- FindTransferAnchors(reference = reference_map,query=seur,k.filter = NA,reference.neighbors = "refdr.annoy.neighbors",reference.assay = "refAssay",query.assay = "refAssay",reference.reduction = "refDR",normalization.method = "SCT",features = intersect(rownames(x = reference_map), VariableFeatures(object = seur)),dims = 1:50,n.trees = 20,mapping.score.k = 100)

print("Move over annotation!")
refdata <- lapply(X = "subclass", function(x) {
  reference_map[[x, drop = TRUE]]
  })
  names(x = refdata) <- "subclass"
  if (FALSE) {
   refdata[["impADT"]] <- GetAssayData(
    object = reference_map[['ADT']],
     slot = 'data')
                   }

seur <- TransferData(reference = reference_map,query=seur,dims = 1:50,anchorset = anchors,refdata = refdata,n.trees = 20,store.weights = TRUE)

return(seur@meta.data)
}




##Read in STARSolo data, given as list of directories (dirs), giving samples names from names of dirs if there, otherwise based off order
##Can run empty dropleta utility if runFilter=T
ReadInSTARSolo<-function(dirs,runFilts=F,numCells=10000)
{
nams=sub("^","samp",1:length(dirs))
if(length(dirs)==length(names(dirs)))
{
nams=names(dirs)
}

print("Read in each matrix")
mats=lapply(1:length(dirs),function(x){
print("Load!")
print(x)
nam=nams[x]
print(nam)
dir=dirs[x]

mat=readMM(paste(dir,"/matrix.mtx",sep=""))
colnames(mat)=sub("^",paste(nam,"_",sep=""),scan(paste(dir,"/barcodes.tsv",sep=""),""))
tab=read.table(paste(dir,"/features.tsv",sep=""),stringsAsFactors=F,sep="\t")
tab["Name"]=tab[,2]
tab[duplicated(tab[,2]),"Name"]=tab[duplicated(tab[,2]),1]
rownames(mat)=tab[,"Name"]

return(mat)
})

print("Combine Matrices")
mat=do.call(cbind,mats)
rm(mats)
print("Done!")
return(mat)

}

if(!interactive())
{
args=commandArgs(trailingOnly=TRUE)

dir=args[1]

prefix=args[2]

MakeSeurat(dir,prefix)
}

