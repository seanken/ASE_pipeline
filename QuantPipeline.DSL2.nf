#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//As written runs for one sample at a time, will modify in later versions
//Requirements:
//1) Conda
//2) Java


//Input files Required input
params.input_vcf=null //vcf used
params.vcf_col="" // sample name in the vcf to use
params.input_dirs=null //comma seperated list of directories containing 10X fastqs, will use all fastqs in those directories, can also be prefixes
params.input_bam=null //If input bam instead of fastq
params.input_bam_bai="${params.input_bam}.bai"
params.outdir=null //directory to output results to
//params.snps //list of SNPs to use
params.cellFile=null //list of cells to use, optional


//Params for STARSolo
params.numCells=3000 //number of cells expected

//Files used in process and other params
params.numThreads=8
params.ref="$projectDir/ref/STAR_ref" //STAR reference to use--for now just human built in, though users can pass their own
params.num10X_version="v3" //version of 10X, allows you to pick the correct list
params.whitelist="$projectDir/ref/whitelist_${params.num10X_version}/whitelist.txt" //10X whitelist
params.gtf="$projectDir/ref/genes.gtf" //gtf used, include standard one but allow user to overright, should match STAR reference
//params.MakeSNPScript="$projectDir/scripts/Make.SNPs.R" //R script used for making SNP counts, probably don't want user to change in most situations
AlleleCountJar="$projectDir/scripts/AlleleCount.jar"
params.UMILen=12
QCScript="$projectDir/scripts/Phased.UMI.QC.R"
makeBeds="$projectDir/scripts/makeBeds.sh"
params.featSTARSolo="GeneFull"


//Check inputs
if(params.input_dirs && params.input_bam)
{
    println("Can only pass one of input_dirs or input_bam")
    System.exit(0)
}

if(!params.input_dirs && !params.input_bam)
{
    println("Need to pass one of input_dirs or input_bam")
    System.exit(0)
}

if(!params.input_vcf || !params.vcf_col)
{
    println("Need to pass input_vcf and vcf_col")
    System.exit(0)
}

if(!params.outdir)
{
    println("Need to pass output directory (outdir)")
    System.exit(0)
}

if(params.cellFile && !params.input_bam)
{
    println("If cell file need bam file")
    System.exit(0)
}


println("Run pipeline")

workflow{
    //Start by prepping the data
    newvcf=PrepVCF(params.vcf_col,params.input_vcf)
    //bamFile=params.input_bam
    
 
    //println(inputDirs)


    if(params.cellFile) //If want to select certain cells does so, updates bam file
    {
        GetCells(params.input_bam,params.input_bam_bai,params.cellFile,params.vcf_col,params.numThreads)
    }

    bamFile= params.cellFile ? GetCells.out : params.input_bam
    //println(bamFile)

    //if(!params.cellFile && params.input_bam) //If input is bam file processes to fastq, updates inputDir
    //{
    //    ProcessInputBam(params.input_bam)
    //}
    //if(params.cellFile && params.input_bam)
    //{
    //    ProcessInputBam(GetCells.out)
    //}
    if(params.input_bam)
    {
        ProcessInputBam(bamFile)
    }

    inputDirs=params.input_dirs ? params.input_dirs : ProcessInputBam.out

    faPath=GetFastqPath(inputDirs)

    //Run STARSolo
    starSoloDir=RunSTARSolo(faPath,newvcf,params.numCells,params.ref,params.whitelist,params.numThreads,params.UMILen,params.featSTARSolo)
    
    //Get ASE information
    countsGeneLevel=CountAlleles(starSoloDir,AlleleCountJar)
    SNPLevel(params.gtf,newvcf,makeBeds)
    
    //Get QC information
    PhasedUMI_QC(QCScript,countsGeneLevel,starSoloDir,params.featSTARSolo)
}







//Gets cells of interest, only if multiple samples mixed or only want to run on subsample
//Not always Run
process GetCells
{
    input:
    path "input.bam"
    path "input.bam.bai"
    path "cells.txt"
    env vcf_col 
    env numThreads 

    output:
    path "output.bam"

    '''
    grep ${vcf_col}\$ cells.txt | awk '{print \$1}' | sed 's/$/\toutput/g' > cells.sinto.txt
    sinto filterbarcodes -b input.bam -c cells.sinto.txt -p $numThreads 
    '''

    //sed 's/$/\toutput/g' cells.txt > cells.sinto.txt
}



//Take a bam file from CellRanger and makes into fastq files
//only run if using bam instead of fastq
process ProcessInputBam
{

    //publishDir "${params.outdir}", mode: 'rellink'

    input:
    path "input.bam"
    
    
    output:
    env pathRet

    '''
    cellranger bamtofastq input.bam FA_DS
    pathRet=$(ls -d $PWD/FA_DS/s*)/
    '''

}


//this step processes the VCF into the form you need it (assumes already phased)
process PrepVCF
{

    input:
    env vcf_col 
    path "input.vcf.gz" 

    output:
    path "new.vcf" 


    '''
    tabix -p vcf input.vcf.gz
    bcftools view -H -O v -s $vcf_col input.vcf.gz | grep -v "0|0" | grep -v "1|1" | grep -v "\\.|\\."  > new.vcf
    '''


}

//Gets path to Fastqs
process GetFastqPath
{
    input:
    env directs 

    output:
    env fastqs 


    '''
    echo $directs
    toSearch=$(echo $directs | sed 's/,/*_L00*_R2*fastq.gz /g'| sed 's/$/*_L00*_R2*fastq.gz/g')
    fastq1=$(ls -m \$toSearch | tr -d '[:space:]')

    toSearch=$(echo $directs | sed 's/,/*_L00*_R1*fastq.gz /g'| sed 's/$/*_L00*_R1*fastq.gz/g')
    fastq2=$(ls -m \$toSearch | tr -d '[:space:]')

    fastqs=$(echo \$fastq1 \$fastq2)
    '''

}


//Runs StarSolo with the inputs given, giving expression info and bams that will be used to calculate allelic expression later
process RunSTARSolo
{

    publishDir "${params.outdir}/STARSolo", mode: 'rellink'

    input:
    env fastqs 
    path "new.vcf"
    env numCells
    env ref 
    path "whitelist.txt" 
    env numThreads 
    env UMILen 
    env feat 

    output:
    path "output" 


    '''
    mkdir output
    echo STAR --genomeDir $ref --readFilesIn $fastqs --soloType CB_UMI_Simple --soloCBwhitelist  whitelist.txt --soloUMIlen $UMILen --soloUMIfiltering MultiGeneUMI_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM vG vA vG --outSAMtype BAM SortedByCoordinate --soloCellFilter EmptyDrops_CR $numCells 0.99 10 45000 90000 500 0.01 20000 0.01 10000 --runThreadN $numThreads --outFileNamePrefix output/results --readFilesCommand zcat --varVCFfile new.vcf --waspOutputMode SAMtag --soloFeatures $feat --limitOutSJcollapsed 10000000 --limitIObufferSize=1500000000 --limitBAMsortRAM 60000000000 --outFilterScoreMin 30 --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 
    STAR --genomeDir $ref --readFilesIn $fastqs --soloType CB_UMI_Simple --soloCBwhitelist  whitelist.txt --soloUMIlen $UMILen --soloUMIfiltering MultiGeneUMI_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM vG vA vG --outSAMtype BAM SortedByCoordinate --soloCellFilter EmptyDrops_CR $numCells 0.99 10 45000 90000 500 0.01 20000 0.01 10000 --runThreadN $numThreads --outFileNamePrefix output/results --readFilesCommand zcat --varVCFfile new.vcf --waspOutputMode SAMtag --soloFeatures $feat --limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 1500000000 --limitBAMsortRAM 60000000000 --outFilterScoreMin 30 --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 
    '''

}



//This process runs a java based method to count number of alleles, giving results as gene level
process CountAlleles
{
    publishDir "${params.outdir}/AlleleCounts", mode: 'rellink'

    input:
    path "output" 
    path "AlleleCount.jar" 

    output:
    path "counts.txt" 

    '''
    java -jar AlleleCount.jar output/resultsAligned.sortedByCoord.out.bam counts.txt
    '''

}




//Extends allele counts to SNP level--creates bed file that will be used later fo create SNP matrix (otherwise matrix too big)
process SNPLevel
{
    publishDir "${params.outdir}/SNPLevelCounts", mode: 'rellink'

    input:
    path "genes.gtf"
    path "new.vcf"
    path "makeBeds.sh" //Makes gene and SNP beds

    output:
    path "comb.bed" 
    path "snps.bed" 




    '''
    source makeBeds.sh genes.gtf new.vcf
    bedtools intersect -a snps.bed -b gene.bed -wa -wb | awk '{print $4"\t"$9"\t"$5}' | sort | uniq > comb.bed
    '''


}


//Calculates some basic QC about the phased data
process PhasedUMI_QC
{

    publishDir "${params.outdir}/QCPhase", mode: 'move'

    input:
    path "Phased.UMI.QC.R" 
    path "counts.txt" 
    path "output" 
    env feat 

    output:
    path "hist.ratio.pdf" 
    path "Basic.QC.txt" 
    path "UMI.counts.by.gene.txt"

    //conda 'r-tidyr=1.1.3 r-matrix'


    '''
    Rscript Phased.UMI.QC.R counts.txt output $feat
    '''



}



