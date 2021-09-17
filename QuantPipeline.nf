#!/usr/bin/env nextflow

//As written runs for one sample at a time, will modify in later versions
//Requirements:
//1) Conda
//2) Java


//Input files Required input
params.input_vcf //vcf used
params.vcf_col // sample name in the vcf to use
params.input_dirs //comma seperated list of directories containing 10X fastqs, will use all fastqs in those directories
params.outdir //directory to output results to


//Params for STARSolo
params.numCells=3000 //number of cells expected

//Files used in process and other params
params.numThreads=8
params.ref="$projectDir/ref/STAR_ref" //STAR reference to use--for now just human built in, though users can pass their own
params.num10X_version="v3" //version of 10X, allows you to pick the correct list
params.whitelist="$projectDir/ref/whitelist_${params.num10X_version}/whitelist.txt" //10X whitelist
params.use_conda=0 //decide if want to download dependencies with conda. If not must be on path
params.gtf="$projectDir/ref/genes.gtf" //gtf used, include standard one but allow user to overright, should match STAR reference
params.MakeSNPScript="$projectDir/scripts/Make.SNPs.R" //R script used for making SNP counts, probably don't want user to change in most situations
params.AlleleCountJar="$projectDir/scripts/AlleleCount.jar"
params.UMILen=12
params.QCScript="scripts/Phased.UMI.QC.R"
params.makeBeds="$projectDir/scripts/makeBeds.sh"

//this step processes the VCF into the form you need it (assumes already phased)
process PrepVCF
{

input:
env vcf_col from params.vcf_col
path input_vcf, stageAs:"input.vcf.gz" from params.input_vcf

output:
path "new.vcf" into new_vcf_ch



conda 'bcftools'

'''
bcftools view -H -O v -s $vcf_col input.vcf.gz | grep -v "0|0" | grep -v "1|1" > new.vcf
'''


//bcftools view -H -O v -s $vcf_col input.vcf.gz > new.vcf

}

//Gets path to Fastqs
process GetFastqPath
{
input:
env directs from params.input_dirs

output:
env fastqs into fastq_ch


'''
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

publishDir "${params.outdir}/STARSolo"

input:
env fastqs from fastq_ch
path vcf, stageAs:"new.vcf" from new_vcf_ch
env numCells from params.numCells
env ref from params.ref
env fastqs from fastq_ch
path whitelist, stageAs: "whitelist.txt" from params.whitelist
env numThreads from params.numThreads
env UMILen from params.UMILen

output:
path "output" into STAR_Dir

//if(params.use_conda==1)
//{
conda 'star=2.7.9a'
//}

'''
mkdir output
STAR --genomeDir $ref --readFilesIn $fastqs --soloType CB_UMI_Simple --soloCBwhitelist  whitelist.txt --soloUMIlen $UMILen --soloUMIfiltering MultiGeneUMI --soloCBmatchWLtype 1MM_multi_pseudocounts --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM vG vA vG --outSAMtype BAM SortedByCoordinate --soloCellFilter CellRanger2.2 $numCells 0.99 10 --runThreadN $numThreads --outFileNamePrefix output/results --readFilesCommand zcat --varVCFfile new.vcf --waspOutputMode SAMtag --soloFeatures GeneFull
'''

}



//This process runs a java based method to count number of alleles, giving results as gene level
process CountAlleles
{
publishDir "${params.outdir}/AlleleCounts"

input:
path output, stageAs:"output" from STAR_Dir
path vcf, stageAs:"new.vcf" from new_vcf_ch 
path AlleleCountJar, stageAs:"AlleleCount.jar" from params.AlleleCountJar

output:
path "counts.txt" into gene_counts_ch

'''
java -jar AlleleCount.jar output/resultsAligned.sortedByCoord.out.bam counts.txt
'''

}




//Extends allele counts to SNP level
process SNPLevel
{
publishDir "${params.outdir}/SNPLevelCounts"

input:
path gtf, stageAs:"genes.gtf" from params.gtf
path vcf, stageAs:"new.vcf" from new_vcf_ch
path counts, stageAs:"counts.txt" from gene_counts_ch
path MakeSNPScript, stageAs:"Make.SNPs.R" from params.MakeSNPScript
path makeBeds, stageAs:"makeBeds.sh" from params.makeBeds //Makes gene and SNP beds

output:
path "SNP.counts.txt" into SNPS_allele_ch


conda 'bedtools=2.30.0 r-tidyr=1.1.3'




'''
source makeBeds.sh genes.gtf new.vcf
bedtools intersect -a snps.bed -b gene.bed -wa -wb | awk '{print $4"\t"$9"\t"$5}' | sort | uniq > comb.bed
Rscript Make.SNPs.R comb.bed counts.txt SNP.counts.txt
'''

}


process PhasedUMI_QC
{

publishDir "${params.outdir}/SNPLevelCounts"

input:
path QCScript, stageAs:"Phased.UMI.QC.R" from params.QCScript
path counts, stageAs:"counts.txt" from gene_counts_ch
path output, stageAs:"output" from STAR_Dir

output:
path "hist.ratio.pdf" into QC_pdf
path "Basic.QC.txt" into QC_tot
path "UMI.counts.by.gene.txt" into QC_counts

conda 'r-tidyr=1.1.3 r-matrix'


'''
Rscript Phased.UMI.QC.R counts.txt output
'''



}

