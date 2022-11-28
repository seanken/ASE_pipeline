params.bam //the bam file from STARSolo
params.peaks //the peak file from Sierra
params.outdir //the name of the output directory
params.textToBed="$projectDir/scripts/peakToBed.R" //R script to turn peak file to bed file
params.countJar="$projectDir/scripts/AlleleCount.jar" //The Jar file for counting


process peaksToBed //Takes in the peak file and makes a bed file
{
input:
path peaks, stageAs:"peaks.txt" from params.peaks
path textToBed, stageAs:"peakToBed.R" from params.textToBed

output:
path "peaks.sort.bed" into peak_bed

'''
Rscript peakToBed.R peaks.txt peak.bed
bedtools sort -i peak.bed > peaks.sort.bed
'''

}

process AnnotateBam //Use the bed from the peaksToBed step to annotate read overlapping peaks
{

input:
path bam, stageAs:"input.bam" from params.bam 
path bed, stageAs:"peaks.bed" from peak_bed

output:
path "annotated.bam" into annot_bam

'''
bedtools tag -names -s -tag PK -i input.bam -files peaks.bed > annotated.bam
'''

}

process GetCounts //Uses annotated bam to get ASE counts 
{
publishDir "${params.outdir}", mode: 'move'

input:
path bam, stageAs:"annotated.bam" from annot_bam
path countJar, stageAs:"CountSierra.jar" from params.countJar

output:
path "counts.txt" into counts

'''
java -jar CountSierra.jar annotated.bam counts.txt PK
'''

}
