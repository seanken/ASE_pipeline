#! /bin/bash

#$ -cwd
#$ -e ErrFiles/star.err
#$ -o ErrFiles/star.log
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=20g
#$ -l h_rt=120:00:00
#$ -l os="RedHat7"

#$ -pe smp 8
#$ -binding linear:8
#$ -R y
#$ -t 1-12
#$ -tc 20

source /broad/software/scripts/useuse



SEEDFILE=samps.txt

input_dirs=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}')
vcf_col=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $2}')
nam=$vcf_col

workdir=samp_$nam

mkdir $workdir

cd $workdir

nextflow=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version2.0/ASE_pipeline/nextflow
pipeline=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version2.0/ASE_pipeline/QuantPipeline.nf 
vcf=/stanley/levin_asap_storage/612-eqtl/GenotypeData/Clean_vcf/Combined/comb.no.chr.vcf.gz
outdir=output
mkdir output

$nextflow $pipeline --input_vcf $vcf --vcf_col $vcf_col --input_dirs $input_dirs --outdir $outdir

