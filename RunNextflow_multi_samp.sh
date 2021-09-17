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


#export CLASSPATH=/ahg/regevdata/projects/612-eqtl/eQTL_testing/STARSolo/Code:/ahg/regevdata/projects/612-eqtl/eQTL_testing/Code_for_Kiran/Phasing/HTSJDK/htsjdk/build/libs/htsjdk-2.21.1-3-g3df5a35-SNAPSHOT.jar

SEEDFILE=samps.txt

input_dirs=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}')
vcf_col=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $2}')
nam=$vcf_col

workdir=samp_$nam

mkdir $workdir

cd $workdir

nextflow=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version1.0/nextflow
pipeline=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version1.0/QuantPipeline.nf 
vcf=/stanley/levin_asap_storage/612-eqtl/GenotypeData/Clean_vcf/Combined/comb.no.chr.vcf.gz
outdir=output
mkdir output




$nextflow $pipeline --input_vcf $vcf --vcf_col $vcf_col --input_dirs $input_dirs --outdir $outdir --use_cond 1 -resume

