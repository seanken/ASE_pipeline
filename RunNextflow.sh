export CLASSPATH=/ahg/regevdata/projects/612-eqtl/eQTL_testing/STARSolo/Code:/ahg/regevdata/projects/612-eqtl/eQTL_testing/Code_for_Kiran/Phasing/HTSJDK/htsjdk/build/libs/htsjdk-2.21.1-3-g3df5a35-SNAPSHOT.jar

nextflow=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version1.0/nextflow
pipeline=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version1.0/QuantPipeline.nf 
vcf=/ahg/regevdata/projects/612-eqtl/eQTL_testing/OutsideData/ASD/Code/WES/Phasing/Eagle2/Test/comb.vcf.gz
vcf_col=10
input_dirs=/stanley/levin_asap_storage/612-eqtl/eQTL_testing/OutsideData/ASD/Data/SingleCell/reads/samp_5278_PFC_10x/
input_dirs=/broad/hptmp/ssimmons/QTL/TestPipeline/TestFA/
outdir=output
mkdir output

$nextflow $pipeline --input_vcf $vcf --vcf_col $vcf_col --input_dirs $input_dirs --outdir $outdir --use_cond 1 -resume

