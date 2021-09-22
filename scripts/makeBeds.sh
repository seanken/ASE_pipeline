gtf=$1
vcf=$2

grep -v \# genes.gtf | awk '{if($3=="gene"){print $1"\t"$4"\t"$5"\t"$14}}' | sed 's/\"//g' | sed 's/;//g' > gene.bed
awk '{print $1"\t"$2-1"\t"$2+1"\t"$3}' new.vcf > snps.v1.bed
awk '{print $10}' new.vcf | awk -F ':' '{print $1}' > snps.v2.bed
paste snps.v1.bed snps.v2.bed > snps.bed
rm snps.v1.bed
rm snps.v2.bed
