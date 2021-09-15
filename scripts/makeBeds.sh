gtf=$1
vcf=$2

grep -v \# genes.gtf | awk '{if($3=="gene"){print $1"\t"$4"\t"$5"\t"$14}}' | sed 's/\"//g' | sed 's/;//g' > gene.bed
awk '{print $1"\t"$2-1"\t"$2+1"\t"$1"_"$2"_"$4"_"$5"\t"$10}' new.vcf | awk -F ':' '{print $1}' > snps.bed
