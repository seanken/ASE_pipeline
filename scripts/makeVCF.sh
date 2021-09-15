input=$1
vcf_col=$2
output=$3

zcat $input | grep -v \#  | awk -v col=$vcf_col '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$col}' | grep -v 0\|0 | grep -v 1\|1 > $output
