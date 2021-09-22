from pysam import VariantFile

##
##Reads in a vcf than keeps the SNPs with at least numHet hetero sites
##Writes the results to outfil
##
def GetUsefulVCF(vcf_file,numHet,outfil):
	print("Start reading in")
	vcf=VariantFile(vcf_file)
	fil=open(outfil,"w")
	i=0
	for rec in vcf.fetch():
		i=i+1
		if i % 10000==0:
			print(i);
		geno = [s['GT'] for s in rec.samples.values()] ##Get genotypes, taken from github help page
		#print(geno[1:10])
		#print(rec.id)
		numHet_cur=len([l for l in geno if sum(l)==1])
		#print(numHet_cur)
		if numHet_cur<numHet:
			continue;
		fil.write(rec.id+"\n")
	vcf.close()
	fil.close()
		
