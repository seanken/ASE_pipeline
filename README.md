# Allele Specific Expression Pipeline For 10X Data

This pipeline is built as part of the ASAP project to take single nuclei data and a phased vcf file and, using that data, get allele specific expression information for the single nuclei data using STARSolo. It requires Java, nextflow, and Conda. Detailed running instructions will be added soon.


## Set Up

You need to have Java installed (v8 or later) and some version of Conda, as well as nextflow (https://www.nextflow.io/). To install this, download it from git hub with:
```
git clone https://github.com/seanken/ASE_pipeline.git
cd ASE_pipeline
gunzip ref/whitelist_v3/whitelist.txt.gz
```


You are then ready to run!

## Requirements

The pipeline is set up so that all you need to run it is Conda, nextflow, and java as mentioned above. It is, however, possible to run the pipeline without conda. In that case, however, the number of requirements greatly increases. In that case, in addition to Java and nextflow, you will need the following installed/on the PATH:
1) bcftools
2) bedtools
3) R with tidyr, dply, and Matrix installed
4) STAR v2.7.8a (note the code does not work with earlier (due to changes in the way GN/GX is calculated for GeneFull) or later versions (due to changes in --limitOutSJcollapsed), though the script can be easily modified to do so)

A bash script to install these with conda is in scripts/makeConda.sh. Though this still requires conda, it allows you to create one conda environment instead of one per sample processed.

If you want to run the annotation step also need to have Seurat v4.0 installed and Azimuth v.0.4.3. See https://satijalab.org/seurat/articles/install.html for details on installing Seurat, for Azimuth see their github.

We are planning on creating a docker container for this pipeline as well, but have not done so yet.



## Generating Reference

In order to run this pipeline you need a STAR reference and a matching gtf file. See STAR documentation for details.

## Phased Genotype Data

This pipeline requires phased genotype data in vcf, vcf.gz, bcf, or bcf.gz form. We have used Eagle2 for our phasing with an online portal (will add link).

## Fastq data

This pipeline requires 10X fastq files. Must be formatted as described by 10X on their website, so should be:

```
name_L00*_R1*fastq.gz
name_L00*_R2*fastq.gz
```

Note the pipeline assumes that the name do not contain L00 or R1/R2, and name of read1 and read2 are the same except for the R1/R2 part.

Instead of passing the fastq directly, directory names are passed. Can pass multiple directories (seperated by commas). All the fastq files in these directories should be from the same sample (so 10X channel and individual).

## How To Run Pipeline

```
nextflow=/path/to/nextflow
pipeline=/path/to/QuantPipeline.nf

$nextflow $pipeline [options]
```

The pipeline has a mix of optional and required options.

## Common Errors

When running the pipeline with conda activated, we have sometimes gotten the error:

```
.command.run: line 92: /bin/activate: No such file or directory
```

for the PrepVCF step. This seems to be because 'conda info --json' throw an error (as does 'conda info --envs'). We found the issue for us is with ~/.conda, and we simply removed that directory. Note take this advide with a grain of salt --it worked for us and how we have our conda directory set up, but might not work in all cases and might cause issues with other conda environments.

