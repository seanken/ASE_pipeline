# Allele Specific Expression Pipeline For 10X Fata

This pipeline is built as part of the ASAP project to take single nuclei data and a phased vcf file and, using that data, get allele specific expression information for the single nuclei data using STARSolo. It requires Java, nextflow, and Conda. Detailed running instructions will be added soon.


# Set Up

You need to have Java installed (v8 or later) and some version of Conda, as well as nextflow (https://www.nextflow.io/). To install this, download it from git hub with:
```
git clone https://github.com/seanken/ASE_pipeline.git
cd ASE_pipeline
gunzip ref/whitelist_v3/whitelist.txt.gz
```


You are then ready to run!

# Generating Reference

In order to run this pipeline you need a STAR reference and a matching gtf file. See STAR documentation for details.

# Phased Genotype Data

This pipeline requires phased genotype data in vcf, vcf.gz, bcf, or bcf.gz form. We have used Eagle2 for our phasing with an online portal (will add link).

# Fastq data

This pipeline requires 10X fastq files. Must be formatted as described by 10X on their website, so should be:

```
name_L00*_R1*fastq.gz
name_L00*_R2*fastq.gz
```

Note the pipeline assumes that the name do not contain L00 or R1/R2, and name of read1 and read2 are the same except for the R1/R2 part.

Instead of passing the fastq directly, directory names are passed. Can pass multiple directories (seperated by commas). All the fastq files in these directories should be from the same sample (so 10X channel and individual).

# How To Run Pipeline

```
nextflow=/path/to/nextflow
pipeline=/path/to/QuantPipeline.nf

$nextflow $pipeline [options]
```

The pipeline has a mix of optional and required options.


