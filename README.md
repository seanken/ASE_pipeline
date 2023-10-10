# Allele Specific Expression Pipeline For 10X Data

This pipeline is built as part of the ASAP project to take single nuclei data and a phased vcf file and, using that data, get allele specific expression information for the single nuclei data using STARSolo. It requires Java, nextflow, and Conda. Detailed running instructions will be added soon.


## Set Up

You need to have Java installed (tested with v11 or later) and some version of Conda (or, if not, install all packages independently), as well as nextflow (https://www.nextflow.io/). To install this, download it from git hub with:
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
4) STAR v2.7.10 (does not work with earlier versions, might work with newer versions but not tested)

A bash script to install these with conda is in scripts/makeConda.star10.sh. 

We are planning on creating a docker and/or singularity container for this pipeline as well, but have not done so yet.


## Generating Reference

In order to run this pipeline you need a STAR reference and a matching gtf file. See STARSolo documentation for details on generating such a reference.

## Phased Genotype Data

This pipeline requires phased genotype data in vcf, vcf.gz, bcf, or bcf.gz form. We have used Eagle2 for our phasing with an online portal, but alternative approaches whould work fine.

## Fastq data

This pipeline requires 10X fastq files. Must be formatted as described by 10X on their website, so should be:

```
name_L00*_R1*fastq.gz
name_L00*_R2*fastq.gz
```

Note the pipeline assumes that the name do not contain L00 or R1/R2, and name of read1 and read2 are the same except for the R1/R2 part.

Instead of passing the fastq directly, directory names and/or prefixes are passed. Can pass multiple directories (seperated by commas). If passing directories will us all fastqs in each directory, the directory names must end in a /. If prefixes will use all fastq of the form $prefix/*. 

## How To Run Pipeline

```
nextflow=/path/to/nextflow
pipeline=/path/to/QuantPipeline.DSL2.nf

$nextflow $pipeline [options]
```

The pipeline has a mix of optional and required options.

## Options

A list of options:

`--input_vcf:` The input vcf, assumed to be compressed (.gz). This is required.

`--vcf_col:` The name of the sample corresponding to this indiviudal in the vcf. This is required.

`--input_dirs:` A comma seperated list of fastq locations from the single cell/nuclei experiment. Can either pass the prefix (in which case will use all fastq of the form $prefix*) or the directory the fastq are in (in which case all fastq in those directories will be considered). Note if passing directories each directory name must end in /. Required unless input_bam is passed.

`--input_bam:` A bam file from a CellRanger piepline. An alternative to the input_dirs option. In general not meant to be used unless demultiplexing samples is needed.

`--input_bam_bai:` If input_bam is used this points to the associated index. Assumes this is the input bam name with .bai added if not specified. Only used if input_bam is passed.

`--outdir:` The name of the output directory. This is required.

`--sortRam:` The amount of ram used by STAR for sorting. Sometimes required if the pipeline crashes.

`--cellFile:` Only used/needed if demultiplexing and if input_bam is used. A tab seperated file, the first column being the cell barcodes, the second being the sample each cell barcode is assigned to. The pipeline will use all cells with a sample name matching the vcf_col argument. Should include the -1 suffix in cell barcodes. 

`--numCells:` The number of expected cells/nuclei. The number can be based on information from the experiment or estimates from other tools (CellRanger, etc). 3000 by default.

`--numThreads:` The number of threads for STAR to use. 8 by default.

`--ref:` The STAR reference to use. By default uses a reference in ref/STAR_ref, though this must be built before use if it is left as the default.

`--num10X_version:` Specifies the 10X version. Not required. Options are v2, v3, and visium, v3 by default. Only used for IDing the correct whitelist.

`--whitelist:` List of possible cell/nuclei barcodes. By default will use the built in list for 10X, based on the num10X_version argument. 

`--gtf:` The gtf to use. Reuqired. By default points to gtf in ref/genes.gtf if it exists. Should match STAR reference.

`--featSTARSolo:` Tells STARSolo how to quantify. By default uses GeneFull (so introns and exons are used), though can use Gene as well (only exons).

## Using with multiplexed samples

This pipeline by default is meant to be run in the setting where there are cells from one individual in a given 10X channel. In many modern experiments, however, many samples from different individuals are mixed together to allow for demultiplexing using genotype, hashing, or similiar approaches. 

In that case one needs to install two addition packages, CellRanger and sinto (https://timoast.github.io/sinto/). One can then run CellRanger on the channel and their tool of choice on the ouptut to get a list of cells per individual. Then all one has to do is run the above pipeline with the --input_bam and --cellFile options. Working on alternative approaches as well that might be cleaner. 

## Older pipeline versions

There are multiple other versions of the pipeline present. The QuantPipeline.no.conda.nf was the one used for the associated publication. 

