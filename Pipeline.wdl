version 1.0

## Allele-Specific Expression Pipeline
## Converted from Nextflow to WDL
## This was built using Claude Sonnet 4.5 in VSCode with github copilot
## to transform from nextflow to WDL, and then modified by hand by
## Sean Simmomns. 
## Requires: Java, bcftools, tabix, STAR, samtools, bedtools, R, though will change with docker or the like
## To do: 
## -add docker images for each task and add CPU/memory control (build them?)
## -Figure out how to handle downloads from google cloud and add default bucket with jars of interest
## -Figure out tests and test
## -Break code up
## -modify to be better
## -Build ATAC version (?)
##
workflow ASE_Pipeline {
    input {
        # Required inputs
        File input_vcf
        String vcf_col
        Array[File]? input_dirs_R1
        Array[File]? input_dirs_R2
        String? input_dir
        File? input_bam
        File? input_bam_bai ####Add way to assume bai is like bam but with .bai extension
        String? gcs_project
        #String outdir ###Not sure used?

        
        # Optional inputs
        File? cellFile 
        Boolean noQC = false
        
        # STAR parameters
        Int numCells = 3000
        String num10X_version = "v3"
        
        # Resource parameters
        Int numThreads = 8
        String sortRam = "60000000000"
        Int UMILen = 12
        String featSTARSolo = "GeneFull"
        
        # Reference files
        File ref_dir="gs://ase-methods-dev-wb-sparkly-blueberry-3616/refs_moma/STAR_ref" 
        File whitelist="gs://ase-methods-dev-wb-sparkly-blueberry-3616/refs_moma/whitelist_"+num10X_version+"/whitelist.txt"
        File gtf="gs://ase-methods-dev-wb-sparkly-blueberry-3616/refs_moma/genes.gtf"

        # Scripts used by the pipeline
        File AlleleCountJar="gs://ase-methods-dev-wb-sparkly-blueberry-3616/scripts_pipe_moma/AlleleCount.jar"
        File AlleleSNPCountJar="gs://ase-methods-dev-wb-sparkly-blueberry-3616/scripts_pipe_moma/AlleleCountSNP.jar"
        File QCScript="gs://ase-methods-dev-wb-sparkly-blueberry-3616/scripts_pipe_moma/Phased.UMI.QC.R"
        File makeBeds="gs://ase-methods-dev-wb-sparkly-blueberry-3616/scripts_pipe_moma/makeBeds.sh"
        File splitBam="gs://ase-methods-dev-wb-sparkly-blueberry-3616/scripts_pipe_moma/SplitFile.py" 
        
    }
    
    # Validate inputs
    if (!defined(input_dirs_R1) && !defined(input_dirs_R2) && !defined(input_dir) && !defined(input_bam)) {
        call Error as Error_missing_input { input: msg = "Need to pass one of input_dirs_R1/input_dirs_R2, input_dir, or input_bam" }
    }
    
    if ((defined(input_dirs_R1) || defined(input_dirs_R2) || defined(input_dir)) && defined(input_bam)) {
        call Error as Error_both_inputs { input: msg = "Can only pass one of input FASTQ options (input_dirs_R1/input_dirs_R2 or input_dir) or input_bam" }
    }

    if (defined(input_dirs_R1) != defined(input_dirs_R2)) {
        call Error as Error_unpaired_fastq_inputs { input: msg = "Need both input_dirs_R1 and input_dirs_R2" }
    }
    
    if (defined(cellFile) && !defined(input_bam)) {
        call Error as Error_cellfile_without_bam { input: msg = "If cell file need bam file" }
    }
    
    # Prepare VCF
    call PrepVCF {
        input:
            vcf = input_vcf,
            vcf_col = vcf_col
    }
    
    # Get cells if specified #### Is select first needed? Can we remove due to if statement?
    if (defined(cellFile) && defined(input_bam) && defined(input_bam_bai) ) {
        #bai=select_first([input_bam_bai, (input_bam + ".bai")])
        call GetCells {
            input:
                bam = select_first([input_bam]),
                bai = select_first([input_bam_bai]),
                cells = select_first([cellFile]),
                splitBam = splitBam,
                vcf_col = vcf_col
        }
    }
    
    
    
    # Process input BAM to FASTQ if needed
    if (defined(input_bam)) {
        File bam_to_use = select_first([GetCells.output_bam, input_bam])
        
        call ProcessInputBam {
            input:
                bam = bam_to_use
        }
    }

    if (defined(input_dir) && !defined(input_dirs_R1) && !defined(input_dirs_R2)) {
        call FindFastqsByPrefix {
            input:
                input_dir = select_first([input_dir]),
                gcs_project = select_first([gcs_project])
        }
    }

    Array[File] read1_fastqs = select_first([input_dirs_R1, FindFastqsByPrefix.read1_fastqs, ProcessInputBam.read1_fastqs, []])
    Array[File] read2_fastqs = select_first([input_dirs_R2, FindFastqsByPrefix.read2_fastqs, ProcessInputBam.read2_fastqs, []])
    

    # Run STAR Solo
    call RunSTARSolo {
        input:
            read1_fastqs = read1_fastqs,
            read2_fastqs = read2_fastqs,
            vcf = PrepVCF.new_vcf,
            numCells = numCells,
            ref_dir = ref_dir,
            whitelist = whitelist,
            numThreads = numThreads,
            UMILen = UMILen,
            feat = featSTARSolo,
            sortRam = sortRam
    }
    
    # Count alleles at gene level
    call CountAlleles {
        input:
            star_output = RunSTARSolo.output_dir,
            jar = AlleleCountJar
    }
    
    # SNP level processing
    call SNPLevel {
        input:
            gtf = gtf,
            vcf = PrepVCF.new_vcf,
            makeBeds = makeBeds
    }
    
    # Count alleles at SNP level
    call CountAllelesSNP {
        input:
            star_output = RunSTARSolo.output_dir,
            jar = AlleleSNPCountJar
    }
    
    # Isoform analysis if requested
    ##HAS BEEN REMOVED, MIGHT ADD IN LATER!

    
    # QC if requested
    if (!noQC) {
        call PhasedUMI_QC {
            input:
                script = QCScript,
                counts = CountAlleles.counts,
                star_output = RunSTARSolo.output_dir,
                feat = featSTARSolo,
                counts_snp = CountAllelesSNP.counts
        }
    }
    
    ##Add STARSolo output dir to outputs?
    output {
        File gene_counts = CountAlleles.counts
        File snp_counts = CountAllelesSNP.counts
        File? qc_hist = PhasedUMI_QC.hist
        File? qc_basic = PhasedUMI_QC.basic_qc
        File? qc_umi = PhasedUMI_QC.umi_counts
    }
}

task Error {
    input {
        String msg
    }
    command {
        echo "ERROR: ${msg}" >&2
        exit 1
    }
    runtime {
        docker: "ubuntu:20.04"
    }
}

##
##Finds and localizes FASTQ files in GCS from a shared prefix, then separates into R1/R2 arrays
##
task FindFastqsByPrefix {
    input {
        String input_dir
        String gcs_project
    }

    command <<<
        set -e
        mkdir -p localized_fastqs/R1 localized_fastqs/R2

        gsutil -u "~{gcs_project}" ls "~{input_dir}*R1*fastq.gz" > r1_uris.txt
        gsutil -u "~{gcs_project}" ls "~{input_dir}*R2*fastq.gz" > r2_uris.txt

        if [ ! -s r1_uris.txt ]; then
            echo "No R1 FASTQ files found for prefix: ~{input_dir}" >&2
            exit 1
        fi

        if [ ! -s r2_uris.txt ]; then
            echo "No R2 FASTQ files found for prefix: ~{input_dir}" >&2
            exit 1
        fi

        r1_count=$(wc -l < r1_uris.txt)
        r2_count=$(wc -l < r2_uris.txt)
        if [ "$r1_count" -ne "$r2_count" ]; then
            echo "Mismatched FASTQ counts for prefix: ~{input_dir} (R1=$r1_count, R2=$r2_count)" >&2
            exit 1
        fi

        gsutil -m cp $(cat r1_uris.txt) localized_fastqs/R1/
        gsutil -m cp $(cat r2_uris.txt) localized_fastqs/R2/
    >>>

    output {
        Array[File] read1_fastqs = glob("localized_fastqs/R1/*")
        Array[File] read2_fastqs = glob("localized_fastqs/R2/*")
    }

    runtime {
        memory: "16 GB"
        cpu: 1
        docker: "google/cloud-sdk:latest"
    }
}


##
##Extracts reads from input BAM that correspond to cells of interest
##
task GetCells {
    input {
        File bam
        File bai
        File cells
        File splitBam
        String vcf_col
    }
    
    command <<<
        grep ~{vcf_col}$ ~{cells} | awk '{print $1}' | sed 's/$/\toutput/g' > cells.sinto.txt
        python ~{splitBam} ~{bam} cells.sinto.txt
    >>>
    
    output {
        File output_bam = "output.bam"
    }
    
    runtime {
        cpu: 1
        docker: "intelliseqngs/pysam"
        memory: "16 GB"
    }
}

##
##Transforms an input BAM into FASTQ files using cellranger bamtofastq
##
task ProcessInputBam {
    input {
        File bam
    }
    
    command <<<
        cellranger bamtofastq ~{bam} FA_DS
    >>>
    
    output {
        Array[File] read1_fastqs = glob("FA_DS/*/*_R1*fastq.gz")
        Array[File] read2_fastqs = glob("FA_DS/*/*_R2*fastq.gz")
    }
    
    runtime {
        memory: "32 GB"
        docker: "litd/docker-cellranger:v8.0.1"
        cpu: 4
    }
}

##
##Prepares VCF by subsetting to only heterozygous variants for the sample of interest
##
task PrepVCF {
    input {
        File vcf
        String vcf_col
    }
    ##Modify to use bcftools to filter genotypes
    command <<<
        tabix -p vcf ~{vcf}
        bcftools view -H -O v -s ~{vcf_col} -g het ~{vcf} > new.vcf
    >>>
    
    output {
        File new_vcf = "new.vcf"
    }
    
    runtime {
        memory: "30 GB"
        cpu: 1
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    }
}




##
##Processes input data with STARSolo to get quantificaition, and add information about allelic origin of reads
##to output BAM, uses WASP for avoiding reference bias
##
task RunSTARSolo {
    input {
        Array[File] read1_fastqs
        Array[File] read2_fastqs
        File vcf
        Int numCells
        File ref_dir
        File whitelist
        Int numThreads
        Int UMILen
        String feat
        String sortRam
    }
    
    command <<<
        mkdir output
        STAR --genomeDir ~{ref_dir} \
            --readFilesIn ~{sep=',' read1_fastqs} ~{sep=',' read2_fastqs} \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist ~{whitelist} \
            --soloUMIlen ~{UMILen} \
            --soloUMIfiltering MultiGeneUMI_CR \
            --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM vG vA vG \
            --outSAMtype BAM SortedByCoordinate \
            --soloCellFilter EmptyDrops_CR ~{numCells} 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \
            --runThreadN ~{numThreads} \
            --outFileNamePrefix output/results \
            --readFilesCommand zcat \
            --varVCFfile ~{vcf} \
            --waspOutputMode SAMtag \
            --soloFeatures ~{feat} \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize 1500000000 1500000000 \
            --limitBAMsortRAM ~{sortRam} \
            --outFilterScoreMin 30 \
            --soloUMIdedup 1MM_CR \
            --clipAdapterType CellRanger4
    >>>
    
    output {
        File output_dir = "output"
    }
    
    runtime {
        cpu: numThreads
        memory: "90 GB"
        docker: "cumulusprod/starsolo:2.7.10b"
    }
}

##
##Quantify ASE at Gene level, using reads in gene and using phasing in VCF
##
task CountAlleles {
    input {
        File star_output
        File jar
    }
    
    command <<<
        java -jar ~{jar} ~{star_output}/resultsAligned.sortedByCoord.out.bam counts.txt
    >>>
    
    output {
        File counts = "counts.txt"
    }
    
    runtime {
        memory: "40 GB"
        cpu: 1
        docker: "openjdk:11"
    }
}

##
##Quantify ASE at SNP level, just using reads overlapping that SNP
##
task CountAllelesSNP {
    input {
        File star_output
        File jar
    }
    
    command <<<
        java -jar ~{jar} ~{star_output}/resultsAligned.sortedByCoord.out.bam counts.txt
    >>>
    
    output {
        File counts = "counts.txt"
    }
    
    runtime {
        memory: "40 GB"
        cpu: 1
        docker: "openjdk:11"
    }
}

##
##Create bed files for SNPs with info about phasing
##
task SNPLevel {
    input {
        File gtf
        File vcf
        File makeBeds
    }
    
    command <<<
        source ~{makeBeds} ~{gtf} ~{vcf}
        bedtools intersect -a snps.bed -b gene.bed -wa -wb | \
            awk '{print $4"\t"$9"\t"$5}' | sort | uniq > comb.bed
    >>>
    
    output {
        File comb_bed = "comb.bed"
        File snps_bed = "snps.bed"
    }
    
    runtime {
        memory: "20 GB"
        cpu: 1
        docker: "biocontainers/bedtools:v2.29.2-1-deb_cv1"
    }
}

##
##Gets QC realted to the UMI phasing, returns as tables and plots
##
task PhasedUMI_QC {
    input {
        File script
        File counts
        File star_output
        String feat
        File counts_snp
    }
    
    command <<<
        Rscript ~{script} ~{counts} ~{star_output} ~{feat} ~{counts_snp}
    >>>
    
    output {
        File hist = "hist.ratio.pdf"
        File basic_qc = "Basic.QC.txt"
        File umi_counts = "UMI.counts.by.gene.txt"
    }
    
    runtime {
        memory: "30 GB"
        cpu: 1
        docker: "rocker/tidyverse:4"
    }
}

