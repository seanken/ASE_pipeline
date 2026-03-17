version 1.0

workflow RunMonopogen {
    input {
        File? input_bam
        File? input_bam_index
        Array[File]? input_fastqs
        File? cellranger_reference_tar_gz
        File ref_fasta="gs://ase-methods-dev-wb-sparkly-blueberry-3616/ref/geno_locations.vcf.gz"
        File ref_fasta_index="gs://ase-methods-dev-wb-sparkly-blueberry-3616/ref/geno_locations.vcf.gz.tbi"
        File germline_vcf="gs://ase-methods-dev-wb-sparkly-blueberry-3616/ref/genome.fa"
        File germline_vcf_index="gs://ase-methods-dev-wb-sparkly-blueberry-3616/ref/genome.fa.fai"
        String sample_name="Sample"
    }

    Boolean has_bam_inputs = defined(input_bam) && defined(input_bam_index)
    Boolean has_fastq_inputs = defined(input_fastqs) && defined(cellranger_reference_tar_gz)
    Boolean has_partial_bam_inputs = defined(input_bam) != defined(input_bam_index)
    Boolean has_partial_fastq_inputs = defined(input_fastqs) != defined(cellranger_reference_tar_gz)

    call ValidateInputs {
        input:
            has_bam_inputs = has_bam_inputs,
            has_fastq_inputs = has_fastq_inputs,
            has_partial_bam_inputs = has_partial_bam_inputs,
            has_partial_fastq_inputs = has_partial_fastq_inputs
    }

    if (!defined(input_bam) && !defined(input_bam_index) && defined(input_fastqs) && defined(cellranger_reference_tar_gz)) {
        call RunCellRanger {
            input:
                input_fastqs = select_first([input_fastqs]),
                cellranger_reference_tar_gz = select_first([cellranger_reference_tar_gz]),
                sample_name = sample_name,
                validation_token = ValidateInputs.validation_report
        }
    }

    call MonopogenPreProcess {
        input:
            input_bam = select_first([input_bam, RunCellRanger.output_bam]),
            input_bam_index = select_first([input_bam_index, RunCellRanger.output_bam_index]),
            validation_token = ValidateInputs.validation_report,
    }

    call MonopogenGermline {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            germline_vcf = germline_vcf,
            germline_vcf_index = germline_vcf_index,
            preprocess_output = MonopogenPreProcess.preprocess_output,
            sample_name = sample_name
    }

    output {
        File germline_output = MonopogenGermline.germline_output
        File bam_for_monopogen = select_first([input_bam, RunCellRanger.output_bam])
        File bam_index_for_monopogen = select_first([input_bam_index, RunCellRanger.output_bam_index])
    }
}

task ValidateInputs {
    input {
        Boolean has_bam_inputs
        Boolean has_fastq_inputs
        Boolean has_partial_bam_inputs
        Boolean has_partial_fastq_inputs
    }

    command <<<
        set -euo pipefail

        if [[ "~{has_partial_bam_inputs}" == "true" ]]; then
            echo "Invalid inputs: provide both input_bam and input_bam_index together." >&2
            exit 1
        fi

        if [[ "~{has_partial_fastq_inputs}" == "true" ]]; then
            echo "Invalid inputs: provide both input_fastqs and cellranger_reference_tar_gz together." >&2
            exit 1
        fi

        if [[ "~{has_bam_inputs}" != "true" && "~{has_fastq_inputs}" != "true" ]]; then
            echo "Invalid inputs: provide either (input_bam + input_bam_index) or (input_fastqs + cellranger_reference_tar_gz)." >&2
            exit 1
        fi

        echo "Input validation passed." > validation_report.txt
    >>>

    output {
        File validation_report = "validation_report.txt"
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1 GB"
        disks: "local-disk 10 SSD"
        cpu: 1
    }
}

task RunCellRanger {
    input {
        Array[File] input_fastqs
        File cellranger_reference_tar_gz
        String sample_name
        File validation_token
        Int mem_gb = 80
        Int disk_gb = 200
        Int cpu = 8
    }

    command <<<
        set -euo pipefail

        mkdir -p fastqs
        for fq in ~{sep=' ' input_fastqs}; do
            ln -s "${fq}" "fastqs/$(basename "${fq}")"
        done

        tar -xzf ~{cellranger_reference_tar_gz}
        ref_dir=$(tar -tzf ~{cellranger_reference_tar_gz} | head -n 1 | cut -d'/' -f1)

        cellranger count \
            --id cellranger_out \
            --transcriptome "${ref_dir}" \
            --fastqs fastqs \
            --sample ~{sample_name} \
            --localcores ~{cpu} \
            --localmem ~{mem_gb}
    >>>

    output {
        File output_bam = "cellranger_out/outs/possorted_genome_bam.bam"
        File output_bam_index = "cellranger_out/outs/possorted_genome_bam.bam.bai"
    }

    runtime {
        docker: "10xgenomics/cellranger:latest"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_gb} SSD"
        cpu: cpu
    }
}

task MonopogenPreProcess {
    input {
        File input_bam
        File input_bam_index
        File validation_token
        Int mem_gb = 80
        Int disk_gb = 100
        Int cpu = 1
    }

    command <<<
        path="/opt/Monopogen"
        python3  ${path}/src/Monopogen.py preProcess \
            -b ~{input_bam} \
            -o res_preprocess \
            -a ${path}/apps \
            -t ~{cpu}
        tar -czvf res_preprocess.tar.gz res_preprocess
    >>>

    output {
        File preprocess_output = "res_preprocess.tar.gz"
    }

    runtime {
        docker: "seanken/monopogen:latest"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_gb} SSD"
        cpu: cpu
    }
}

task MonopogenGermline {
    input {
        File ref_fasta
        File ref_fasta_index
        File germline_vcf
        File germline_vcf_index
        File preprocess_output
        String sample_name
        Int mem_gb = 80
        Int disk_gb = 100
        Int cpu = 1
    }

    command <<<
        path="/opt/Monopogen"
        tar -xf ~{preprocess_output}
        python3  ${path}/src/Monopogen.py germline \
            -g ~{ref_fasta} \
            -p ~{germline_vcf} \
            -o out \
            -a ${path}/apps \
            -t ~{cpu}
    >>>

    output {
        File germline_output = "out"
    }

    runtime {
        docker: "seanken/monopogen:latest"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_gb} SSD"
        cpu: cpu
    }
}
