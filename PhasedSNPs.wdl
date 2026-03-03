version 1.0

workflow RunMonopogen {
    input {
        File input_bam
        File input_bam_index
        File ref_fasta
        File ref_fasta_index
        File germline_vcf
        File germline_vcf_index
        String sample_name
    }

    call MonopogenPreProcess {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
    }

    call MonopogenGermline {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            germline_vcf = germline_vcf,
            germline_vcf_index = germline_vcf_index,
            preprocess_output = MonopogenPreProcess.preprocess_output,
    }

    output {
        File germline_output = MonopogenGermline.germline_output
    }
}

task MonopogenPreProcess {
    input {
        File input_bam
        File input_bam_index
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
    >>>

    output {
        File preprocess_output = "res_preprocess"
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
        Int mem_gb = 80
        Int disk_gb = 100
        Int cpu = 1
    }

    command <<<
        path="/opt/Monopogen"
        python3  ${path}/src/Monopogen.py germline \
            -g ~{ref_fasta} \
            -p ~{germline_vcf} \
            -o ~{preprocess_output} \
            -a ${path}/apps \
            -t ~{cpu}
    >>>

    output {
        File germline_output = "~{sample_name}_germline"
    }

    runtime {
        docker: "seanken/monopogen:latest"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_gb} SSD"
        cpu: cpu
    }
}
