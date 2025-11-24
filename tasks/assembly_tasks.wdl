version 1.0

import "https://raw.githubusercontent.com/CDPHE-bioinformatics/wdl-shared/dba3e70cee747617bacbd0312d1de2f6b0731de3/version_capture_tasks.wdl" as version_capture

task align_reads_bwa {
    input {
        File fastq_1
        File fastq_2
        File ref
        String sample_name
    }

    String docker = "quay.io/broadinstitute/viral-core:2.2.3"

    command <<<
        bwa 2>&1 | awk '/Version/{print $2}' | tee VERSION_BWA
        samtools --version | awk '/samtools / {print $2}' | tee VERSION_SAMTOOLS
        bwa index -p reference.fasta -a is ~{ref}
        bwa mem -t 2 reference.fasta ~{fastq_1} ~{fastq_2} | \
        samtools sort | \
        samtools view -u -h -F 4 -o ./~{sample_name}_aln.sorted.bam
        samtools index ./~{sample_name}_aln.sorted.bam
    >>>

    output {
        VersionInfo bwa_version_info = object {
            software: "bwa",
            docker: docker,
            version: read_string("VERSION_BWA")
        }

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION_SAMTOOLS")
        }

        File out_bam = "${sample_name}_aln.sorted.bam"
        File out_bamindex = "${sample_name}_aln.sorted.bam.bai"
        String assembler_version = read_string("VERSION_BWA")
    }

    runtime {
        cpu: 2
        memory: "2G"
        disks: "local-disk 2 HDD"
        docker: docker
    }
}

task trim_primers_ivar {
    input {
        File primers
        File bam
        String sample_name
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    command {
        ivar version | awk '/version/ {print $3}' | tee VERSION_IVAR
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION_SAMTOOLS
        ivar trim -e -i ${bam} -b ${primers} -p ${sample_name}_trim.bam
        samtools sort ${sample_name}_trim.bam -o ${sample_name}_trim.sort.bam
        samtools index ${sample_name}_trim.sort.bam
    }

    output {
        VersionInfo ivar_version_info = object {
            software: "ivar",
            docker: docker,
            version: read_string("VERSION_IVAR")
        }

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION_SAMTOOLS")
        } 

        File trim_bam = "${sample_name}_trim.bam"
        File trimsort_bam = "${sample_name}_trim.sort.bam"
        File trimsort_bamindex = "${sample_name}_trim.sort.bam.bai"
    }

    runtime {
        cpu: 6
        memory: "16G"
        disks: "local-disk 1 HDD"
        maxRetries: 2
        bootDiskSizeGb: 10
        docker: docker 
    }
}

task call_variants_ivar {
    input {
        String sample_name
        File ref
        File gff
        File bam
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    command <<<
        ivar version | awk '/version/ {print $3}' | tee VERSION_IVAR
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION_SAMTOOLS
        
        samtools faidx ~{ref}
        samtools mpileup -A -aa -d 600000 -B -Q 30 -q 30 -f ~{ref} ~{bam} | \
        ivar variants -p ~{sample_name}_variants -q 30 -t 0.6 -m 10 -r ~{ref} ~{if defined(gff) then "-g " + gff else ""}
    >>>

    output {
        VersionInfo ivar_version_info = object {
            software: "ivar",
            docker: docker,
            version: read_string("VERSION_IVAR")
        }

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION_SAMTOOLS")
        }

        File var_out = "${sample_name}_variants.tsv"
    }

    runtime {
        cpu: 6
        memory: "16G"
        disks: "local-disk 1 HDD"
        maxRetries: 2
        bootDiskSizeGb: 10
        docker: docker
    }
}

task call_consensus_ivar {
    input {
        String sample_name
        File ref
        File bam
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    command <<<
        ivar version | awk '/version/ {print $3}' | tee VERSION_IVAR
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION_SAMTOOLS
        samtools faidx ~{ref}
        samtools mpileup -A -aa -d 600000 -B -Q 30 -q 30 -f ~{ref} ~{bam} | \
        ivar consensus -p ~{sample_name}_consensus -q 30 -t 0.6 -m 10
    >>>

    output {
        VersionInfo ivar_version_info = object {
            software: "ivar",
            docker: docker,
            version: read_string("VERSION_IVAR")
        }

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION_SAMTOOLS")
        }

        File consensus_out = "${sample_name}_consensus.fa"
    }

    runtime {
        cpu: 6
        memory: "16G"
        disks: "local-disk 1 HDD"
        maxRetries: 2
        bootDiskSizeGb: 10
        docker: docker
    }
}