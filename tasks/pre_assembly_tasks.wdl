version 1.0

import "https://raw.githubusercontent.com/CDPHE-bioinformatics/wdl-shared/b59cb189af2149f00ac0ad04eb3e0813d1cc3971/version_capture_tasks.wdl" as vc

task filter_reads_seqyclean {
    input {
        File contam
        String sample_name
        File fastq_1
        File fastq_2
    }

    String docker = "staphb/seqyclean:1.10.09"

    command <<<
        seqyclean -h | awk '/Version/ {print $2}' | tee VERSION
        seqyclean -minlen 25 -qual 30 30 -gz -1 ~{fastq_1} -2 ~{fastq_2} -c ~{contam} -o ~{sample_name}_clean
    >>>

    output {
        VersionInfo seqyclean_version_info = object {
            software: "seqyclean",
            docker: docker,
            version: read_string("VERSION")
        }

        File cleaned_1 = "${sample_name}_clean_PE1.fastq.gz"
        File cleaned_2 = "${sample_name}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: docker
    }
}

task assess_quality_fastqc {
    input {
        File fastq_1
        File fastq_2
    }

    String docker = "staphb/fastqc:0.11.9"
    String fastq1_name = basename(basename(basename(fastq_1, ".gz"), ".fastq"), ".fq")
    String fastq2_name = basename(basename(basename(fastq_2, ".gz"), ".fastq"), ".fq")

    command <<<
        fastqc --version | awk '/FastQC/ {print $2}' | tee VERSION
        fastqc --outdir "$PWD" ~{fastq_1} ~{fastq_2}
    >>>

    output {
        VersionInfo fastqc_version_info = object {
            software: "fastqc",
            version: read_string("VERSION"),
            docker: docker
        }

        File fastqc1_html = "${fastq1_name}_fastqc.html"
        File fastqc1_zip = "${fastq1_name}_fastqc.zip"
        File fastqc2_html = "${fastq2_name}_fastqc.html"
        File fastqc2_zip = "${fastq2_name}_fastqc.zip"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: docker
    }
}

