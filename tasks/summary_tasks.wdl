version 1.0

task concatenate_consensus {
    input {
        Array[File] renamed_consensus
    }

    command <<<
        cat ~{sep=" " renamed_consensus} > concatenate_assemblies.fasta
    >>>

    output {
        File cat_fastas = "concatenate_assemblies.fasta"
    }

    runtime {
        cpu: 2
        memory: "6G"
        disks: "local-disk 4 HDD"
        docker: "ubuntu:focal"
    }
}

task summarize_results {
    input {
        String workflow_version
        Array[String] sample_name
        File concat_seq_results_py
        Array[File] cov_out
        String project_name
        String assembler_version
        File workbook_path
    }

    command <<<
        python3 ~{concat_seq_results_py} \
            --workflow_version "~{workflow_version}" \
            --sample_name_array "~{write_lines(sample_name)}" \
            --workbook_path "~{workbook_path}" \
            --cov_out_files "~{write_lines(cov_out)}" \
            --assembler_version "~{assembler_version}" \
            --project_name "~{project_name}" 
    >>>

    output {
        File sequencing_results_csv = "~{project_name}_sequencing_results.csv"
    }

    runtime {
        cpu: 2
        memory: "6G"
        disks: "local-disk 4 HDD"
        docker: "biocontainers/pandas:1.5.1_cv1"
    }
}

task transfer_outputs {
    input {
        String out_dir
        File cat_fastas
        File sequencing_results_csv
    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<
    #remove the mdkir part and add gsutil -m before cp for gcp 
        mkdir -p ~{outdirpath}/multifasta/
        mkdir -p ~{outdirpath}/summary_results/

        cp ~{cat_fastas} ~{outdirpath}/multifasta/
        cp ~{sequencing_results_csv} ~{outdirpath}/summary_results/
    >>>

    runtime {
        cpu: 2
        memory: "6G"
        disks: "local-disk 4 HDD"
        docker: "theiagen/utility:1.0"
    }
}
