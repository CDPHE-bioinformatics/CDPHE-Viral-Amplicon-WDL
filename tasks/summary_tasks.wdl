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
    meta {
        description: "Transfers files generated in the summary workflow."
    }
    input {
        String out_dir
        File cat_fastas
        File sequencing_results_csv
        File version_capture_file
    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<
        gsutil -m cp ~{cat_fastas} ~{outdirpath}/multifasta/
        gsutil -m cp ~{sequencing_results_csv} ~{outdirpath}/summary_results/
        gsutil -m cp ~{version_capture_file} ~{outdirpath}/versions/

        TRANSFER_DATE=$(date)
        echo "$TRANSFER_DATE" | tee TRANSFER_DATE
    >>>
    
    output {
        String transfer_date = read_string("TRANSFER_DATE")
    }

    runtime {
        cpu: 2
        memory: "6G"
        disks: "local-disk 4 HDD"
        docker: "theiagen/utility:1.0"
    }
}
