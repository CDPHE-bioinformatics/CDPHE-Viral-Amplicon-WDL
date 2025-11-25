version 1.0

import "https://raw.githubusercontent.com/CDPHE-bioinformatics/wdl-shared/b59cb189af2149f00ac0ad04eb3e0813d1cc3971/version_capture_tasks.wdl" as version_capture

task calc_bam_stats_samtools {
    input {
        String sample_name
        File bam
        File bai
    }

    String docker = "staphb/samtools:1.16"

    command <<<
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION
        samtools flagstat ~{bam} > ~{sample_name}_flagstat.txt
        samtools stats ~{bam} > ~{sample_name}_stats.txt
        samtools coverage -m -o ~{sample_name}_coverage_hist.txt ~{bam}
        samtools coverage -o ~{sample_name}_coverage.txt ~{bam}
    >>>

    output {
        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION")
        }

        File flagstat_out = "${sample_name}_flagstat.txt"
        File stats_out = "${sample_name}_stats.txt"
        File covhist_out = "${sample_name}_coverage_hist.txt"
        File cov_out = "${sample_name}_coverage.txt"
    }

    runtime {
        cpu: 2
        memory: "2G"
        disks: "local-disk 1 HDD"
        docker: docker
    }
}

task rename_fasta {
    input {
        String sample_name
        File fasta
    }

    command <<<
        sed 's/>.*/~{sample_name}/' ~{fasta} > ~{sample_name}_consensus_renamed.fa
    >>>

    output {
        File renamed_consensus = "${sample_name}_consensus_renamed.fa"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: "theiagen/utility:1.0"
    }
}

task transfer_outputs {
    meta {
        description: "Transfers files generated in the assembly workflow."
    }

    input {
        String out_dir

        File filtered_reads_1
        File filtered_reads_2
        File seqyclean_summary

        File fastqc_raw1_html
        File fastqc_raw1_zip
        File fastqc_raw2_html
        File fastqc_raw2_zip
        File fastqc_clean1_html      
        File fastqc_clean1_zip        
        File fastqc_clean2_html     
        File fastqc_clean2_zip

        File trimsort_bam
        File trimsort_bamindex

        File variants

        File consensus

        File flagstat_out
        File stats_out
        File covhist_out
        File cov_out

        File renamed_consensus
        File version_capture_file
    }

    String out_dir_path = sub(out_dir, "/$", "") # remove trailing slash

    command <<<
        gsutil -m cp ~{filtered_reads_1} ~{out_dir_path}/seqyclean/
        gsutil -m cp ~{filtered_reads_2} ~{out_dir_path}/seqyclean/
        gsutil -m cp ~{seqyclean_summary} ~{out_dir_path}/seqyclean/
                       
        gsutil -m cp ~{fastqc_raw1_html} ~{out_dir_path}/fastqc/
        gsutil -m cp ~{fastqc_raw1_zip} ~{out_dir_path}/fastqc/
        gsutil -m cp ~{fastqc_raw2_html} ~{out_dir_path}/fastqc/
        gsutil -m cp ~{fastqc_raw2_zip} ~{out_dir_path}/fastqc/
                       
        gsutil -m cp ~{trimsort_bam} ~{out_dir_path}/alignments/
        gsutil -m cp ~{trimsort_bamindex} ~{out_dir_path}/alignments/
                       
        gsutil -m cp ~{variants} ~{out_dir_path}/variants/
                       
        gsutil -m cp ~{consensus} ~{out_dir_path}/assemblies/
                       
        gsutil -m cp ~{flagstat_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{stats_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{covhist_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{cov_out} ~{out_dir_path}/bam_stats/
                       
        gsutil -m cp ~{renamed_consensus} ~{out_dir_path}/assemblies/
        gsutil -m cp ~{version_capture_file} ~{out_dir_path}/summary_results/
                       
        TRANSFER_DATE=$(date)
        echo "$TRANSFER_DATE" | tee TRANSFER_DATE
    >>>

    output {
        String transfer_date = read_string("TRANSFER_DATE")
    }

    runtime {
        cpu: 2
        memory: "2G"
        disks: "local-disk 4 HDD"
        docker: "theiagen/utility:1.0"
    }
}
