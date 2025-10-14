version 1.0

import "../tasks/summary_tasks.wdl"
import "../tasks/version_capture_tasks.wdl"

workflow viral_amp_illumina_pe_summary {
    input {
        Array[String] sample_name
        Array[File?] renamed_consensus
        Array[File?] cov_out # cov as in coverage
        Array[File?] percent_cvg_csv
        Array[File?] nextclade_csv
        Array[String]? out_dir_array
        Boolean transfer_results = true
        Array[String] project_name_array
        Array[String?] assembler_version_array
        Array[File] workbook_path_array
        Array[String] workflow_version
        Array[String] workflow_version_und

        File concat_seq_results_py
    }
    #private declarations
    String version_capture_docker = 'ariannaesmith/cdphe_wdl_version_capture:v0.1.0'
    String workflow_name = 'viral_amp_illumina_pe_summary'
    

    String project_name = project_name_array[0]
    File workbook_path = workbook_path_array[0]
    String assembler_version = select_all(assembler_version_array)[0]
    String out_dir_path = sub(out_dir_array[0], "/$", "") # remove trailing slash
    
    String wf_version = select_first(workflow_version)
    String wf_version_und = select_first(workflow_version_und)

     call version_capture.workflow_metadata as w_meta {
        input:
             docker = version_capture_docker,
             workflow_name = workflow_name
             workflow_version = workflow_version

    } 

    call version_capture_tasks.workflow_version_capture {
        input:
    }

    call summary_tasks.concatenate_consensus as concatenate_consensus {
        input:
            renamed_consensus = select_all(renamed_consensus)
    }

    call summary_tasks.summarize_results as summarize_results {
      input:
        workflow_version = workflow_version_capture.workflow_version,
        sample_name = sample_name,
        concat_seq_results_py = concat_seq_results_py,
        nextclade_csv = select_all(nextclade_csv),
        cov_out = select_all(cov_out),
        percent_cvg_csv = select_all(percent_cvg_csv),
        project_name = project_name,
        assembler_version= assembler_version,
        workbook_path = workbook_path
    }

    if (transfer_results) {
        call summary_tasks.transfer_outputs as transfer_outputs {
            input:
                out_dir = "~{out_dir_path}/summary_results/assembly/~{version_capture.workflow_version_path}",
                cat_fastas = concatenate_consensus.cat_fastas,
                sequencing_results_csv = summarize_results.sequencing_results_csv
        }
    }

    call version_capture.capture_versions as version_cap {
        input:
            version_array = version_array,
            workflow_name = workflow_name,
            workflow_version = workflow_version_und,
            project_name = project_name,
            analysis_date = w_meta.analysis_date,
            docker = version_capture_docker
    }

    output {
        String workflow_version = workflow_version_capture.workflow_version

        File cat_fastas = concatenate_consensus.cat_fastas

        File sequencing_results_csv = summarize_results.sequencing_results_csv
    }
}
