version 1.0

import "../tasks/summary_tasks.wdl"
import "../tasks/transfer_task.wdl" as transfer_task
import "https://raw.githubusercontent.com/CDPHE-bioinformatics/wdl-shared/dba3e70cee747617bacbd0312d1de2f6b0731de3/version_capture_tasks.wdl" as version_capture


workflow viral_amp_illumina_pe_summary {
    input {
        Array[String] sample_name
        Array[File?] renamed_consensus
        Array[File?] cov_out # cov as in coverage
        Array[String] out_dir_array
        Array[String] project_name_array
        Array[String?] assembler_version_array
        Array[File] workbook_path_array
        Array[String] workflow_version
        Array[String] workflow_version_und
        Boolean transfer_results = true


        File concat_seq_results_py
    }
    #private declarations
    String version_capture_docker = 'ariannaesmith/cdphe_wdl_version_capture:v1.0.0'
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
            workflow_name = workflow_name,
            workflow_version = workflow_version[0]
    } 

    call summary_tasks.concatenate_consensus as concatenate_consensus {
        input:
            renamed_consensus = select_all(renamed_consensus)
    }

    call summary_tasks.summarize_results as summarize_results {
      input:
        workflow_version = workflow_version[0],
        sample_name = sample_name,
        concat_seq_results_py = concat_seq_results_py,
        cov_out = select_all(cov_out),
        project_name = project_name,
        assembler_version= assembler_version,
        workbook_path = workbook_path
    }

    if (transfer_results) {
      call summary_tasks.transfer_outputs as transfer_outputs {
          input:
              out_dir = "~{out_dir_path}/~{wf_version_und}",
              cat_fastas = concatenate_consensus.cat_fastas,
              sequencing_results_csv = summarize_results.sequencing_results_csv,
      }
    }

    call version_capture.capture_versions as capture_versions {
        input:
            version_array = [w_meta.version_info],
            workflow_name = workflow_name,
            workflow_version = wf_version_und,
            project_name = project_name,
            analysis_date = w_meta.analysis_date,
            docker = version_capture_docker
    }

    output {
        String wf_version_output = w_meta.version_info.version
        File cat_fastas = concatenate_consensus.cat_fastas
        File sequencing_results_csv = summarize_results.sequencing_results_csv
        String? transfer_date = transfer_outputs.transfer_date 
    }
}
