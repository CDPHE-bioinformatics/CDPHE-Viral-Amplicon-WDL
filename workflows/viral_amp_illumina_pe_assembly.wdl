version 1.0

import "../tasks/pre_assembly_tasks.wdl"
import "../tasks/assembly_tasks.wdl"
import "../tasks/post_assembly_tasks.wdl"
import "https://raw.githubusercontent.com/CDPHE-bioinformatics/wdl-shared/b59cb189af2149f00ac0ad04eb3e0813d1cc3971/version_capture_tasks.wdl" as version_capture

struct VersionInfo {
    String software
    String docker
    String version
}

workflow viral_amp_illumina_pe_assembly {
    input {
        String project_name
        String sample_name
        String viral_amp_nextclade_organism_id
        File fastq_1
        File fastq_2
        File contam_fasta
        Boolean transfer_results = true
        String? out_dir
        String analysis_date

        File viral_amp_primer_bed
        File viral_amp_ref_fasta
        File? viral_amp_ref_gff

        File calc_percent_coverage_py
    }

    #private declarations
    String version_capture_docker = 'ariannaesmith/cdphe_wdl_version_capture:v0.1.0'
    String ubuntu_docker = "ubuntu:jammy-20240627.1"
    String workflow_name = 'viral_amp_illumina_pe_assembly'
    String workflow_version = 'v1.0.0'
    String workflow_version_und = sub(workflow_version, "\\.", "_")

    String out_dir_path = sub(out_dir, "/$", "") # remove trailing slash
    
    call version_capture.workflow_metadata as w_meta {
        input:
             docker = ubuntu_docker,
             workflow_name = workflow_name
             workflow_version = workflow_version

    }

    ## Pre-assembly tasks    

    call pre_assembly_tasks.filter_reads_seqyclean as filter_reads {
        input:
            contam = contam_fasta,
            sample_name = sample_name,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2
    }

    call pre_assembly_tasks.assess_quality_fastqc as assess_quality {
        input:
           fastq_1 = fastq_1,
           fastq_2 = fastq_2
    }

    call pre_assembly_tasks.align_reads_bwa as align_reads {
        input:
            sample_name = sample_name,
            ref = viral_amp_ref_fasta,
            fastq_1 = filter_reads.cleaned_1,
            fastq_2 = filter_reads.cleaned_2
    }

    ##Assembly tasks

    call assembly_tasks.trim_primers_ivar as trim_primers {
        input:
            sample_name = sample_name,
            primers = viral_amp_primer_bed,
            bam = align_reads.out_bam
    }

    call assembly_tasks.call_variants_ivar as call_variants {
        input:
            sample_name = sample_name,
            ref = viral_amp_ref_fasta,
            gff = viral_amp_ref_gff,
            bam = trim_primers.trimsort_bam
    }

    call assembly_tasks.call_consensus_ivar as call_consensus {
        input:
            sample_name = sample_name,
            ref = viral_amp_ref_fasta,
            bam = trim_primers.trimsort_bam
    }
    
    ##Post assembly tasks

    call post_assembly_tasks.calc_bam_stats_samtools as calc_bam_stats {
        input:
            sample_name = sample_name,
            bam = trim_primers.trimsort_bam,
            bai = trim_primers.trimsort_bamindex
    }

    call post_assembly_tasks.rename_fasta as rename_fasta {
        input:
            sample_name = sample_name,
            fasta = call_consensus.consensus_out
    }

    call post_assembly_tasks.calc_percent_coverage as calc_percent_coverage {
        input:
            sample_name = sample_name,
            fasta = rename_fasta.renamed_consensus,
            reference_file = viral_amp_ref_fasta,
            calc_percent_coverage_py = calc_percent_coverage_py
    }

    call post_assembly_tasks.call_clades_nextclade as call_clades {
        input:
            sample_name = sample_name,
            renamed_consensus = rename_fasta.renamed_consensus,
            organism_id = viral_amp_nextclade_organism_id
    }

    ##VersionInfo structs
    
    VersionInfo seqyclean_version_info = {
        "software": "seqyclean",
        "docker": "staphb/seqyclean:1.10.09",
        "version": "v1.10.09"
    }

    VersionInfo fastqc_version_info = {
        "software": "fastqc",
        "docker": "staphb/fastqc:0.12.1",
        "version": "v0.12.1"
    }

    VersionInfo bwa_version_info = {
        "software": "bwa",
        "docker": "staphb/bwa:0.7.17",
        "version": "v0.7.17"
    }

    VersionInfo samtools_version_info = {
        "software": "samtools",
        "docker": "staphb/samtools:1.19",
        "version": "v1.19"
    }

    VersionInfo ivar_version_info = {
        "software": "ivar",
        "docker": "andersenlabapps/ivar:1.3.1",
        "version": "v1.3.1"
    }

    VersionInfo nextclade_version_info = {
        "software": "nextclade",
        "docker": "ghcr.io/nextstrain/nextclade:3.7.0",
        "version": "v3.7.0"
    }

    VersionInfo version_capture_version_info = {
        "software": "version_capture",
        "docker": version_capture_docker,
        "version": "v0.1.0"
    }

    #call version_capture.workflow_version_capture {
     #   input:
    #}

    Array[VersionInfo] version_array = [
        w_meta.version_info,
        filter_reads.seqyclean_version_info,
        assess_quality.fastqc_version_info,
        align_reads.bwa_version_info,
        align_reads.samtools_version_info,
        call_consensus.ivar_version_info,
        call_consensus.samtools_version_info,
        calc_bam_stats.samtools_version_info,
        call_clades.nextclade_version_info
    ]

    call version_capture.capture_versions as capture_versions {
        input:
            version_array = version_array,
            workflow_name = workflow_name,
            workflow_version = workflow_version,
            project_name = "",
            analysis_date = w_meta.analysis_date,
            docker = version_capture_docker
    }

    call version_capture.task_version_capture as task_version_capture {
        input:
            version_array = version_array,
            workflow_name = "viral_amp_illumina_pe_assembly",
            workflow_version_path = workflow_version_capture.workflow_version_path,
            project_name = project_name,
            analysis_date = workflow_version_capture.analysis_date,
            version_capture_py = version_capture_py
    }

    if (transfer_results) {
      call post_assembly_tasks.transfer_outputs as transfer_outputs {
          input:
              out_dir = "~{out_dir_path}/~{workflow_version_underscores}",
              filtered_reads_1 = filter_reads.cleaned_1,
              filtered_reads_2 = filter_reads.cleaned_2,
              seqyclean_summary = filter_reads.seqyclean_summary,
              fastqc_raw1_html = assess_quality.fastqc1_html,
              fastqc_raw1_zip = assess_quality.fastqc1_zip,
              fastqc_raw2_html = assess_quality.fastqc2_html,
              fastqc_raw2_zip = assess_quality.fastqc2_zip,
              trimsort_bam = trim_primers.trimsort_bam,
              trimsort_bamindex = trim_primers.trimsort_bamindex,
              variants = call_variants.var_out,
              consensus = call_consensus.consensus_out,
              flagstat_out = calc_bam_stats.flagstat_out,
              stats_out = calc_bam_stats.stats_out,
              covhist_out = calc_bam_stats.covhist_out,
              cov_out = calc_bam_stats.cov_out,
              renamed_consensus = rename_fasta.renamed_consensus,
              version_capture_file = [capture_versions.output_file],
      }
    }

    output {
        #String workflow_version = workflow_version_capture.workflow_version
        String workflow_version = w_meta.version_info.version
        String wf_version_und = workflow_version_und

        File filtered_reads_1 = filter_reads.cleaned_1
        File filtered_reads_2 = filter_reads.cleaned_2
        File seqyclean_summary = filter_reads.seqyclean_summary

        File fastqc_raw1_html = assess_quality.fastqc1_html
        File fastqc_raw1_zip = assess_quality.fastqc1_zip
        File fastqc_raw2_html = assess_quality.fastqc2_html
        File fastqc_raw2_zip = assess_quality.fastqc2_zip

        File out_bam = align_reads.out_bam
        File out_bamindex = align_reads.out_bamindex
        String assembler_version = align_reads.assembler_version

        File trim_bam = trim_primers.trim_bam
        File trimsort_bam = trim_primers.trimsort_bam
        File trimsort_bamindex = trim_primers.trimsort_bamindex

        File variants = call_variants.var_out

        File consensus = call_consensus.consensus_out

        File flagstat_out = calc_bam_stats.flagstat_out
        File stats_out = calc_bam_stats.stats_out
        File covhist_out = calc_bam_stats.covhist_out
        File cov_out = calc_bam_stats.cov_out

        File renamed_consensus = rename_fasta.renamed_consensus

        File percent_cvg_csv = calc_percent_coverage.percent_cvg_csv

        File nextclade_csv = call_clades.nextclade_csv
        File nextclade_json = call_clades.nextclade_json

<<<<<<< HEAD
        File version_capture_file = capture_versions.output_file
        String transfer_date_assembly = transfer_outputs.transfer_date
=======
        File version_capture_file = task_version_capture.version_capture_file
        String? transfer_date_assembly = transfer_outputs.transfer_date
>>>>>>> 85ecf0320c3b2844f848fae9e3da70fd4a7aec08
    }
}
