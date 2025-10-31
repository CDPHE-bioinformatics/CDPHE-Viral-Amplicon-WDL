version 1.0

import "https://raw.githubusercontent.com/CDPHE-bioinformatics/wdl-shared/dba3e70cee747617bacbd0312d1de2f6b0731de3/version_capture_tasks.wdl" as version_capture

import "../tasks/transfer_task.wdl" as transfer_task

workflow viral_amp_wwt_variant_calling {

    input {

        Array[File] trimsort_bam
        Array[String] sample_name
        Array[String] out_dir_array
        Boolean overwrite = true
        Boolean transfer_results = true
        Array[String] project_name_array 
        Array[String] freyja_pathogen
        Array[String] workflow_version
        Array[String] workflow_version_und

        # reference files/workspace data
        File reference_genome
        File? reference_gff

    }
    #private declarations
    String version_capture_docker = 'ariannaesmith/cdphe_wdl_version_capture:v1.0.0'
    String workflow_name = 'viral_amp_wwt_variant_calling'
    String wf_version = select_first(workflow_version)
    String wf_version_und = select_first(workflow_version_und)
    String out_dir_path = sub(out_dir_array[0], "/$", "") # remove trailing slash

    # secret variables
    String project_name = project_name_array[0]
    #String out_dir = select_first([out_dir_array])[0]
    String pathogen = select_first(freyja_pathogen)

    scatter (id_bam in zip(sample_name, trimsort_bam)) {

        call variant_calling {
            input:
                bam = id_bam.right,
                ref = reference_genome,
                ref_gff = reference_gff,
                sample_name = id_bam.left
        }

        call freyja_demix {
            input:
                variants = variant_calling.variants,
                depth = variant_calling.depth,
                sample_name = id_bam.left,
                freyja_pathogen = pathogen
        }
        
        call mutations_tsv {
            input:
                variants = variant_calling.variants,
                sample_name = id_bam.left,
                project_name = project_name
        }
    }

    call freyja_aggregate {
        input:
            demix = select_all(freyja_demix.demix)
    }

    call combine_mutations_tsv {
        input:
            mutations = mutations_tsv.mutations
    }

    call version_capture.workflow_metadata as w_meta {
        input:
            docker = version_capture_docker,
            workflow_name = workflow_name,
            workflow_version = workflow_version[0]

    }

    call version_capture.capture_versions as version_cap {
        input:
            version_array = [],
            workflow_name = workflow_name,
            workflow_version = wf_version_und,
            project_name = project_name,
            analysis_date = w_meta.analysis_date,
            docker = version_capture_docker
    }
    
    SubdirsToFiles subdirs_to_files = object { subdirs_to_files: [
        ("viral_amp_wwt_variant_calling/freyja",
            flatten([
                variant_calling.variants,
                variant_calling.depth,
                select_all(freyja_demix.demix)
            ])
        ),
        ("viral_amp_wwt_variant_calling", [
            combine_mutations_tsv.combined_mutations_tsv,
            freyja_aggregate.demix_aggregated
        ])
    ]}

    if (transfer_results) {
        call transfer_task.transfer as transfer_set_results {
            input:
                out_dir = "~{out_dir_path}/~{wf_version_und}",
                overwrite = overwrite,
                cpu = 8,
                subdirs_to_files = subdirs_to_files
        }
    }

    output {
        Array[File] variants = variant_calling.variants
        Array[File] depth = variant_calling.depth
        Array[File] demix = select_all(freyja_demix.demix)
        File demix_aggregated = freyja_aggregate.demix_aggregated
        File combined_mutations_tsv = combine_mutations_tsv.combined_mutations_tsv
        String? transfer_date_viral_amp_wwt_variant_calling = transfer_set_results.transfer_date
    }
}

task variant_calling {
    input {
        File bam
        File ref
        File? ref_gff
        String sample_name
    }

    command <<<

    # grab ivar and samtools versions
    ivar version | awk '/version/ {print $3}' | tee VERSION_ivar
    samtools --version | awk '/samtools/ {print $2}' | tee VERSION_samtools

    samtools mpileup -A -aa -d 600000 -B -Q 20 -q 0 -f ~{ref} ~{bam} | tee >(cut -f1-4 > ~{sample_name}_depth.tsv) | \
    ivar variants -p ~{sample_name}_variants.tsv -q 20 -t 0.0 -r ~{ref} -g ~{ref_gff}
    
    >>>

    output {
        File variants = "~{sample_name}_variants.tsv"
        File depth = "~{sample_name}_depth.tsv"
        String samtools_version_andersenlabapps = read_string("VERSION_samtools")
        String ivar_version = read_string("VERSION_ivar")
    }

     runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    } 
}

task freyja_demix {
    input {
        String sample_name
        File variants
        File depth
        String freyja_pathogen
    }

    command <<<

        freyja --version | awk '{print $NF}' | tee VERSION
        # $NF refers to the last field split by white spaces


        freyja demix --eps 0.01 --covcut 10 --pathogen ~{freyja_pathogen} --depthcutoff 10 ~{variants} ~{depth} --output ~{sample_name}_demixed.tsv
    >>>

    output {
        File? demix = "${sample_name}_demixed.tsv"
        String freyja_version = read_string("VERSION")
    }

    runtime {
        docker: "staphb/freyja:2.0.1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 200 SSD"
        continueOnReturnCode: [0, 1]
    }
}

task mutations_tsv {
    input {
        String sample_name
        String project_name
        File variants
    }

    command <<<
        #add columns with sample_name and project_name
        paste ~{variants} <(yes ~{sample_name} | head -n $(cat ~{variants} | wc -l)) <(yes ~{project_name} | head -n $(cat ~{variants} | wc -l)) > ~{sample_name}_mutations.tsv
        sed -i -e '1s/REGION/ref_genome/' -e '1s/POS/position/' -e '1s/REF/ref_nucl/' -e '1s/ALT/alt_nucl/' -e '1s/REF_DP/ref_depth/' -e '1s/REF_QUAL/ref_qual/' -e '1s/REF_CODON/ref_codon/' -e '1s/REF_AA/ref_aa/' -e '1s/ALT_DP/alt_depth/' -e '1s/ALT_QUAL/alt_qual/' -e '1s/ALT_CODON/alt_codon/' -e '1s/ALT_AA/alt_aa/' -e '1s/ALT_FREQ/alt_freq/' -e '1s/TOTAL_DP/total_depth/' -e '1s/PVAL/pval/' -e '1s/PASS/pass/' -e '1s/GFF_FEATURE/gff_feature/' -e '1s/~{sample_name}/sample_name/' -e '1s/~{project_name}/project_name/' ~{sample_name}_mutations.tsv
    >>>

    output {
        File mutations = "${sample_name}_mutations.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 500 HDD"
    }
}

task freyja_aggregate {
    input {
        Array[File] demix
    }

    command <<<
        mkdir demix_outputs
        mv ~{sep=' ' demix} demix_outputs/
        freyja aggregate demix_outputs/ --output demix_aggregated.tsv
    >>>

    output {
        File demix_aggregated = "demix_aggregated.tsv"
    }

    runtime {
        docker: "staphb/freyja:2.0.1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 200 SSD"
    }
}

task combine_mutations_tsv {
    input {
        Array[File] mutations
    }

    command <<<
        # combine the coutns and frequency files for all samples into one
        awk 'FNR==1 && NR!=1{next;}{print}' ~{sep=' ' mutations} >> combined_mutations.tsv
    >>>

    output {
        File combined_mutations_tsv = "combined_mutations.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 500 HDD"
    }
}
