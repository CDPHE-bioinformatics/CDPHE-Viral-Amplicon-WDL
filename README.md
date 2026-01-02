# CDPHE_Viral_Amp_WDL Workflows

## Disclaimer
**Next generation sequencing and bioinformatic and genomic analysis at the Colorado Department of Public Health and Environment (CDPHE) is not CLIA validated at this time. These workflows and their outputs are not to be used for diagnostic purposes and should only be used for public health action and surveillance purposes. CDPHE is not responsible for the incorrect or inappropriate use of these workflows or their results.**

<br/>

## Overview

The following documentation describes the Colorado Department of Public Health and Environment's workflows for the assembly and analysis of viral amplicon whole genome and targeted sequencing data on GCP's Terra.bio platform. Workflows are written in WDL and can be imported into a Terra.bio workspace through dockstore (see Setup section below: https://dockstore.org/).

Our viral reference-based assembly workflows are highly adaptable and facilitate the assembly and analysis of tiled amplicon based sequencing data of viral samples. The workflows can accommodate various amplicon primer schemes sequenced using Illumina paired-end sequencing. These workflows can be applied to both (1) whole-genome tiled amplicon approaches and (2) targeted sequencing approaches using known primer locations to amplify specific regions of the viral genome.

<br/>

### Available Workflows

<br/>

|Workflow Name | Description |
|--------------|-------------|
| ``viral_amp_illumina_pe_assembly`` | Performs reference-based assembly of viral genomes or targeted regions from Illumina paired-end amplicon data. |
| ``viral_amp_illumina_pe_summary`` | Generates summary statistics and quality metrics from assembled viral genomes or targeted regions. |
| ``viral_amp_wwt_variant_calling`` | Uses Freyja to estimate relative lineage abundances and variant composition from wastewater samples. |

### Process

The first, ``viral_amp_illumina_pe_assembly``, takes raw paired-end Illumina reads and performs quality control, contamination filtering, primer trimming, reference-guided assembly, variant calling, and consensus genome generation. Intermediate files and consensus sequences can be transferred to a designated Google Cloud bucket(GCP) for storage and downstream access (optional).

Next, the ``viral_amp_illumina_pe_summary`` workflow aggregates the outputs from multiple samples, concatenates consensus sequences, and generates a comprehensive sequencing results report (including coverage metrics), while also organizing the outputs into a versioned results directory.

Finally, if wastewater samples are being analyzed, the ``viral_amp_wwt_variant_calling workflow`` applies Freyja to estimate relative lineage abundances in wastewater samples, accounting for the mixed nature of viral populations present. 


## Setup

<br/>

### Import Workflow from Dockstore
To use the workflow on the Terra.bio platform, first you will need to import the workflow from Dockstore. All workflows can be found under our dockstore organization called CDPHE-bioinformatics.
1. Go to dockstore (https://dockstore.org/).
2. Along the top search bar click on Organizations and search for "CDPHE".
3. Select the workflow.
4. On the right hand side of the workflow description, select "Launch with Terra".
5. Select the Destination workspace and select "Import". 
6. The workflow will now be displayed as a card under your workflows tab in your Terra workspace. 

<br/>

### Workspace Data
Prior to running any of the workflows, you must set up the Terra workspace data with the correct reference files and custom python scripts. Python scripts can be found in the ``scripts`` directory. Workspace variables are named using the following format ``{organism}_{description}_{file_type}``, except for the primer bed files which are named as ``{description}_{file_type}``. Reference files and python scripts should be saved in a GCP bucket. The GCP bucket path to the file will serve as the "value" when adding data to the terra workspace data table. Alternatively, reference files and scripts can be uploaded and saved as workspace files on Terra.bio, and the backend Terra.bio GCP bucket address can be used.

To add data to the terra workspace data:
1. Navigate to the Data tab in your Terra workspace.
2. In the left hand list of data tables, under "Other Data" select "Workspace Data".
3. Click on the "+" button in the lower right hand corner of the workspace data table. 
4. Fill in the "Key" column with the workspace variable name, the "Value" column with GCP bucket path to the file and the "Description" column with a brief description if desired. 
5. Once complete hit the check mark to the right.

Below is a data table detailing the workspace data you will need to set up in order to run the Viral workflows. 

| workspace variable name | workflow | file name | description |
|-------------------------|------------|----------------|-----------------|
| ``adapters_and_contaminants`` | ``viral_amp_illumina_pe_assembly`` | Adapters_plus_PhiX_174.fasta | Adapters and PhiX contaminant sequences removed during FASTQ cleaning and filtering using SeqyClean. Thanks to Erin Young at Utah Public Health Laboratory for providing this file! |
| ``<your_primers_bed>`` | ``viral_amp_illumina_pe_assembly`` | <your_primers>.bed | Primer BED file for the amplicon set used during genome or targeted assembly. |
| ``<your_ref>_fasta`` | ``viral_amp_illumina_pe_assembly`` | <your_ref>.fasta | Reference genome FASTA that correpsonds to your primer bed file to use for read mapping. |
| ``<your_ref>_gff`` | ``viral_amp_illumina_pe_assembly`` | <your_ref>.gff | Genome annotation file (GFF) corresponding to the reference fasta used for viral amplicon assembly. If you are running the wastewater workflow. This reference fasta need to match the reference used by Freyja. |
| ``concat_seq_results_py`` | ``viral_amp_illumina_pe_assembly`` | concat_seq_results.py | Python script used to concatenate sequence metrics and results from the Viral Amplicon assembly workflow. |


## viral_amp_illumina_pe_assembly

File: viral_amp_illumina_pe_assembly.wdl

This workflow was developed for the assembly of Illumina paired-end read data. The workflow accepts "sample" as the root entity type. The workflow will:

1. Human Read Scrubbing: Perform human read scrubbing on the fastq files using [hostile](https://github.com/bede/hostile) if `scrub_reads` is set to true.
2. Quality Filtering & Trimming: Use Seqyclean to quality filter and trim raw fastq files (or human read scrubbed fastq files, if `scrub_reads` was set to true).
    - Seqyclean parameters include a minimum read length set to 25 bp and quality trimming set to a minimum Phred quality score of 30.
3. Quality Assessment: Run FastQC on both the raw and cleaned reads
4. Read Alignment: Align reads to the reference genome using bwa and then sort the bam by coordinates using Samtools
5. Primer Trimming: Use iVar trim to trim primer regions and then sort the trimmed bam by coordinates using Samtools
6. Variant Calling: Use iVar variants to call variants from the trimmed and sorted bam
   - iVar variants parameters include a minimum quality score set to 30, a minimum variant base frequency set to 0.6 and a minimum read depth set to 10.
7. Consensus Calling: Use iVar consensus to call the consensus genome sequence from the trimmed and sorted bam
   - iVar consensus parameters include a minimum quality score set to 30, a minimum variant base frequency set to 0.6 and a minimum read depth set to 10.
8. Coverage Statistics: Use Samtools flagstat, stats, and coverage to output statistics from the bam

```mermaid
---
title: viral_amp_illumina_pe_assembly
---
graph LR
    A[Raw Paired-End Reads] --> B[Human Read Scrubbing<br/>Hostile]
    B --> C[Quality Filtering<br/>SeqyClean]
    A & C --> D[Quality Assessment<br/>FastQC]
    C --> E[Read Alignment<br/>BWA]
    E --> F[Primer Trimming<br/>iVar]
    F --> G[Variant Calling<br/>iVar]
    F --> H[Coverage Stats<br/>Samtools]
    F --> I[Consensus Calling<br/>iVar]
    I --> J[FASTA Header Formatting]
    D & F & G & H & I & J --> K[Assembly Workflow Outputs]
    K --> L[Cloud Storage Transfer<br/>optional task]
```
### Inputs

#### 1. Terra Data Table

The terra data table must include the following columns as listed below. Columns noted as required, are required for running viral_amp_illumina_pe_assembly.wdl. Optional columns are not necessary for the assembly workflow but must be present if `transfer_results` is set to true or for running the SC2_lineage_calling_and results.wdl and Transfer workflows described below under `Lineage Calling Workflows` and `Transfer Workflows`, respectively.

| column header      | description                                                                                      |
| ------------------ | ------------------------------------------------------------------------------------------------ |
| `entity:sample_id` | column with the list of sample names. (e.g. `entity:test-measles-0005-miseq_id`)                             |
| `fastq_1`          | The google bucket path to the R1 fastq file.                                                            |
| `fastq_2`          | The google bucket path to the R2 fastq file.                                                            |
| `run_date`         | The date on which the sequencing run performed.                                                       |
| `project_name`     | The bane of your sequencing run (e.g. project name on BaseSpace).                                                      |
| `out_dir`          | (optional; out_dir is needed if `transfer_results` = true) User defined google bucket for where the files will be transferred.           |
| `workbook_path`    | (optional; required for summary workflow)                                                        |
| `project_name`     | (optional; required for summary workflow)                                                        |

#### 2. Terra Workspace Data

See [setup](setup.md).

#### 3. Setting Up the Workflow Inputs

Below is a summary of the workflow input variables along with the syntax used for the attribute column when setting up the workflow to run on Terra.bio:

| workflow variable           | attribute (input syntax into workflow)                         |
| --------------------------- | ---------------------------------------------------------------|
| `contam_fasta`              | workspace.adapters_and_contaminants                            |
| `analysis_date`             | this.run_date                                                  |
| `viral_amp_ref_fasta`       | workspace.viral_amp_*_ref_fasta                                |
| `viral_amp_ref_gff`         | workspace.viral_amp_*_ref_gff                                  |
| `fastq_1`                   | this.fastq_1                                                   |
| `fastq_2`                   | this.fastq_2                                                   |
| `out_dir`                   | this.out_dir                                                   |
| `overwrite`                 | `true` or `false`                                              |
| `viral_amp_primer_bed`      | workspace.viral_amp_*_primer_bed                               |
| `project_name`              | this.project_name                                              |
| `sample_name`               | this.{entity_name}\_id                                         |
| `transfer_results`          | `true` or `false`                                              |
| `scrub_reads`               | `true` or `false`                                              |
| `scrub_genome_index`        | workspace.hostile_human_t2t_hla_bt2 (if using read scrubbing)  |


### Outputs

| WDL task name            | software/program                           | variable name                | description                                                                                                             |
| ------------------------ | ------------------------------------------ | ---------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| hostile                  | hostile                                    | `human_reads_removed`      | integer                                                                                                                   |
| hostile                  | hostile                                    | `human_reads_removed_proportion` | floating-point number                                                                                               |
| hostile                  | hostile                                    | `fastq1_scrubbed`          | file
| hostile                  | hostile                                    | `fastq2_scrubbed`          | file
| seqyclean                | seqyclean                                  | `filtered_reads_1`         | file                                                                                                                    |
| seqyclean                | seqyclean                                  | `filtered_reads_2`         | file                                                                                                                    |
| seqyclean                | seqyclean                                  | `seqyclean_summary`        | file                                                                                                                    |
| fastqc as fastqc_raw     | fastqc                                     | `fastqc_raw1_html`         | file                                                                                                                    |
| fastqc as fastqc_raw     | fastqc                                     | `fastqc_raw1_zip`          | file                                                                                                                    |
| fastqc as fastqc_raw     | fastqc                                     | `fastqc_raw2_html`         | file                                                                                                                    |
| fastqc as fastqc_raw     | fastqc                                     | `fastqc_raw2_zip`          | file                                                                                                                    |
| fastqc as fastqc_cleaned | fastqc                                     | `fastqc_clean1_html`       | file                                                                                                                    |
| fastqc as fastqc_cleaned | fastqc                                     | `fastqc_clean1_zip`        | file                                                                                                                    |
| fastqc as fastqc_cleaned | fastqc                                     | `fastqc_clean2_html`       | file                                                                                                                    |
| fastqc as fastqc_cleaned | fastqc                                     | `fastqc_clean2_zip`        | file                                                                                                                    |
| align_reads              | bwa and samtools                           | `out_bam`                  | file                                                                                                                    |
| align_reads              | bwa and samtools                           | `out_bamindex`             | file                                                                                                                    |
| ivar trim                | ivar trim and samtools                     | `trim_bam`                 | file                                                                                                                    |
| ivar trim                | ivar trim and samtools                     | `trimsort_bam`             | file                                                                                                                    |
| ivar trim                | ivar trim and samtools                     | `trimsort_bamindex`        | file                                                                                                                    |
| ivar variants            | ivar variants                              | `variants`                 | vcf file formatted as a tsv                                                                                             |
| ivar consensus           | ivar consensus                             | `consensus`                | fasta file of consensus genome, Ns are called in places with less than 10 bp read depth.                                |
| bam_stats                | samtools flagstat, stats, percent_coverage | `flagstat_out`             | file                                                                                                                    |
| bam_stats                | samtools flagstat, stats, percent_coverage | `stats_out`                | file                                                                                                                    |
| bam_stats                | samtools flagstat, stats, percent_coverage | `covhist_out`              | file                                                                                                                    |
| bam_stats                | samtools flagstat, stats, percent_coverage | `cov_out`                  | file                                                                                                                    |                                                        
| rename_fasta             | N/A                                        | `renamed_consensus`        | fasta file; consensus genome sequence with the fasta header renamed to remove extraneous information added by iVar                        |
| version_capture          | version_capture.py                         | `version_capture_file` |  file                                                                                                    |
| transfer                 | gsutil                                     | `transfer_date_assembly`   | String                                                                                                                  |


## viral_amp_illumina_pe_summary

File: viral_amp_illumina_pe_summary.wdl

The workflow accepts "sample_set" as the root entity type and uses the data table from the assembly workflow. 
This workflow summarizes results from viral amplicon Illumina paired-end sequencing assembly runs. It processes outputs from the assembly workflow to generate consolidated reports and transfers results to Google Cloud Storage (if `transfer_results` = true). The workflow performs the following operations:

1. Concatenate Consensus Sequences - Combines all consensus sequences from the sequencing run into a single multi-FASTA file
2. Summarize Results - Aggregates sequencing metrics, coverage statistics, and metadata into a comprehensive CSV report
3. Transfer Outputs - Uploads summary files and concatenated sequences to a user-defined Google Cloud Storage bucket (if `transfer_results` = true)
4. Capture Version Information - Records workflow and tool versions for reproducibility


```mermaid
---
title: viral_amp_illumina_pe_summary workflow
---
graph LR
    A[Assembly Workflow Outputs] --> B[Sequence Concatenation]
    A --> C[Results Summarization]
    C --> D[Summary Report]
    B & D --> E[Cloud Storage Transfer<br/>optional task]
```

### Inputs

Below is a summary of the workflow input variables along with the syntax used for the attribute column when setting up the workflow to run on Terra.bio. For the attributes, the "this.sample{terra_data table_name}s." syntax tells Terra to pull the variable from the sample-level terra data table. These variables were either in the original terra data table as inputs for the assembly workflow or added as outputs during the assembly workflow (see reference based assembly workflow inputs and outputs sections for more details). The "workspace." syntax tells Terra to pull the variable from the terra workspace data. Workspace data is described in the `Getting Started` drop down menu above.

| workflow variable                                | terra attribute (input syntax into workflow)                          |     |
| ------------------------------------------------ | --------------------------------------------------------------------- | --- |
| `concat_seq_results_py`                          | workspace.viral_amp_concat_results_py                                 |     |
| `cov_out`                                        | this.sample{terra_data_table_name}s.cov_out                           |     |
| `out_dir_array`                                  | this.sample{terra_data_table_name}s.out_dir                           |     |
| `transfer_results`                               | `true` or `false`                                                     |     |
| `project_name_array`                             | this.sample{terra_data_table_name}s.project_name                      |     |
| `renamed_consensus`                              | this.sample{terra_data_table_name}s.renamed_consensus                 |     |
| `sample_name`                                    | this.sample{terra_data_table_name}s.sample{terra_data_table_name}\_id |     |
| `workbook_path_array`                            | this.sample{terra_data_table_name}s.workbook_path                     |     |
| `assembler_version_array`                        | this.sample{terra_data_table_name}s.assembler_version                 |     | 
| `workflow_version`                               | this.sample{terra_data_table_name}s.wf_version                        |     |
| `workflow_version_und`                           | this.sample{terra_data_table_name}s.wf_version_und                    |     |

### Outputs

This workflow generates several output files that can be transferred to the user defined user google bucket defined by this.sample{terra_data_table_name}s.out_dir. The table below details each output.

| output variable name                          | file_name                                                               | description                                    | google bucket path (if transferred)                  |
| --------------------------------------------- | ----------------------------------------------------------------------- | ---------------------------------------------- | ------------------------------------------------- |
| `cat_fastas`                                  | `concatenate_assemblies.fasta`                                          | all consensus sequences from assembly in a single fasta file  | `gs://{user_defined_gcp_bucket}/multifasta/`|
| `sequencing_results_csv`                      | `{seq_run}_sequencing_results.csv`                                      | summary of the sequencing metrics for each sequence generated from the `concat_seq_results.py` script. | `gs://{user_defined_gcp_bucket}/summary_results/` |
| `transfer_date`                               | N/A                                                                     | date transfer task was run (if `transfer_results` = true)                     | N/A                                               |
| `wf_version_output`                           | `version_capture_{workflow_name}_{project_name}_{workflow_version}.csv` | version capture CSV file                       | `gs://{user_defined_gcp_bucket}/summary_results/` |



## viral_amp_wwt_variant_calling

File: viral_amp_wwt_variant_calling.wdl

The workflow accepts "sample_set" as the root entity type and uses the data table from the assembly workflow. This workflow uses Freyja (https://andersen-lab.github.io/Freyja/index.html#) to determine which viral lineages are likely present in the sample and estimates abundances
The workflow performs the following operations:

1. Variant Calling: Use iVar variants to call variants from the trimmed and sorted bam from the assembly workflow with adjusted parameter to detect lower frequnecy variants needed for wastewater analysis.
    - iVar variants parameters include a minimum quality score set to 20, a minimum variant base frequency set to 0.01 and a minimum read depth set to 10.
2. Lineage deconvolution: Uses Freyja demix to determine which viral lineages are likely present in the sample and estimates abundances.
3. Mutations: takes the variant calling results and generates of file of all mutaions detected for each sample
4. Aggregate results: Uses Freyja aggregate to pull the Freyja demix outputs into a single file.
5. Combine mutations: takes the mutations files for all samples and aggregeates them into a single file
3. Transfer Outputs - (if `transfer_results` = true) Uploads summary files and concatenated sequences to a user-defined Google Cloud Storage bucket
4. Capture Version Information - Records workflow and tool versions for reproducibility

```mermaid
---
title: viral_amp_wwt_variant_calling
---
graph LR
    A[Assembly Workflow Outputs<br/>bam files] --> B[Variant Calling<br/>ivar]
    B --> C[Lineage Deconvolution<br/>Freyja]
    C --> D[Aggregate Results<br/>Freyja]
    B --> E[Generate Mutations Table]
    E --> F[Aggregate Mutations Table]
    D & F --> G[Version Capture]
    G --> H[Cloud Storage Transfer<br/>optional task]
```
### Inputs

Below is a summary of the workflow input variables along with the syntax used for the attribute column when setting up the workflow to run on Terra.bio.

| workflow variable                        | attribute (input syntax into workflow)                                |
| ---------------------------------------- | --------------------------------------------------------------------- |
| `freyja_pathogen`                        | eg: ["MEASLES"], must match Freyja's documentation                    |
| `reference_genome`                       | workspace.viral_amp_*_ref_fasta                                       |
| `reference_gff` (optional)               | workspace.viral_amp_*_ref_gff                                         |
| `out_dir_array`                          | this.{terra_data_table_names}s.out_dir                                |
| `overwrite`                              | `true` or `false`                                                     |
| `project_name_array`                     | `this.{terra_data_table_name}s.project_name                           |
| `sample_name`                            | this.{terra_data_table_name}s.{terra_data_table_name}_id              |
| `trimsort_bam`                           | this.{terra_data_table_name}s.trimsort_bam                            |
| `workflow_version`                       | this.{terra_data_table_name}s.wf_version                              |
| `workflow_version_und`                   | this.{terra_data_table_name}s.wf_version_und                          |
| `transfer_results`                       | `true` or `false`                                                     |

### Outputs

This workflow generates several output files which can if transferred to the user defined google bucket as defined by a string (e.g. "gs://covid_terra/NEXSEQ_101/terra_outputs"), if `transfer_results` = true. The table below details each output. 

| output variable name                  | file_name                                                               | description                                            | google bucket path                                                            |
| ------------------------------------- | ----------------------------------------------------------------------- | ------------------------------------------------------ | ----------------------------------------------------------------------------- |
| `combined_mutations_tsv`              | `combined_mutations.tsv`                                                | mutation counts for all samples                        | `gs://{user_defined_gcp_bucket}/viral_amp_wwt_variant_calling/`                 |
| `variants`                            | `{sample_name}_variants.tsv`                                            | generated for each sample; output from variant_calling task    | `gs://{user_defined_gcp_bucket}/viral_amp_wwt_variant_calling/freyja/`          |
| `depth`                               | `{sample_name}_depth.tsv`                                               | generated for each sample; output from variant_calling task | `gs://{user_defined_gcp_bucket}/viral_amp_wwt_variant_calling/freyja/`          |
| `demix`                               | `{sample_name}_demixed.tsv`                                             | generated for each sample; output from freyja demix    | `gs://{user_defined_gcp_bucket}/viral_amp_wwt_variant_calling/freyja/`          |
| `demix_aggregated`                    | `demix_aggregated.tsv`                                                  |  aggregation of demix outputs into a single file; output from freyja aggregate                                                      | `gs://{user_defined_gcp_bucket}/viral_amp_wwt_variant_calling/`                 |
| `transfer_date_wwt_variant_calling`   | N/A                                                                     | date the workflow was run                              | N/A                                                                           |
| `version_capture_wwt_variant_calling` | `version_capture_{workflow_name}_{project_name}_{workflow_version}.csv` | software versions used in analysis                     | `gs://{user_defined_gcp_bucket}/summary_results/`                             |
