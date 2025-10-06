# CDPHE_Viral_Amp_WDL Workflows

## Disclaimer
**Next generation sequencing and bioinformatic and genomic analysis at the Colorado Department of Public Health and Environment (CDPHE) is not CLIA validated at this time. These workflows and their outputs are not to be used for diagnostic purposes and should only be used for public health action and surveillance purposes. CDPHE is not responsible for the incorrect or inappropriate use of these workflows or their results.**

<br/>

## Overview

The following documentation describes the Colorado Department of Public Health and Environment's workflows for the assembly and analysis of whole genome sequencing data of Viral Amplicon on GCP's Terra.bio platform. Workflows are written in WDL and can be imported into a Terra.bio workspace through dockstore (see Setup section below: https://dockstore.org/).

Our Viral whole genome reference-based assembly workflows are highly adaptable and facilitate the assembly and analysis of tiled amplicon based sequencing data of Viral samples. The workflows can accommodate various amplicon primer schemes including Artic V3, Artic V4, Artic V4.1, Artic V5.3.2 and Measles (BRAZ/WHO), as well as different sequencing technology platforms including both Illumina and Oxford Nanopore Technology (ONT). 

<br/>

## Workflows

Below is a list of available and maintained workflows and a brief description of the workflow. A full description of each workflow can be found on each workflow's readme page. 

<br/>

|Workflow Name | Description |
|--------------|-------------|
| ``viral_amp_wwt_illumina_pe_assembly`` | Performs reference-based assembly of viral genomes from Illumina paired-end amplicon data. |
| ``viral_amp_wwt_illumina_pe_summary`` | Generates summary statistics and quality metrics from assembled viral genomes. |
| ``viral_amp_wwt_variant_calling`` | Uses Freyja to estimate relative lineage abundances and variant composition from mixed viral samples (e.g., wastewater). |

## Workflow

```mermaid
graph TD
    subgraph Assembly
    A1(Raw Reads) --> B1
    A1 --> C1
    B1[filter_reads_seqyclean] --> D1
    B1 --> L1
    C1[assess_quality_fastqc] --> L1
    D1[align_reads_bwa] --> E1
    D1 --> L1
    E1[trim_primers_ivar] --> F1
    E1 --> G1
    E1 --> H1
    F1[call_variants_ivar] --> L1
    G1[call_consensus_ivar] --> I1
    G1 --> L1
    H1[calc_bam_stats_samtools] --> L1
    I1[rename_fasta] --> J1
    I1 --> K1
    I1 --> L1
    J1[calc_percent_coverage] --> L1
    K1[call_clades_nextclade] --> L1
    L1([Assembly Files]) --> M1
    M1[transfer_outputs] --> N1
    N1{{Cloud Bucket}}
    end
```

```mermaid
graph TD
    subgraph Summary
    A2(Assembly Files) --> B2
    B2[concatenate_consensus] --> D2
    A2 --> C2
    C2[summarize_results] --> D2
    D2([Summary Files]) --> E2
    E2[transfer_outputs] --> F2
    F2{{Cloud Bucket}}
    end
```

```mermaid
graph TD
    subgraph Variant_Calling
    A3(Assembly Files) --> B3
    A3 --> C3
    B3[align_reads_bwa] --> D3
    B3 --> L3
    C3[trim_primers_ivar] --> D3
    D3[call_variants_freyja] --> E3
    D3 --> F3
    E3[aggregate_lineage_abundances] --> G3
    F3[generate_variant_tables] --> G3
    G3[version_capture_viral_amp_variant_calling] --> H3
    G3 --> L3
    H3([Variant Calling Files]) --> I3
    I3[transfer_outputs] --> J3
    J3{{Cloud Bucket}}
    end



## Process

<br/>

### Wastewater Viral sequence assembly and variant calling
Processing of wastewater viral sequencing data involves three coordinated workflows (Figure 1). The first, ``viral_amp_wwt_illumina_pe_assembly``, takes raw paired-end Illumina reads and performs quality control, contamination filtering, primer trimming, reference-guided assembly, variant calling, and consensus genome generation. Intermediate files and consensus sequences are then transferred to a designated Google Cloud bucket(GCP) for storage and downstream access.

Next, the ``viral_amp_wwt_illumina_pe_summary`` workflow aggregates the outputs from multiple samples, concatenates consensus sequences, and generates a comprehensive sequencing results report (including coverage statistics and clade assignments), while also organizing the outputs into a versioned results directory.

Finally, the ``viral_amp_wwt_variant_calling workflow`` applies Freyja to estimate relative lineage abundances in wastewater samples, accounting for the mixed nature of viral populations present. 

<br/>

Figure 1. High level overview of workflow process for clinical and wastewater Viral samples.


![Viral Amp high level overview workflow diagram](./docs/img/SC2_overview_workflow_diagram.png "High level overview of Viral Amp workflow")

<br/>


## Setup

<br/>

### Input Workflow from Dockstore
To use the workflow on the Terra platform, first you will need to import the workflow from Dockstore. All workflows can be found under our dockstore organization called CDPHE-bioinformatics.
1. Go to dockstore (https://dockstore.org/).
2. Along the top search bar click on Organizations and search for "CDPHE".
3. Select the workflow.
4. On the right hand side of the workflow description, select "Launch with Terra".
5. Select the Destination workspace and select "Import". 
6. The workflow will now be displayed as a card under your workflows tab in your Terra workspace. 

<br/>

### Workspace Data
Prior to running any of the workflows, you must set up the Terra workspace data with the correct reference files and custom python scripts. The reference files can be found in this repository in the ``data/workspace_data`` directory. Python scripts can be found in the ``scripts`` directory. Workspace variables are named using the following format ``{organism}_{description}_{file_type}``, except for the primer bed files which are named as ``{description}_{file_type}``. Reference files and python scripts should be copied from this repo into a GCP bucket. The GCP bucket path to the file will serve as the "value" when adding data to the terra workspace data table. 

To add data to the terra workspace data:
1. Navigate to the Data tab in your Terra workspace.
2. In the left hand list of data tables, under "Other Data" select "Workspace Data".
3. Click on the "+" button in the lower right hand corner of the workspace data table. 
4. Fill in the "Key" column with the workspace variable name, the "Value" column with GCP bucket path to the file and the "Description" column with a brief description if desired. 
5. Once complete hit the check mark to the right.

Below is a data table detailing the workspace data you will need to set up in order to run the Viral workflows. 

| workspace variable name | workflow | file name | description |
|---------------------------|-----------------|--------------------|-----------------|
| ``adapters_and_contaminants`` | ``viral_amp_wwt_illumina_pe_assembly`` | Adapters_plus_PhiX_174.fasta | Adapters and PhiX contaminant sequences removed during FASTQ cleaning and filtering using SeqyClean. Thanks to Erin Young at Utah Public Health Laboratory for providing this file! |
| ``calc_percent_coverage_py`` | ``viral_amp_wwt_illumina_pe_assembly`` | calc_percent_coverage.py | Python script used in the Viral Amplicon assembly workflow to calculate percent genome coverage from consensus sequences. |
| ``k2_standard_8gb`` | ``viral_amp_wwt_variant_calling`` | k2_standard_08gb_20230605.tar.gz | Kraken2 standard database (8 GB version) used for taxonomic classification during variant calling. |
| ``viral_amp_braz_primer_bed`` | ``viral_amp_wwt_illumina_pe_assembly`` | measles_braz_primers.bed | Primer BED file for the Measles Braz amplicon set used during genome assembly. |
| ``viral_amp_braz_ref_fasta`` | ``viral_amp_wwt_illumina_pe_assembly`` | AF266290_A_Zagreb_vax.fasta | Reference genome FASTA for Measles Braz amplicon assembly workflow. Matches the Zagreb vaccine strain used for assembly. |
| ``viral_amplicon_braz_ref_gff`` | ``viral_amp_wwt_illumina_pe_assembly`` | AF266290_A_Zagreb_vax.gff | Genome annotation file (GFF) corresponding to the Braz reference FASTA used for Measles viral amplicon assembly. |
| ``viral_amp_imap_primer`` | ``viral_amp_wwt_illumina_pe_assembly`` | measles_artic_v-1-0-0_primers.bed | Primer BED file for the Measles ARTIC v1.0.0 (IMAP) primer set used for assembly. |
| ``viral_amp_imap_ref_fasta`` | ``viral_amp_wwt_illumina_pe_assembly`` | NC_001498-1_measles_reference_genome.fasta | Reference genome FASTA for Measles IMAP amplicon assembly. Matches the whole genome reference used by Freyja. |
| ``viral_amp_imap_ref_gff`` | ``viral_amp_wwt_illumina_pe_assembly`` | NC_001498-1_measles_reference_genome.gff | Genome annotation file (GFF) for the Measles IMAP reference genome FASTA. |
| ``viral_amp_concat_seq_results_py`` | ``viral_amp_wwt_illumina_pe_assembly`` | concat_seq_results.py | Python script used to concatenate sequence metrics and results from the Viral Amplicon assembly workflow. |
| ``viral_amp_version_capture_py`` | ``viral_amp_wwt_illumina_pe_assembly`` | version_capture.py | Generates version capture output files for documenting software versions used in the assembly workflow. |
| ``viral_amp_version_capture_variant_calling_py`` | ``viral_amp_wwt_variant_calling`` | version_capture_viral_amp_variant_calling.py | Generates version capture output files for documenting software versions used in the variant calling workflow. |
