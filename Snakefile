"""
Snakefile for single-cell endoderm differentiation project

Author: Davis McCarthy
Affiliation: EMBL-EBI
Study: Single-cell endoderm differentiation project
Date: Sunday 12 June 2016
Run: snakemake --jobs 3000 --latency-wait 30 --cluster 'bsub -q short -R "rusage[mem=16000]" -M 16000 -o ./snake_logs -e ./snake_logs'
Latest modification:
  - todo

STUDY ID: 3963
STUDY and PROJECT TITLE: Single Cell RNAseq at various stages of HiPSCs differentiating toward definitive endoderm and endoderm derived lineages.
Study name abbreviation: SC-RNAseq-definitive endoderm
PROJECT ID: 2010
HMDMC number (ethical approval): 13/042, 15_074
Project Cost Code: S0901

STUDY ID: 4262
STUDY and PROJECT TITLE: Single cell RNA-seq and bisulfite-seq at various stares of HiPSCs differentiating toward definitive endoderm derived lineages
PROJECT ID: 2218
HMDMC number (ethical approval): 13/042, 15_074
Project Cost Code: S0901

Raw data from Sanger first downloaded to:
data_raw/{diff_expt}/run_{run}/cram

Raw data from Sanger sequencing core with some additions then uploaded to:
/nfs/research2/hipsci/drop/hip-drop/incoming/stegle

From there it is moved into the HipSci tracked data area by Ian Streeter:
/nfs/research2/hipsci/drop/hip-drop/tracked/

Transient files and analyses, and working files to share can be put in the scratch directory:
/nfs/research2/hipsci/drop/hip-drop/scratch/
"""

import glob

configfile: "config.yaml"

## path to fasta file for reference transcriptome
fasta = expand('{basedir}{fasta}', basedir=config['references_dir'], fasta=config['human_fasta'])
## path to kallisto index file
kallisto_idx = expand('{basedir}{index}', basedir=config['references_dir'], index=config['kallisto_idx'])

## define commands
cramtools_cmd = 'java -jar /nfs/software/stegle/cramtools-3.0.jar'
kallisto_cmd = config['kallisto_cmd']
read_kallisto_to_scesets_cmd = 'src/preprocessing/read_kallisto_to_scesets.R'
kallisto_idx = os.path.join(config['references_dir'], config['kallisto_idx'])

## parameter objects and samples
DIFF_EXPT = ['diff_1']
RUNS = ['run_19776']
SAMPLES = glob.glob("data_raw/*/*/cram/*.cram")
SAMPLES = [os.path.basename(w).replace('.cram', '') for w in SAMPLES]

## targets
bam_files = expand('data_raw/{diff_expt}/{run}/bam/{sample}.bam', diff_expt=DIFF_EXPT, run=RUNS, sample=SAMPLES)
kallisto_results = expand('data_raw/{diff_expt}/{run}/quant_kallisto/{sample}/abundance.tsv', diff_expt=DIFF_EXPT, run=RUNS, sample=SAMPLES)
fastqc_html_reports = expand('reports/fastqc/{diff_expt}/{run}/{sample}_fastqc.html', diff_expt=DIFF_EXPT, run=RUNS, sample=SAMPLES)
scesets = expand('data_processed/{diff_expt}/scesets_{run}_kallisto_preqc.RData', diff_expt=DIFF_EXPT, run=RUNS)
scater_first_html_reports = expand('reports/first_qc/{diff_expt}_{run}_first_qc.html', diff_expt=DIFF_EXPT, run=RUNS)

rule all:
    input:
        bam_files, kallisto_idx, kallisto_results,
        scesets, scater_first_html_reports, fastqc_html_reports,

rule build_kallisto_index:
    input:
        fasta
    output:
        kallisto_idx
    shell:
        '{kallisto_cmd} index -i {output} {input}'


rule cram2fastq:
    input:
        "data_raw/{diff_expt}/{run}/cram/{sample}.cram"
    output:
        "data_raw/{diff_expt}/{run}/fastq/{sample}_1.fastq",
        "data_raw/{diff_expt}/{run}/fastq/{sample}_2.fastq"
    params:
        prefix="data_raw/{diff_expt}/{run}/fastq/{sample}"
    shell:
        '{cramtools_cmd} fastq --input-cram-file {input} '
        '--fastq-base-name {params.prefix} '

rule cram2bam:
    input:
        "data_raw/{diff_expt}/{run}/cram/{sample}.cram"
    output:
        "data_raw/{diff_expt}/{run}/bam/{sample}.bam"
    shell:
        '{cramtools_cmd} bam --input-cram-file {input} '
        '--output-bam-file {output} '


rule fastqc_reports:
    input:
        bam_files
    output:
        'reports/fastqc/{diff_expt}/{run}/{sample}_fastqc.html'
    params:
        output_dir="reports/fastqc/{diff_expt}/{run}/"
    shell:
        'fastqc -o {params.output_dir} {input}'


rule kallisto_quant:
    input:
        kidx=kallisto_idx,
        fq1='data_raw/{diff_expt}/{run}/fastq/{sample}_1.fastq',
        fq2='data_raw/{diff_expt}/{run}/fastq/{sample}_2.fastq'
    output:
        'data_raw/{diff_expt}/{run}/quant_kallisto/{sample}/abundance.tsv'
    params:
        folder="data_raw/{diff_expt}/{run}/quant_kallisto/{sample}/"
    shell:
        '{kallisto_cmd} quant -i {input.kidx} '
        '-o {params.folder} '
        '--bias {input.fq1} {input.fq2}'


rule kallisto_to_sceset:
    input:
        kallisto_results
    output:
        'data_processed/{diff_expt}/scesets_{run}_kallisto_preqc.RData',
        'data_processed/{diff_expt}/scesets_{run}_kallisto_preqc.feather'
    params:
        input_dir='data_raw/{diff_expt}/{run}/quant_kallisto',
        output_dir='data_processed/{diff_expt}/ '
    shell:
        'Rscript {read_kallisto_to_scesets_cmd} '
        '--input_dir {params.input_dir} '
        '--output_dir {params.output_dir} '
        '--scesets_out scesets_{wildcards.run}_kallisto_preqc.RData'


rule rough_qc:
    input:
        'data_processed/{diff_expt}/scesets_{run}_kallisto_preqc.RData'
    output:
        'reports/first_qc/{diff_expt}_{run}_first_qc.html'
    shell:
        'Rscript src/R/compile_report.R -i {input} -o {output} '
        '--template src/Rmd/rough_qc_template.Rmd '
