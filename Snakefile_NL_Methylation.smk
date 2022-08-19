# remember to make logs_slurm directory
#"conda create --name environment.yml" to create new enviroment
##"conda activate "environment name"" to activate environment
#Snakemake -s Snakefile_NL_Methylation.smk --profile slurm
##By default, searches Home/.config/snakemake/ for /slurm for config.yaml file


#Create sample list
SAMPLES, = glob_wildcards("combined_fastq/{sample}.fastq.gz")

fast5_dir = "/home/derick.k.thompson/NL_Methylation/fast5_pass"
rn6_mmi = "/home/groupdirs/genomic_resources/genomic_resources/UCSC/Rat/minimap_index_rn6/rn6.mmi"
rn6_fa = "/home/groupdirs/genomic_resources/genomic_resources/UCSC/Rat/minimap_index_rn6/rn6.fa"
summary_fofn = "/home/derick.k.thompson/NL_Methylation/sequencing_summary.fofn"

#Define final target files
rule all:
    input:
        expand("meth_calls/{sample}_meth_calls.tsv", sample = SAMPLES),
        expand("combined_fastq/{sample}.fastq.gz.index.fai", sample = SAMPLES),
        expand("sorted_bams/{sample}_sorted.bam.bai", sample = SAMPLES),
        expand("mapping_stats/{sample}_flagstat.txt", sample = SAMPLES),
        expand("methylation_frequency/{sample}_methylation_frequency.tsv", sample=SAMPLES)
        

#Index fast5/fastq
rule nano_index:
    input:
        fq = "combined_fastq/{sample}.fastq.gz", 
        summary = "seq_summary/sequencing_summary_{sample}.txt"
    output:
        fai = ("combined_fastq/{sample}.fastq.gz.index.fai"),
        readdb = ("combined_fastq/{sample}.fastq.gz.index.readdb")
    resources: time_min=40320, mem_mb=190000, cpus=24
    params: partition = "talon"
    benchmark:
        "benchmarks/{sample}_nano_index.benchmark.txt"
    shell:
        "nanopolish index -d {fast5_dir} -s {input.summary} {input.fq}"
#Mapping reads, pipe to samtools
rule minimap:
    input:
        fq = ("combined_fastq/{sample}.fastq.gz"),
        fai = ("combined_fastq/{sample}.fastq.gz.index.fai")
    output:
        bam = ("sorted_bams/{sample}_sorted.bam")
    resources: time_min=40320, mem_mb=190000, cpus=24
    params: partition = "talon"
    benchmark:
        "benchmarks/{sample}_minimap.benchmark.txt"
    shell:
        "minimap2 -t {resources.cpus} -a -x map-ont {rn6_mmi} {input.fq} | samtools view -@ {resources.cpus} -S -b | samtools sort -@ {resources.cpus} -o {output.bam}"
#Minimap options: -a output to SAM; -x map-ont is default, sets error rate for noisy long reads. 
#Samtools options: -S if input is SAM; -b output to BAM


rule sam_index:
    input:
        "sorted_bams/{sample}_sorted.bam"
    output:
        "sorted_bams/{sample}_sorted.bam.bai"
    resources: time_min=2880, mem_mb=190000, cpus=12
    params: partition = "talon"
    benchmark:
        "benchmarks/{sample}_sam_index.benchmark.txt"
    shell:
        "samtools index -@ {resources.cpus} {input} {output}"

rule flagstat:
    input:
        "sorted_bams/{sample}_sorted.bam"
    output:
        "mapping_stats/{sample}_flagstat.txt"
    resources: time_min=2880, mem_mb=50000, cpus=12
    params: partition = "talon"
    benchmark:
        "benchmarks/{sample}_flagstat.benchmark.txt"
    shell:
        "samtools flagstat -@ {resources.cpus} {input} > {output}"

rule nano_meth_calling:
    input:
        fq = ("combined_fastq/{sample}.fastq.gz"),
        bam = ("sorted_bams/{sample}_sorted.bam"),
        bai = ("sorted_bams/{sample}_sorted.bam.bai"),
        readdb = ("combined_fastq/{sample}.fastq.gz.index.readdb")
    output:
        "meth_calls/{sample}_meth_calls.tsv"
    resources: time_min=40320, mem_mb=190000, cpus=24
    params: partition = "talon"
    benchmark:
        "benchmarks/{sample}_nano_meth_calling.benchmark.txt"
    shell:
        "nanopolish call-methylation -t {resources.cpus} -r {input.fq} -b {input.bam} -g {rn6_fa} > {output}"

#Meth freq using nanopolish script
rule call_freq:
    input:
        rules.nano_meth_calling.output
    output:
        "methylation_frequency/{sample}_methylation_frequency.tsv"
    resources: time_min=2880, mem_mb=100000, cpus=8
    params: partition = "talon"
    benchmark:
        "benchmarks/{sample}_call_freq.benchmark.txt"
    shell:
        "calculate_methylation_frequency.py -s {input} > {output}"
