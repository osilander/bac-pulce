##############################################################
### The input for this Snakefile is a fasta file of reads from
### a sequencing run (likely barcoded)
### There are X steps:
### (1) map reads to genome from non focal organism
### (2) get sorted .bam
### (3) map reads to focal genome
### (4) get sorted .bam
### (5) get depth over focal organism
### (6) output file with total reads, totals mapped to each organism
### (7) output pdf of histogram of mapped reads
### 
### October 5 2019 Olin Silander
###############################################################

#shell.executable("/bin/bash")

BC, = glob_wildcards("./data/{barcode}.fastq")
GENOMES = ["Ecoli_MG1655"]

rule all:
    input:
        #expand("examples/results/{barcode}.{genome}.depth.pdf", barcode=BC, genome=GENOMES),
        expand("results/{barcode}.{genome}.stats.txt", barcode=BC, genome=GENOMES),
        expand("results/{barcode}.{genome}.cut_sites.pdf", barcode=BC, genome=GENOMES),
        expand("results/{barcode}.{genome}.cumulative_cut_sites.pdf", barcode=BC, genome=GENOMES),	

rule map_reads:
    input:
        reads="data/{barcode}.fastq",
        genome="/data/genome_data/{genome}.fasta"
    output: "results/{barcode}.{genome}.mapped.sorted.bam"
    benchmark: "results/benchmarks/{barcode}.{genome}.map.bench.txt"
    log: "results/logs/{barcode}.{genome}.map.log"
    #conda: "envs/minimap2.yaml"
    shell: "minimap2 -ax map-ont {input.genome} {input.reads} | samtools sort -o {output} - 2> {log}"

rule get_depth:
    input: "results/{barcode}.{genome}.mapped.sorted.bam"
    output: "results/{barcode}.{genome}.depth.txt"
    #params:
    benchmark: "results/benchmarks/{barcode}.{genome}.depth.bench.txt"
    log: "results/logs/{barcode}.{genome}.depth.log"
    #conda: "envs/samtools.yaml"
    shell: "samtools depth {input} > {output} 2> {log}"

rule plot_depth:
    input: "results/{barcode}.{genome}.depth.txt"
    output: "results/{barcode}.{genome}.depth.pdf"
    shell: "R --slave --no-restore --file=scripts/plot_depth.R --args {input} {output}"

rule get_cuts:
    input: "results/{barcode}.{genome}.mapped.sorted.bam"
    output: "results/{barcode}.{genome}.cut_sites.txt"
    #params:
    benchmark: "results/benchmarks/{barcode}.{genome}.cuts.bench.txt"
    log: "results/logs/{barcode}.{genome}.cuts.log"
    #conda: "envs/samtools.yaml"
    shell: "samtools view {input} | cut -f4 > {output} 2> {log}"

rule plot_cut_sites:
    input: "results/{barcode}.{genome}.cut_sites.txt"
    output: "results/{barcode}.{genome}.cut_sites.pdf"
    shell: "R --slave --no-restore --file=scripts/hist_cuts.R --args {input} {output}"

rule plot_cumulative:
    input: "results/{barcode}.{genome}.cut_sites.txt"
    output: "results/{barcode}.{genome}.cumulative_cut_sites.pdf"
    shell: "R --slave --no-restore --file=scripts/cumulative_cuts.R --args {input} {output}"

rule get_stats:
    input: "results/{barcode}.{genome}.mapped.sorted.bam"
    output: "results/{barcode}.{genome}.stats.txt"
    #params:
    benchmark: "results/benchmarks/{barcode}.{genome}.stats.bench.txt"
    log: "results/logs/{barcode}.{genome}.stats.log"
    #conda: "envs/samtools.yaml"
    shell: "samtools flagstat {input} > {output} 2> {log}"
