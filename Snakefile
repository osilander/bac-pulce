##############################################################
### The input for this Snakefile is a fastq file of reads from
### a bac-pulce sequencing run and a genome to map to
### There are 5 basic steps:
### (1) map reads to focal genome(s)
### (2) get sorted .bam
### (3) get depth over focal organism
### (4) output file with total reads, totals mapped to each organism
### (5) output pdf of histogram of mapped reads of the high coverage areas
### 
### October 5 2019 Olin Silander
###############################################################

#shell.executable("/bin/bash")

BC, = glob_wildcards("./data/{barcode}.fastq")
#GENOMES, = glob_wildcards("/data/genome_data/flye/{genome}.fasta")
GENOMES = ["P_fluorescens_SBW25"]
#GENOMES = ["CCip28-2019-12-20_avon_flye","CCip20-2019-12-20_avon_flye","L4Cip1-2019-12-20_avon_flye","L3Cip2-2019-12-20_avon_flye","L3Cip2-circul-2019-12-20_avon_flye","CCip22-2019-12-20_avon_flye","DChl3-2019-12-20_avon_flye","DTBX1-circul-2019-12-20_avon_flye","DTBX1-2019-12-20_avon_flye","CChl2-2019-12-20_avon_flye"]

rule all:
    input:
        expand("results/{barcode}.{genome}.depth.pdf", barcode=BC, genome=GENOMES),
        expand("results/{barcode}.{genome}.stats.txt", barcode=BC, genome=GENOMES),
        expand("results/{barcode}.{genome}.cut_sites.txt", barcode=BC, genome=GENOMES),
        expand("results/{barcode}.{genome}.paf", barcode=BC, genome=GENOMES),
        expand("results/{barcode}.{genome}.rev.depth.txt", barcode=BC, genome=GENOMES),
        expand("results/{barcode}.{genome}.fwd.depth.txt", barcode=BC, genome=GENOMES),
        #expand("results/{barcode}.{genome}.cumulative_cut_sites.pdf", barcode=BC, genome=GENOMES),	

rule map_reads:
    input:
        reads="data/{barcode}.fastq",
        genome="/data/genome_data/{genome}.fasta"
    output: "results/{barcode}.{genome}.mapped.sorted.bam"
    benchmark: "results/benchmarks/{barcode}.{genome}.map.bench.txt"
    log: "results/logs/{barcode}.{genome}.map.log"
    #conda: "envs/minimap2.yaml"
    shell: "minimap2 -ax map-ont --secondary=no {input.genome} {input.reads} | samtools sort -o {output} - 2> {log}"

rule map_paf:
    input:
        reads="data/{barcode}.fastq",
        genome="/data/genome_data/{genome}.fasta"
    output: "results/{barcode}.{genome}.paf"
    benchmark: "results/benchmarks/{barcode}.{genome}.paf.bench.txt"
    log: "results/logs/{barcode}.{genome}.paf.log"
    #conda: "envs/minimap2.yaml"
    shell: "minimap2 -x map-ont --secondary=no {input.genome} {input.reads} > {output} 2> {log}"

rule get_depth:
    input: "results/{barcode}.{genome}.mapped.sorted.bam"
    output: "results/{barcode}.{genome}.depth.txt"
    #params:
    benchmark: "results/benchmarks/{barcode}.{genome}.depth.bench.txt"
    log: "results/logs/{barcode}.{genome}.depth.log"
    #conda: "envs/samtools.yaml"
    shell: "samtools depth {input} > {output} 2> {log}"

rule get_rev_depth:
    input: "results/{barcode}.{genome}.mapped.sorted.bam"
    output: "results/{barcode}.{genome}.rev.depth.txt"
    #conda: "envs/samtools.yaml"
    shell: "samtools view -bh -f 16 {input} | samtools depth - > {output}"

rule get_fwd_depth:
    input: "results/{barcode}.{genome}.mapped.sorted.bam"
    output: "results/{barcode}.{genome}.fwd.depth.txt"
    #conda: "envs/samtools.yaml"
    shell: "samtools view -bh -f 0 {input} | samtools view -bh -F 16 - | samtools depth - > {output}"

rule plot_depth:
    input: "results/{barcode}.{genome}.depth.txt"
    output: "results/{barcode}.{genome}.depth.pdf"
    shell: "R --slave --no-restore --file=scripts/plot_depth_window.R --args {input} {output}"

rule get_cuts:
    input: "results/{barcode}.{genome}.mapped.sorted.bam"
    output: "results/{barcode}.{genome}.cut_sites.txt"
    #params:
    benchmark: "results/benchmarks/{barcode}.{genome}.cuts.bench.txt"
    log: "results/logs/{barcode}.{genome}.cuts.log"
    #conda: "envs/samtools.yaml"
    shell: "samtools view {input} | cut -f1,2,3,4 > {output} 2> {log}"

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
