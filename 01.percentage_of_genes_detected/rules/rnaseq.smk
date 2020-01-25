

##############################################################################
## Kallisto (pseudo-alignment) analysis for transcriptome and custom databases
##############################################################################

# rule estimate_transcript_abundance_using_kallisto:
#     input:
#         index = WORKING_DIR + "index/kallisto_index.kidx", 
#         fq1 = WORKING_DIR + "{sample}_R1_trimmed_{percentage}_percent.fq",
#         fq2 = WORKING_DIR + "{sample}_R2_trimmed_{percentage}_percent.fq"
#     output:
#         RESULT_DIR + "kallisto/{sample}_abundance_{percentage}_percent.tsv"
#     message:"computing {wildcards.sample} abundances using kallisto"
#     threads: 10
#     params:
#         sampleName      = "{sample}",
#         outDir          = "results/kallisto/{sample}/",
#         fragmentLength  = str(config["kallisto"]["fragment-length"]),
#         sd              = str(config["kallisto"]["sd"]),
#         bootstrap       = str(config["kallisto"]["bootstrap"])
#     run:
#         if sample_is_single_end(wildcards.sample):
#             shell("mkdir -p results/kallisto/; \
#             kallisto quant -i {input.index} -o {params.outDir} \
#             --single -l {params.fragmentLength} -s {params.sd} \
#             -b {params.bootstrap} \
#             --threads {threads} \
#             {input.fq1}")
#         else:
#             shell("mkdir -p results/kallisto/; \
#             kallisto quant -i {input.index} -o {params.outDir} \
#             -b {params.bootstrap} \
#             --threads {threads} \
#             {input.fq1} {input.fq2}")




#########################
# RNA-Seq read alignement
#########################

rule index:
    input:
        fasta = config["refs"]["genome"]
    output:
        [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing {input} fasta file"
    params:
        WORKING_DIR + "genome/genome"
    threads: 10
    shell:
        "hisat2-build -p {threads} {input} {params} --quiet"

rule hisat_mapping:
    input:
        fq1 = WORKING_DIR + "{sample}_R1_trimmed_{percentage}_percent.fq",
        fq2 = WORKING_DIR + "{sample}_R2_trimmed_{percentage}_percent.fq",
        index = [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    output:
        bams  = temp(WORKING_DIR + "mapped/{sample}_{percentage}_percent.bam"),
        sum   = RESULT_DIR + "logs/{sample}_{percentage}_percent_sum.txt",
        met   = RESULT_DIR + "logs/{sample}_{percentage}_met.txt"
    params:
        indexName = WORKING_DIR + "genome/genome",
        sampleName = "{sample}"
    message:
        "mapping {wildcards.sample} reads (subset {wildcards.percentage}%) to genome."
    threads: 10
    run:
        if sample_is_single_end(params.sampleName):
            shell("hisat2 -p {threads} --summary-file {output.sum} \
                --met-file {output.met} -x {params.indexName} \
                -U {input.fq1} | samtools view -Sb -F 4 -o {output.bams}")
        else:
            shell("hisat2 -p {threads} --summary-file {output.sum} \
                --met-file {output.met} -x {params.indexName} \
                -1 {input.fq1} -2 {input.fq2} | samtools view -Sb -F 4 -o {output.bams}")



#########################################
# Get table containing the raw counts
#########################################

rule create_counts_table:
    input:
        bams = expand(WORKING_DIR + "mapped/{sample}_{percentage}_percent.bam", sample = SAMPLES, percentage = PERCENTAGES),
        gff  = config["refs"]["gff"]
    output:
        RESULT_DIR + "counts.txt"
    message:"Creating count results "
    conda:
        "envs/subread.yaml"
    shell:
        "featureCounts -O -t mRNA -g ID -F 'gtf' -a {input.gff} -o {output} {input.bams}"
