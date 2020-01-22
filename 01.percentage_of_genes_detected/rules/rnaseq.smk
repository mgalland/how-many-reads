


###############################
## Rules
###############################

#######################
## create sleuth object
#######################

rule run_sleuth:
    input:
        expand(RESULT_DIR + "kallisto/{samples}/abundance.tsv",samples=SAMPLES)
    output:
        normalized_counts = RESULT_DIR + "abundance_tidy.tsv"
    params:
        sample_file              = config["samples"],
        input_directory          = RESULT_DIR + "kallisto",
        output_directory         = RESULT_DIR
    threads: 10
    shell:
        "Rscript --vanilla scripts/sleuth_analysis.R -i {params.input_directory} "
        "-s {params.sample_file} "
        "-c {threads} "
        "-o {params.output_directory} "


##############################################################################
## Kallisto (pseudo-alignment) analysis for transcriptome and custom databases
##############################################################################

rule estimate_transcript_abundance_using_kallisto:
    input:
        index = WORKING_DIR + "index/kallisto_index.kidx",
        #reads = get_trimmed, 
        fq1 = WORKING_DIR + "{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "{sample}_R2_trimmed.fq.gz"
    output:
        RESULT_DIR + "kallisto/{sample}/abundance.tsv"
    message:"computing {wildcards.sample} abundances using kallisto"
    threads: 10
    params:
        sampleName      = "{sample}",
        outDir          = "results/kallisto/{sample}/",
        fragmentLength  = str(config["kallisto"]["fragment-length"]),
        sd              = str(config["kallisto"]["sd"]),
        bootstrap       = str(config["kallisto"]["bootstrap"])
    run:
        if sample_is_single_end(wildcards.sample):
            shell("mkdir -p results/kallisto/; \
            kallisto quant -i {input.index} -o {params.outDir} \
            --single -l {params.fragmentLength} -s {params.sd} \
            -b {params.bootstrap} \
            --threads {threads} \
            {input.fq1}")
        else:
            shell("mkdir -p results/kallisto/; \
            kallisto quant -i {input.index} -o {params.outDir} \
            -b {params.bootstrap} \
            --threads {threads} \
            {input.fq1} {input.fq2}")


rule create_kallisto_index:
    input:
        fasta = config["kallisto"]["fasta_ref"]
    output:
         WORKING_DIR + "index/kallisto_index.kidx"
    params:
        WORKING_DIR + "index/kallisto_index.kidx"
    message:"creating kallisto index"
    shell:
        "kallisto index --make-unique -i {params} {input};"

################
## Read trimming
################
rule fastp:
    input:
        get_fastq
    output:
        fq1 = WORKING_DIR + "{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html"
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sample_name = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    run:
        if sample_is_single_end(params.sample_name): # single end
            shell("fastp --thread {threads} --html {output.html} --qualified_quality_phred {params.qualified_quality_phred}  --in1 {input} --out1 {output} 2>{log}")
            shell("touch {output.fq2}")
        else:  
            shell("fastp --thread {threads} --html {output.html} --qualified_quality_phred {params.qualified_quality_phred}  --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2} 2>{log}")  




