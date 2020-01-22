#######################
## Decompress for seqtk
#######################
rule decompress:
    input: 
        get_fastq
    output:
        fq1 = temp(WORKING_DIR + "{sample}_R1.fq"),
        fq2 = temp(WORKING_DIR + "{sample}_R2.fq")
    message:"decompressing {wildcards.sample} reads"
    params:
        sample_name = "{sample}"
    run:
        if sample_is_single_end(params.sample_name):
            shell("gunzip {input}")
            shell("touch {output.fq2")
        else:
        	shell("gunzip {input[0]}")
        	shell("gunzip {input[1]}")

################
## Read trimming
################
rule fastp:
    input:
        fq1 = temp(WORKING_DIR + "{sample}_R1.fq"),
        fq2 = temp(WORKING_DIR + "{sample}_R2.fq")
    output:
        fq1 = WORKING_DIR + "{sample}_R1_trimmed.fq",
        fq2 = WORKING_DIR + "{sample}_R2_trimmed.fq",
        html = RESULT_DIR + "fastp/{sample}.html",
        json = temp(WORKING_DIR + "fastp/{sample}.html")
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sample_name = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    run:
        if sample_is_single_end(params.sample_name): # single end
            shell("fastp --thread {threads} --html {output.html} --json {output.json} --qualified_quality_phred {params.qualified_quality_phred}  --in1 {input} --out1 {output} 2>{log}")
            shell("touch {output.fq2}")
        else:  
            shell("fastp --thread {threads} --html {output.html} --json {output.json}  --qualified_quality_phred {params.qualified_quality_phred}  --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2} 2>{log}")  


########
# subset
########
rule subset:
    input: 
        fq1 = WORKING_DIR + "{sample}_R1_trimmed.fq",
        fq2 = WORKING_DIR + "{sample}_R2_trimmed.fq"
    output:
        fq1 = WORKING_DIR + "{sample}_R1_trimmed_{percentage}_percent.fq",
        fq2 = WORKING_DIR + "{sample}_R2_trimmed_{percentage}_percent.fq",
    message: "Extracting {wildcards.percentage} percent of the initial {wildcards.sample} reads"
    params:
        sample_name = "{sample}",
        percentage = "{percentage}"
    run:
        if sample_is_single_end(params.sample_name): # single end
            shell("seqtk sample {input.fq1} {params.percentage} > {output.fq1}")
            shell("touch {output.fq2}")
        else:  
            shell("seqtk sample -s100 {input.fq1} {params.percentage} > {output.fq1}")
            shell("seqtk sample -s100 {input.fq2} {params.percentage} > {output.fq2}")