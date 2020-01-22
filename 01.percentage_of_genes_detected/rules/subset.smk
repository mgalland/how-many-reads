

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


########
# subset
########
rule subset:
    input: 
        fq1 = WORKING_DIR + "{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "{sample}_R2_trimmed.fq.gz"
    output:
        fq1 = WORKING_DIR + "{sample}_R1_trimmed_{percentage}_percent.fq",
        fq2 = WORKING_DIR + "{sample}_R2_trimmed_{percentage}_percent.fq",
    message: "Extracting {wildcards.percentage} percent of the initial {wildcards.sample} reads"
    params:
        sample_name = "{sample}",
        tar_directory = WORKING_DIR + "tar/",
        sample_basename = get_name_without_file_extension("{input.fq1}"),
        params = "{percentage}"
    run:
        if sample_is_single_end(params.sample_name): # single end
            shell("tar -zxvf {input.fq1} --directory {params.tar_directory}")
            shell("seqtk sample {params.tar_directory}{sample_name} {percentage} > {output.fq1}")
            shell("touch {output.fq2}")
        else:  
            shell("fastp --thread {threads} --html {output.html} --qualified_quality_phred {params.qualified_quality_phred}  --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2} 2>{log}")  