rule fastp:
    input:
        R2="fastq/{sequencing_name}/demux/{sample_name}_R2.fastq.gz",
    output:
        R2="fastp/{sequencing_name}/{sample_name}_R2_fastq.gz",
    log:
        "logs/fastp/fastp_{sequencing_name}_{sample_name}.log",
    resources:
        mem_mb=1024,
    threads: 8
    message:
        "Running fastp on the demultiplexed samples."
    shell:
        """
        #fastp -I {input.R2} -O {output.R2} --threads {threads} &> {log}
        mkdir -p fastp/{wildcards.sequencing_name}
        touch {output.R2}
        """
