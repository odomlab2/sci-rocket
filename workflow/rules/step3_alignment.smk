rule fastp:
    input:
        R2="{sequencing_name}/demux_reads/{sample_name}_R2.fastq.gz",
    output:
        R2="{sequencing_name}/demux_reads/fastp_{sample_name}_R2.fastq.gz",
    log:
        "logs/step3_alignment/fastp_{sequencing_name}_{sample_name}.log",
    resources:
        mem_mb=1024,
    threads: 8
    message:
        "Running fastp on the demultiplexed samples."
    shell:
        """
        #fastp -I {input.R2} -O {output.R2} --threads {threads} &> {log}
        touch {output.R2}
        """
