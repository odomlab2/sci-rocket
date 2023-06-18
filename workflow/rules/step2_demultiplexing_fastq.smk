# Rules related to demultiplexing sci-RNA-seq3 data.


# ---- Split R1 and R2 files into smaller files which will be handled in parallel. ----
rule split_R1:
    input:
        "fastq/{sequencing_name}/raw/Undetermined_S0_L001_R1_001.fastq.gz",
    output:
        temp(
            scatter.fastq_split(
                "fastq/{{sequencing_name}}/raw/R1_{scatteritem}.fastq.gz"
            )
        ),
    resources:
        mem_mb=512,
    params:
        out=lambda w: [
            f"-o fastq/{w.sequencing_name}/raw/R1_{i}-of-"
            + str(workflow._scatter["fastq_split"])
            + ".fastq.gz"
            for i in range(1, workflow._scatter["fastq_split"] + 1)
        ],
    message:
        "Generating splitted R1 files."
    shell:
        """
        fastqsplitter -i {input} {params.out} -t 1
        """


rule split_R2:
    input:
        "fastq/{sequencing_name}/raw/Undetermined_S0_L001_R2_001.fastq.gz",
    output:
        temp(
            scatter.fastq_split(
                "fastq/{{sequencing_name}}/raw/R2_{scatteritem}.fastq.gz"
            )
        ),
    resources:
        mem_mb=512,
    params:
        out=lambda w: [
            f"-o fastq/{w.sequencing_name}/raw/R2_{i}-of-"
            + str(workflow._scatter["fastq_split"])
            + ".fastq.gz"
            for i in range(1, workflow._scatter["fastq_split"] + 1)
        ],
    message:
        "Generating splitted R2 files."
    shell:
        """
        fastqsplitter -i {input} {params.out} -t 1
        """


# ---- Demultiplex each splitted R1 and R2 file to generate sample-specific fastq.gz files. ----
rule demultiplex_fastq_split:
    input:
        R1="fastq/{sequencing_name}/raw/R1_{scatteritem}.fastq.gz",
        R2="fastq/{sequencing_name}/raw/R2_{scatteritem}.fastq.gz",
    output:
        temp(directory("fastq/{sequencing_name}/demux_scatter/{scatteritem}/")),
    log:
        "logs/demux/demultiplex_fastq_{sequencing_name}_{scatteritem}.log",
    resources:
        mem_mb=1024,
    params:
        path_samples=config["path_samples"],
        path_barcodes=config["path_barcodes"],
    message:
        "Demultiplexing the scattered .fastq.gz files."
    shell:
        "python3.10 {workflow.basedir}/scripts/demultiplexing_samples.py --sequencing_name {wildcards.sequencing_name} --samples {params.path_samples} --barcodes {params.path_barcodes} --r1 {input[0]} --r2 {input[1]} --out {output} &> {log}"


rule gather_demultiplex_fastq_split:
    input:
        gather.fastq_split("fastq/{{sequencing_name}}/demux_scatter/{scatteritem}/"),
    output:
        R1_discarded="fastq/{sequencing_name}/demux/{sequencing_name}_R1_discarded.fastq.gz",
        R2_discarded="fastq/{sequencing_name}/demux/{sequencing_name}_R2_discarded.fastq.gz",
    message:
        "Combining the scattered demultiplexed results."
    shell:
        """
        touch {output.R1_discarded} {output.R2_discarded}
        """


rule gather_combined_demultiplexed_samples:
    input:
        R1="fastq/{sequencing_name}/demux/{sequencing_name}_R1_discarded.fastq.gz",
        R2="fastq/{sequencing_name}/demux/{sequencing_name}_R2_discarded.fastq.gz",
    output:
        R1=dynamic("fastq/{sequencing_name}/demux/{sample_name}_R1.fastq.gz"),
        R2=dynamic("fastq/{sequencing_name}/demux/{sample_name}_R2.fastq.gz"),
    message:
        "This should be a command for each sample, same discarded R1 and R2 files."
    shell:
        """
        {input.R1} {input.R2} {output.R1} {output.R2}
        """
