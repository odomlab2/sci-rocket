# Rules related to demultiplexing sci-RNA-seq3 data.


# ---- Split R1 and R2 files into smaller files which will be handled in parallel. ----
rule split_R1:
    input:
        "{sequencing_name}/raw_reads/Undetermined_S0_R1_001.fastq.gz",
    output:
        temp(
            scatter.fastq_split(
                "{{sequencing_name}}/raw_reads/R1_{scatteritem}.fastq.gz"
            )
        ),
    threads: 5
    resources:
        mem_mb=1024 * 2,
    params:
        out=lambda w: [
            f"-o {w.sequencing_name}/raw_reads/R1_{i}-of-"
            + str(workflow._scatter["fastq_split"])
            + ".fastq.gz"
            for i in range(1, workflow._scatter["fastq_split"] + 1)
        ],
    message:
        "Generating multiple evenly-sized R1 chunks."
    shell:
        """
        fastqsplitter -i {input} {params.out} -t 1
        """


rule split_R2:
    input:
        "{sequencing_name}/raw_reads/Undetermined_S0_R2_001.fastq.gz",
    output:
        temp(
            scatter.fastq_split(
                "{{sequencing_name}}/raw_reads/R2_{scatteritem}.fastq.gz"
            )
        ),
    threads: 5
    resources:
        mem_mb=1024 * 2,
    params:
        out=lambda w: [
            f"-o {w.sequencing_name}/raw_reads/R2_{i}-of-"
            + str(workflow._scatter["fastq_split"])
            + ".fastq.gz"
            for i in range(1, workflow._scatter["fastq_split"] + 1)
        ],
    message:
        "Generating multiple evenly-sized R2 chunks."
    shell:
        """
        fastqsplitter -i {input} {params.out} -t 1 -c 1
        """


# ---- Demultiplex each splitted R1 and R2 file to generate sample-specific fastq.gz files. ----
rule demultiplex_fastq_split:
    input:
        R1="{sequencing_name}/raw_reads/R1_{scatteritem}.fastq.gz",
        R2="{sequencing_name}/raw_reads/R2_{scatteritem}.fastq.gz",
    output:
        temp(directory("{sequencing_name}/demux_reads_scatter/{scatteritem}/")),
    log:
        "logs/step2_demultiplexing_reads/demultiplex_fastq_split_{sequencing_name}_{scatteritem}.log",
    threads: 1
    resources:
        mem_mb=1024 * 2,
    params:
        path_samples=config["path_samples"],
        path_barcodes=config["path_barcodes"],
    message:
        "Demultiplexing the scattered .fastq.gz files."
    shell:
        "python3.10 {workflow.basedir}/scripts/demux_rocket.py --sequencing_name {wildcards.sequencing_name} --samples {params.path_samples} --barcodes {params.path_barcodes} --r1 {input[0]} --r2 {input[1]} --out {output} &> {log}"


rule gather_demultiplexed_sequencing:
    input:
        gather.fastq_split("{{sequencing_name}}/demux_reads_scatter/{scatteritem}/"),
    output:
        R1_discarded="{sequencing_name}/demux_reads/{sequencing_name}_R1_discarded.fastq.gz",
        R2_discarded="{sequencing_name}/demux_reads/{sequencing_name}_R2_discarded.fastq.gz",
        discarded_log="{sequencing_name}/demux_reads/log_{sequencing_name}_discarded_reads.tsv.gz",
        qc="{sequencing_name}/demux_reads/{sequencing_name}_qc.pickle",
        whitelist_p7="{sequencing_name}/demux_reads/{sequencing_name}_whitelist_p7.txt",
        whitelist_p5="{sequencing_name}/demux_reads/{sequencing_name}_whitelist_p5.txt",
        whitelist_ligation="{sequencing_name}/demux_reads/{sequencing_name}_whitelist_ligation.txt",
        whitelist_rt="{sequencing_name}/demux_reads/{sequencing_name}_whitelist_rt.txt",
    threads: 1
    resources:
        mem_mb=1024 * 2,
    params:
        path_demux_scatter=lambda w: "{sequencing_name}/demux_reads_scatter/".format(
            sequencing_name=w.sequencing_name
        ),
    message:
        "Combining the scattered discarded reads and log and combined the qc pickles."
    shell:
        """
        # Combine pickles.
        python3.10 {workflow.basedir}/scripts/demux_combine_pickles.py --path_demux_scatter {params.path_demux_scatter} --path_out {output.qc}

        # Combine the sequencing-specific R1/R2 discarded reads and logs.
        find ./{wildcards.sequencing_name}/demux_reads_scatter/ -maxdepth 2 -type f -name {wildcards.sequencing_name}_R1_discarded.fastq.gz -print0 | xargs -0 cat > {wildcards.sequencing_name}/demux_reads/{wildcards.sequencing_name}_R1_discarded.fastq.gz
        find ./{wildcards.sequencing_name}/demux_reads_scatter/ -maxdepth 2 -type f -name {wildcards.sequencing_name}_R2_discarded.fastq.gz -print0 | xargs -0 cat > {wildcards.sequencing_name}/demux_reads/{wildcards.sequencing_name}_R2_discarded.fastq.gz
        find ./{wildcards.sequencing_name}/demux_reads_scatter/ -maxdepth 2 -type f -name log_{wildcards.sequencing_name}_discarded_reads.tsv.gz -print0 | xargs -0 cat > {wildcards.sequencing_name}/demux_reads/log_{wildcards.sequencing_name}_discarded_reads.tsv.gz

        # Move the whitelist files.
        find ./{wildcards.sequencing_name}/demux_reads_scatter/ -maxdepth 2 -type f -name {wildcards.sequencing_name}_whitelist_p7.txt -print0 | xargs -0 cat | sort -u > {output.whitelist_p7}
        find ./{wildcards.sequencing_name}/demux_reads_scatter/ -maxdepth 2 -type f -name {wildcards.sequencing_name}_whitelist_p5.txt -print0 | xargs -0 cat | sort -u > {output.whitelist_p5}
        find ./{wildcards.sequencing_name}/demux_reads_scatter/ -maxdepth 2 -type f -name {wildcards.sequencing_name}_whitelist_ligation.txt -print0 | xargs -0 cat | sort -u > {output.whitelist_ligation}
        find ./{wildcards.sequencing_name}/demux_reads_scatter/ -maxdepth 2 -type f -name {wildcards.sequencing_name}_whitelist_rt.txt -print0 | xargs -0 cat | sort -u > {output.whitelist_rt}
        """

rule gather_demultiplexed_samples:
    input:
        gather.fastq_split("{{sequencing_name}}/demux_reads_scatter/{scatteritem}/"),
    output:
        R1="{sequencing_name}/demux_reads/{sample_name}_R1.fastq.gz",
        R2="{sequencing_name}/demux_reads/{sample_name}_R2.fastq.gz",
    threads: 1
    resources:
        mem_mb=1024 * 2,
    message:
        "Combining the sample-specific fastq.fz files."
    shell:
        """
        # Combine files.
        find ./{wildcards.sequencing_name}/demux_reads_scatter/ -maxdepth 2 -type f -name {wildcards.sample_name}_R1.fastq.gz -print0 | xargs -0 cat > {output.R1}
        find ./{wildcards.sequencing_name}/demux_reads_scatter/ -maxdepth 2 -type f -name {wildcards.sample_name}_R2.fastq.gz -print0 | xargs -0 cat > {output.R2}
        """
