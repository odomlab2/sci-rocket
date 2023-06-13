#############
#   Rules related to demultiplexing the initial sequencing reads.
#############


rule demultiplex_samples:
    input:
        lambda w: [
            "{fastq}/{sequencing_name}-LR-67093_R1.fastq.gz".format(
                fastq=config["dir_fastq"], sequencing_name=w.sequencing_name
            ),
            "{fastq}/{sequencing_name}-LR-67093_R2.fastq.gz".format(
                fastq=config["dir_fastq"], sequencing_name=w.sequencing_name
            ),
        ],
    output:
        [
            "demultiplex_fastq/untrimmed/{sequencing_name}_{sample_name}_R1.fq.gz",
            "demultiplex_fastq/untrimmed/{sequencing_name}_{sample_name}_R2.fq.gz",
        ],
    shell:
        "1"
