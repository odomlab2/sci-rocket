#############
#   Rules related to demultiplexing the initial sequencing reads.
#############


rule demultiplex_samples:
    input:
        lambda w: [
            "{fastq}/{sample}_R1.fastq.gz".format(
                fastq=config["dir_fastq"], sample=w.sample
            ),
            "{fastq}/{sample}_R2.fastq.gz".format(
                fastq=config["dir_fastq"], sample=w.sample
            ),
        ],
    output:
        [
            "demultiplex_fastq/untrimmed/{sample}_R1.fq.gz",
            "demultiplex_fastq/untrimmed/{sample}_R2.fq.gz",
        ],
    shell:
        "1"
