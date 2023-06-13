#############
#   Rules related to demultiplexing the initial sequencing reads.
#############

# Demultiplexing an entire sequencing run to generate sample-specific fastq files.
rule demultiplex_fastq:
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
        directory("demultiplex_fastq/untrimmed/{sequencing_name}/")
    shell:
        "1"

rule demultiplex_samples:
    input:
        directory("demultiplex_fastq/untrimmed/{sequencing_name}/")
    output:
        "demultiplex_fastq/untrimmed/{sequencing_name}_{sample_name}_R1.fq.gz"
    shell:
        "1"