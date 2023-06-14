#############
#   Rules related to demultiplexing the initial sequencing reads.
#############

# Demultiplexing an entire sequencing run to generate sample-specific fastq files.
rule demultiplex_fastq:
    input:
        lambda w: [
            "{fastq}/{sequencing_name}_R1.fastq.gz".format(
                fastq=config["dir_fastq"], sequencing_name=w.sequencing_name
            ),
            "{fastq}/{sequencing_name}_R2.fastq.gz".format(
                fastq=config["dir_fastq"], sequencing_name=w.sequencing_name
            ),
        ],
    output:
        directory("demultiplex_fastq/raw/{sequencing_name}/")
    log: 
        "logs/demultiplex_fastq/{sequencing_name}.log"
    params:
        path_samples=config['path_samples'],
        path_barcodes=config['path_barcodes']
    shell:
        "python3.10 {workflow.basedir}/scripts/demultiplexing_samples.py --sequencing_name {wildcards.sequencing_name} --samples {params.path_samples} --barcodes {params.path_barcodes} --r1 {input[0]} --r2 {input[1]} --out {output} &> {log}"

rule demultiplex_samples:
    input:
        "demultiplex_fastq/raw/{sequencing_name}/"
    output:
        R1 = "demultiplex_fastq/fastp/{sequencing_name}_{sample_name}_R1.fq.gz",
        R2 = "demultiplex_fastq/fastp/{sequencing_name}_{sample_name}_R2.fq.gz"
    shell:
        "echo {output.R1} {output.R2}"