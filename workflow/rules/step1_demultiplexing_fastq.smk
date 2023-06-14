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
        directory("demultiplex_fastq/raw/{sequencing_name}/")
    log: 
        "logs/demultiplex_fastq/{sequencing_name}.log"
    params:
        path_samples=config['path_samples'],
        path_barcodes=config['path_barcodes']
    shell:
        "python3.10 {workflow.basedir}/scripts/sciseq_sample_demultiplexing -sequencing_name {wildcards.sequencing_name} -path_samples {params.path_samples} -path_barcodes {params.path_barcodes} -path_r1 {input[0]} -path_r2 {input[1]} -path_out {output} &> {log}"

rule demultiplex_samples:
    input:
        "demultiplex_fastq/raw/{sequencing_name}/"
    output:
        R1 = "demultiplex_fastq/fastp/{sequencing_name}_{sample_name}_R1.fq.gz",
        R2 = "demultiplex_fastq/fastp/{sequencing_name}_{sample_name}_R2.fq.gz"
    shell:
        "{output.R1} + {output.R2}"