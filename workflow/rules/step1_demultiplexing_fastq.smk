#############
#   Rules related to demultiplexing the initial sequencing reads.
#############

def get_bcl2fastq_input(sequencing_name):
    """Return the path to the bcl files for a given sequencing run."""
    return samples_unique.query("sequencing_name == @sequencing_name").path_bcl.values[0]

rule make_fake_samplesheet:
    output:
        sample_sheet = temp("fastq/{sequencing_name}/fake.csv"),
    params:
        path_out = "fastq/{sequencing_name}/raw/"
    shell:
        """
        # Generate fake sample-sheet to allow indexes to be added to R1/R2.
        mkdir -p {params.path_out}
        echo -e "[DATA]\nLane,Sample_ID,Sample_Name,index,index2\n,fake,fake,NNNNNNNNNN,NNNNNNNNNN" > {output.sample_sheet}
        """

rule bcl2fastq:
    input:
        path_bcl = lambda w: get_bcl2fastq_input(w.sequencing_name),
        fake_sample_sheet = "fastq/{sequencing_name}/fake.csv",
    output:
        R1 = "fastq/{sequencing_name}/raw/Undetermined_S0_L001_R1_001.fastq.gz",
        R2 = "fastq/{sequencing_name}/raw/Undetermined_S0_L001_R2_001.fastq.gz",
        R1_fake = temp("fastq/{sequencing_name}/raw/fake_S1_L001_R1_001.fastq.gz"),
        R2_fake = temp("fastq/{sequencing_name}/raw/fake_S1_L001_R2_001.fastq.gz"),
    log:
        "logs/demux/bcl2fastq_{sequencing_name}.log"
    resources:
        mem_mb=10000
    params:
        path_out = "fastq/{sequencing_name}/raw/"
    threads: 20
    shell:
        """
        bcl2fastq \
        -R {input.path_bcl} \
        --sample-sheet {input.fake_sample_sheet} \
        --output-dir {params.path_out} \
        --loading-threads {threads} \
        --processing-threads {threads}   \
        --writing-threads {threads}  \
        --barcode-mismatches 1 \
        --ignore-missing-positions \
        --ignore-missing-controls \
        --ignore-missing-filter \
        --ignore-missing-bcls \
        --minimum-trimmed-read-length 15 \
        --mask-short-adapter-reads 15 &> {log}
        """

# Demultiplexing an entire sequencing run to generate sample-specific fastq files.
rule demultiplex_fastq:
    input:
        lambda w: [
            "fastq/{sequencing_name}/raw/Undetermined_S0_L001_R1_001.fastq.gz".format(
                sequencing_name=w.sequencing_name
            ),
            "fastq/{sequencing_name}/raw/Undetermined_S0_L001_R2_001.fastq.gz".format(
                sequencing_name=w.sequencing_name
            ),
        ],
    output:
        directory("fastq/{sequencing_name}/samples/")
    log: 
        "logs/demux/demultiplex_fastq_{sequencing_name}.log"
    resources:
        mem_mb=10000
    params:
        path_samples=config['path_samples'],
        path_barcodes=config['path_barcodes']
    shell:
        "python3.10 {workflow.basedir}/scripts/demultiplexing_samples.py --sequencing_name {wildcards.sequencing_name} --samples {params.path_samples} --barcodes {params.path_barcodes} --r1 {input[0]} --r2 {input[1]} --out {output} &> {log}"

rule demultiplex_samples:
    input:
        "fastq/{sequencing_name}/samples/"
    output:
        R1 = "fastq/{sequencing_name}/samples/fastp/{sequencing_name}_{sample_name}_R1.fq.gz",
        R2 = "fastq/{sequencing_name}/samples/fastp/{sequencing_name}_{sample_name}_R2.fq.gz"
    shell:
        "touch {output.R1} {output.R2}"