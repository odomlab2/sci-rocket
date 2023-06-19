# Rules related to converting bcl files to fastq files with p5 and p7 indexes within the read name.
def get_bcl2fastq_input(sequencing_name):
    """Return the path to the bcl files for a given sequencing run."""
    return samples_unique.query("sequencing_name == @sequencing_name").path_bcl.values[
        0
    ]


rule make_fake_samplesheet:
    output:
        sample_sheet=temp("fastq/{sequencing_name}/fake.csv"),
    params:
        path_out="fastq/{sequencing_name}/raw/",
    shell:
        """
        # Generate fake sample-sheet to allow indexes to be added to R1/R2.
        mkdir -p {params.path_out}
        echo -e "[DATA]\nLane,Sample_ID,Sample_Name,index,index2\n,fake,fake,NNNNNNNNNN,NNNNNNNNNN" > {output.sample_sheet}
        """


rule bcl2fastq:
    input:
        path_bcl=lambda w: get_bcl2fastq_input(w.sequencing_name),
        fake_sample_sheet="fastq/{sequencing_name}/fake.csv",
    output:
        R1="fastq/{sequencing_name}/raw/Undetermined_S0_L001_R1_001.fastq.gz",
        R2="fastq/{sequencing_name}/raw/Undetermined_S0_L001_R2_001.fastq.gz",
    log:
        "logs/demux/bcl2fastq_{sequencing_name}.log",
    resources:
        mem_mb=10000,
    params:
        path_out="fastq/{sequencing_name}/raw/",
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
