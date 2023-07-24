# Rules related to converting bcl files to fastq files with p5 and p7 indexes within the read name.
def get_bcl2fastq_input(sequencing_name):
    """Return the path to the bcl files for a given sequencing run."""
    return samples_unique.query("sequencing_name == @sequencing_name").path_bcl.values[
        0
    ]


rule make_fake_samplesheet:
    output:
        sample_sheet=temp("{sequencing_name}/raw_reads/fake.csv"),
    params:
        path_out="{sequencing_name}/raw_reads/",
    shell:
        """
        # Generate fake sample-sheet to allow indexes to be added to R1/R2.
        mkdir -p {params.path_out}
        echo -e "[DATA]\nLane,Sample_ID,Sample_Name,index,index2\n,fake,fake,NNNNNNNNNN,NNNNNNNNNN" > {output.sample_sheet}
        """


rule bcl2fastq:
    input:
        path_bcl=lambda w: get_bcl2fastq_input(w.sequencing_name),
        fake_sample_sheet="{sequencing_name}/raw_reads/fake.csv",
    output:
        R1=temp("{sequencing_name}/raw_reads/Undetermined_S0_R1_001.fastq.gz"),
        R2=temp("{sequencing_name}/raw_reads/Undetermined_S0_R2_001.fastq.gz"),
        dReports=temp(directory("{sequencing_name}/raw_reads/Reports")),
        dStats=temp(directory("{sequencing_name}/raw_reads/Stats")),
    log:
        "logs/step1_bcl2fastq/bcl2fastq_{sequencing_name}.log",
    threads: 25
    resources:
        mem_mb=1024 * 40,
    params:
        path_out="{sequencing_name}/raw_reads/",
        extra=config["settings"]["bcl2fastq"],
    shell:
        """
        bcl2fastq \
        {params.extra} \
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
        --no-lane-splitting \
        --minimum-trimmed-read-length 15 \
        --mask-short-adapter-reads 15 &> {log}
        """
