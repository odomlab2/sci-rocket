#   * Convert BCL files to Undetermined.fastq.gz files with p5+p7 in the read-name. *
#
#   1. make_fake_samplesheet:       Generate fake sample-sheet to allow indexes to be added to R1/R2.
#   2. bcl2fastq:                   Convert bcl to fastq with p5 and p7 indexes within the read name.
#############

def get_bcl2fastq_input(sequencing_name, experiment_name):
    """RReturn the path to the bcl file for a given sequencing run."""
    return samples_unique.query("sequencing_name == @sequencing_name & experiment_name == @experiment_name").path_bcl.values[0]
    

rule make_fake_samplesheet:
    output:
        sample_sheet=temp("{experiment_name}/raw_reads/{sequencing_name}/fake.csv"),
    params:
        path_out="{experiment_name}/raw_reads/{sequencing_name}/",
    message: "Generate fake sample-sheet to allow indexes to be added to R1/R2."
    shell:
        """
        # Generate fake sample-sheet to allow indexes to be added to R1/R2.
        mkdir -p {params.path_out}
        echo -e "[DATA]\nLane,Sample_ID,Sample_Name,index,index2\n,fake,fake,NNNNNNNNNN,NNNNNNNNNN" > {output.sample_sheet}
        """


rule bcl2fastq:
    input:
        path_bcl=lambda w: get_bcl2fastq_input(w.sequencing_name, w.experiment_name),
        fake_sample_sheet="{experiment_name}/raw_reads/{sequencing_name}/fake.csv",
    output:
        R1=temp("{experiment_name}/raw_reads/{sequencing_name}/Undetermined_S0_R1_001.fastq.gz"),
        R2=temp("{experiment_name}/raw_reads/{sequencing_name}/Undetermined_S0_R2_001.fastq.gz"),
        dReports=temp(directory("{experiment_name}/raw_reads/{sequencing_name}/Reports")),
        dStats=temp(directory("{experiment_name}/raw_reads/{sequencing_name}/Stats")),
    log:
        "logs/step1_bcl2fastq/bcl2fastq_{experiment_name}_{sequencing_name}.log",
    threads: 40
    resources:
        mem_mb=1024 * 40,
    params:
        path_out="{experiment_name}/raw_reads/{sequencing_name}/",
        extra=config["settings"]["bcl2fastq"],
    conda:
        "envs/sci-rocket.yaml",
    message: "Converting bcl to fastq with p5 and p7 indexes within the read name ({wildcards.experiment_name}: {wildcards.sequencing_name})."
    shell:
        """
        bcl2fastq \
        {params.extra} \
        -R {input.path_bcl} \
        --sample-sheet {input.fake_sample_sheet} \
        --output-dir {params.path_out} \
        --loading-threads 8 \
        --processing-threads 30 \
        --writing-threads 2 &> {log}
        """

def get_sequencing_runs(experiment_name):
    """Return a list of sequencing runs."""
    return samples_unique.query("experiment_name == @experiment_name").sequencing_name.unique().tolist()

rule merge_sequencing_runs:
    input:
        R1=lambda w: expand("{experiment_name}/raw_reads/{sequencing_name}/Undetermined_S0_R1_001.fastq.gz", experiment_name=w.experiment_name, sequencing_name=get_sequencing_runs(w.experiment_name)),
        R2=lambda w: expand("{experiment_name}/raw_reads/{sequencing_name}/Undetermined_S0_R2_001.fastq.gz", experiment_name=w.experiment_name, sequencing_name=get_sequencing_runs(w.experiment_name)),
    output:
        R1="{experiment_name}/raw_reads/Undetermined_S0_R1_001.fastq.gz",
        R2="{experiment_name}/raw_reads/Undetermined_S0_R2_001.fastq.gz",
    threads: 1
    resources:
        mem_mb=1024 * 2,
    params:
        total_sequencing_runs=lambda w: len(get_sequencing_runs(w.experiment_name)),
    message: "Merge sequencing runs ({wildcards.experiment_name})."
    shell:
        """
        # If only one sequencing run, then just move it.
        if [ {params.total_sequencing_runs} -eq 1 ]; then
            mv {input.R1} {output.R1}
            mv {input.R2} {output.R2}
        else
            cat {input.R1} > {output.R1}
            cat {input.R2} > {output.R2}
        fi
        """