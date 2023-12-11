#   * Generate QC-dashboard (sci-dash). *
#
#   1. sci-dash:                Generates the QC-dashboard for a given sequencing run.
#############

# Get the samples for a given sequencing run (and sci-dash).
def getsamples_sequencing(wildcards):
    x = samples_unique[samples_unique["experiment_name"] == wildcards.experiment_name]
    
    files = ["{experiment_name}/alignment/{sample_name}_{species}_Aligned.sortedByCoord.out.bam.bai".format(
        experiment_name=experiment_name,
        sample_name=sample_name,
        species=species,
    )
    for experiment_name, sample_name, species in zip(
        x["experiment_name"],
        x["sample_name"],
        x["species"],
    )]

    return files

rule sci_dash:
    input:
        lambda w: getsamples_sequencing(w),
        qc="{experiment_name}/demux_reads/{experiment_name}_qc.pickle"
    output:
        dash_folder=directory("{experiment_name}/sci-dash/"),
        dash_json="{experiment_name}/sci-dash/js/qc_data.js",
        metrics_hashing="{experiment_name}/hashing/{experiment_name}_hashing_metrics.tsv"
    threads: 1
    resources:
        mem_mb=1024 * 2,
    conda:
        "envs/sci-rocket.yaml",
    message:
        "Generating sci-dashboard report ({wildcards.experiment_name})."
    shell:
        """
        # Generate the sci-dashboard report.
        cp -R {workflow.basedir}/scirocket-dash/* {output.dash_folder}

        # Combine the sample-specific QC and STARSolo metrics.
        python3.10 {workflow.basedir}/rules/scripts/demultiplexing/demux_dash.py \
        --path_out {output.dash_json} \
        --path_pickle {input.qc} \
        --path_star {wildcards.experiment_name}/alignment/ \
        --path_hashing {output.metrics_hashing}

        # Remove all empty (leftover) folders.
        find . -empty -type d -delete
        """
