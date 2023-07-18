#############################################
#  fastp: Trimming adapters and low-quality reads
#############################################


rule trim_fastp:
    input:
        R1="{sequencing_name}/demux_reads/{sample_name}_R1.fastq.gz",
        R2="{sequencing_name}/demux_reads/{sample_name}_R2.fastq.gz",
    output:
        R1=temp("{sequencing_name}/fastp/{sample_name}_R1.fastq.gz"),
        R2=temp("{sequencing_name}/fastp/{sample_name}_R2.fastq.gz"),
        html="{sequencing_name}/fastp/{sample_name}.html",
        json="{sequencing_name}/fastp/{sample_name}.json",
    log:
        "logs/step3_alignment/fastp_{sequencing_name}_{sample_name}.log",
    params:
        extra=config["settings"]["fastp"],
    threads: 10
    resources:
        mem_mb=1024 * 4,
    message:
        "Trimming adapters and low-quality reads with fastp."
    shell:
        "fastp {params.extra} --detect_adapter_for_pe --qualified_quality_phred 15 --dont_eval_duplication --length_required 10 --html {output.html} --json {output.json} --thread {threads} --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} >& {log}"


#############################################
#  STAR: Alignment
#############################################


rule generate_index_STAR:
    output:
        directory("{sequencing_name}/resources/index_star/{species}/"),
    log:
        "logs/step3_alignment/generate_index_STAR_{sequencing_name}_{species}.log",
    threads: 20
    resources:
        mem_mb=1024 * 50,
    message:
        "Generating (or symlinking) STAR indexes."
    params:
        fasta=lambda w: config["species"][w.species]["genome"],
        gtf=lambda w: config["species"][w.species]["genome_gtf"],
        star_index=lambda w: config["species"][w.species]["star_index"],
        length_R2=config["length_R2"],
        extra=config["settings"]["star_index"],
    run:
        if params.star_index:
            shell("ln -s {input.star_index} {output}")
        else:
            shell(
                "STAR {params.extra} --runThreadN {threads} --runMode genomeGenerate --genomeFastaFiles {params.fasta} --genomeDir {output} --sjdbGTFfile {params.gtf} --sjdbOverhang {params.length_R2} >& {log}"
            )


rule starSolo_align:
    input:
        R1="{sequencing_name}/fastp/{sample_name}_R1.fastq.gz",
        R2="{sequencing_name}/fastp/{sample_name}_R2.fastq.gz",
        index="{sequencing_name}/resources/index_star/{species}/",
        whitelist_p7="{sequencing_name}/demux_reads/{sequencing_name}_whitelist_p7.txt",
        whitelist_p5="{sequencing_name}/demux_reads/{sequencing_name}_whitelist_p5.txt",
        whitelist_ligation="{sequencing_name}/demux_reads/{sequencing_name}_whitelist_ligation.txt",
        whitelist_rt="{sequencing_name}/demux_reads/{sequencing_name}_whitelist_rt.txt",
    output:
        BAM=temp(
            "{sequencing_name}/alignment/{sample_name}_{species}_Aligned.sortedByCoord.out.bam"
        ),
        SJ="{sequencing_name}/alignment/{sample_name}_{species}_SJ.out.tab",
        log1="{sequencing_name}/alignment/{sample_name}_{species}_Log.final.out",
        log2=temp("{sequencing_name}/alignment/{sample_name}_{species}_Log.out"),
        log3=temp(
            "{sequencing_name}/alignment/{sample_name}_{species}_Log.progress.out"
        ),
        dir_tmp=temp(
            directory("{sequencing_name}/alignment/{sample_name}_{species}__STARtmp/")
        ),
        dir_solo=directory(
            "{sequencing_name}/alignment/{sample_name}_{species}_Solo.out/"
        ),
        barcodes_raw="{sequencing_name}/alignment/{sample_name}_{species}_Solo.out/GeneFull/raw/barcodes.tsv",
        barcodes_raw_converted="{sequencing_name}/alignment/{sample_name}_{species}_Solo.out/GeneFull/raw/barcodes_converted.tsv",
        barcodes_filtered="{sequencing_name}/alignment/{sample_name}_{species}_Solo.out/GeneFull/filtered/barcodes.tsv",
        barcodes_filtered_converted="{sequencing_name}/alignment/{sample_name}_{species}_Solo.out/GeneFull/filtered/barcodes_converted.tsv",
    log:
        "logs/step3_alignment/star_align_{sequencing_name}_{sample_name}_{species}.log",
    params:
        sampleName="{sequencing_name}/alignment/{sample_name}_{species}_",
        extra=config["settings"]["star"],
        path_barcodes=config["path_barcodes"],
    threads: 30
    resources:
        mem_mb=1024 * 60,
    message:
        "Aligning reads with STAR."
    shell:
        """
        STAR {params.extra} --genomeDir {input.index} --runThreadN {threads} \
        --readFilesIn {input.R2} {input.R1} --readFilesCommand zcat \
        --soloType CB_UMI_Complex --soloCBmatchWLtype Exact \
        --soloCBposition 0_0_0_9 0_10_0_19 0_20_0_29 0_30_0_39 --soloUMIposition 0_40_0_47 \
        --soloCBwhitelist {input.whitelist_p7} {input.whitelist_p5} {input.whitelist_ligation} {input.whitelist_rt} \
        --soloCellFilter CellRanger2.2 --soloFeatures GeneFull --soloMultiMappers Uniform --soloCellReadStats Standard \
        --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outFileNamePrefix {params.sampleName} \
        --outSAMmultNmax 1 --outSAMstrandField intronMotif --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 \
        --outSAMattributes NH HI AS nM NM MD jM jI MC ch XS CR UR GX GN sM CB UB \
        >& {log}

        # Convert the barcodes to the barcode naming scheme.
        python3.10 {workflow.basedir}/scripts/STARSolo_convertBarcodes.py --starsolo_barcodes {output.barcodes_raw} --barcodes {params.path_barcodes} --out {output.barcodes_raw_converted}
        python3.10 {workflow.basedir}/scripts/STARSolo_convertBarcodes.py --starsolo_barcodes {output.barcodes_filtered} --barcodes {params.path_barcodes} --out {output.barcodes_filtered_converted}
        """

#############################################
#  Sambamba: Marking duplicates and indexing
#############################################


rule sambamba_markdup:
    input:
        "{sequencing_name}/alignment/{sample_name}_{species}_Aligned.sortedByCoord.out.bam",
    output:
        "{sequencing_name}/alignment/{sample_name}_{species}_Aligned_sortedByCoord_markDup.bam",
    log:
        "logs/step3_alignment/sambamba_markdup_{sequencing_name}_{sample_name}_{species}.log",
    threads: 8
    resources:
        mem_mb=1024 * 10,
    message:
        "Marking duplicates."
    shell:
        "sambamba markdup -t {threads} {input} {output} >& {log}"


rule sambamba_index:
    input:
        "{sequencing_name}/alignment/{sample_name}_{species}_Aligned_sortedByCoord_markDup.bam",
    output:
        "{sequencing_name}/alignment/{sample_name}_{species}_Aligned_sortedByCoord_markDup.bam.bai",
    log:
        "logs/step3_alignment/sambamba_index_{sequencing_name}_{sample_name}_{species}.log",
    threads: 8
    resources:
        mem_mb=1024 * 2,
    message:
        "Indexing BAM."
    shell:
        "sambamba index -t {threads} {input} {output} >& {log}"


#############################################
# Generate the sci-dashboard report.
#############################################

# Get the samples for a given sequencing run (and sci-dash).
def getsamples_sequencing(wildcards):
    x = samples_unique[samples_unique["sequencing_name"] == wildcards.sequencing_name]
    
    files = ["{sequencing_name}/alignment/{sample_name}_{species}_Aligned_sortedByCoord_markDup.bam.bai".format(
        sequencing_name=sequencing_name,
        sample_name=sample_name,
        species=species,
    )
    for sequencing_name, sample_name, species in zip(
        x["sequencing_name"],
        x["sample_name"],
        x["species"],
    )]

    return files

rule sci_dash:
    input:
        lambda w: getsamples_sequencing(w),
        qc="{sequencing_name}/demux_reads/{sequencing_name}_qc.pickle",
    output:
        dash_folder=directory("{sequencing_name}/sci-dash/"),
        dash_json="{sequencing_name}/sci-dash/js/qc_data.js",
    threads: 1
    resources:
        mem_mb=1024 * 2,
    message:
        "Generating sci-dashboard report."
    shell:
        """
        # Generate the sci-dashboard report.
        cp -R {workflow.basedir}/scirocket-dash/* {output.dash_folder}

        # Combine the sample-specific QC and STARSolo metrics.
        python3.10 {workflow.basedir}/scripts/demux_dash.py --path_out {output.dash_json} --path_pickle {input.qc} --path_star {wildcards.sequencing_name}/alignment/
        """
