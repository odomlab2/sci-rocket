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
        mem_mb=1024 * 10,
    message:
        "Trimming adapters and low-quality reads with fastp."
    shell:
        "fastp {params.extra} --detect_adapter_for_pe --qualified_quality_phred 15 --dont_eval_duplication --length_required 10 --html {output.html} --json {output.json} --thread {threads} --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} >& {log}"


#############################################
#  STAR: Alignment
#############################################


rule generate_index_STAR:
    output:
        directory("resources/index_star/{species}/"),
    log:
        "logs/step3_alignment/generate_index_STAR_{species}.log",
    threads: 20
    resources:
        mem_mb=1024 * 50,
    message:
        "Generating (or symlinking) STAR indexes."
    params:
        fasta=lambda w: config["species"][w.species]["genome"],
        gtf=lambda w: config["species"][w.species]["genome_gtf"],
        star_index=lambda w: config["species"][w.species]["star_index"],
    run:
        if params.star_index:
            shell("ln -s {input.star_index} {output}")
        else:
            shell(
                "STAR --runThreadN {threads} --runMode genomeGenerate --genomeFastaFiles {params.fasta} --genomeDir {output} --sjdbGTFfile {params.gtf} --sjdbOverhang 99 >& {log}"
            )


rule starSolo_align:
    input:
        R1="{sequencing_name}/fastp/{sample_name}_R1.fastq.gz",
        R2="{sequencing_name}/fastp/{sample_name}_R2.fastq.gz",
        index="resources/index_star/{species}/",
        whitelist="{sequencing_name}/demux_reads/{sample_name}_whitelist.txt",
    output:
        BAM=temp("{sequencing_name}/alignment/{sample_name}_{species}_Aligned.sortedByCoord.out.bam"),
        SJ="{sequencing_name}/alignment/{sample_name}_{species}_SJ.out.tab",
        log1="{sequencing_name}/alignment/{sample_name}_{species}_Log.final.out",
        log2=temp("{sequencing_name}/alignment/{sample_name}_{species}_Log.out"),
        log3=temp("{sequencing_name}/alignment/{sample_name}_{species}_Log.progress.out"),
        dir_tmp=temp(directory("{sequencing_name}/alignment/{sample_name}_{species}__STARtmp/")),
        dir_solo=directory("{sequencing_name}/alignment/{sample_name}_{species}_Solo.out/")
    log:
        "logs/step3_alignment/star_align_{sequencing_name}_{sample_name}_{species}.log",
    params:
        sampleName="{sequencing_name}/alignment/{sample_name}_{species}_",
        extra=config["settings"]["star"],
    threads: 30
    resources:
        mem_mb=1024 * 40,
    message:
        "Aligning reads with STAR."
    shell:
        """
        STAR {params.extra} --genomeDir {input.index} --runThreadN {threads} \
        --readFilesIn {input.R2} {input.R1} --readFilesCommand zcat \
        --soloType CB_UMI_Complex --soloCBmatchWLtype Exact --soloCBwhitelist {input.whitelist} --soloCBposition  0_0_0_39 --soloUMIposition 0_40_0_47 \
        --soloCellFilter CellRanger2.2 --soloFeatures GeneFull --soloMultiMappers Uniform --soloCellReadStats Standard \
        --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outFileNamePrefix {params.sampleName} \
        --outSAMmultNmax 1 --outSAMstrandField intronMotif \
        --outSAMattributes NH HI AS nM NM MD jM jI MC ch XS CR UR GX GN sM CB UB \
        >& {log}
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

rule sci_dash:
    input:
        qc="{sequencing_name}/demux_reads/{sequencing_name}_qc.pickle",
        log_star="{sequencing_name}/alignment/{sample_name}_{species}_Log.final.out"
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
        # python3.10 {workflow.basedir}/scripts/demux_dash.py --path_out {output.dash_json} --path_qc {input.qc} --path_star {input.log_star}
        """
