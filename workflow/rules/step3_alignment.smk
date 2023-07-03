#############################################
#  fastp: Trimming adapters and low-quality reads
#############################################


rule trim_fastp:
    input:
        R2="{sequencing_name}/demux_reads/{sample_name}_R2.fastq.gz",
    output:
        R2=temp("{sequencing_name}/fastp/{sample_name}_R2.fastq.gz"),
        html="{sequencing_name}/fastp/{sample_name}.html",
        json="{sequencing_name}/fastp/{sample_name}.json",
    log:
        "logs/step3_alignment/fastp_{sequencing_name}_{sample_name}.log",
    params:
        extra=config["settings"]["fastp"],
    threads: 10
    message:
        "Trimming adapters and low-quality reads with fastp."
    shell:
        "fastp {params.extra} --html {output.html} --json {output.json} --thread {threads} --in1 {input.R2} --out1 {output.R2} >& {log}"


#############################################
#  STAR: Alignment
#############################################


rule generate_index_STAR:
    output:
        directory("resources/index_star/{species}/"),
    log:
        "logs/step3_alignment/generate_index_STAR_{species}.log",
    threads: 8
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


rule star_align:
    input:
        R2="{sequencing_name}/fastp/{sample_name}_R2.fastq.gz",
        index="resources/index_star/{species}/",
    output:
        BAM=temp(
            "{sequencing_name}/alignment/{sample_name}_{species}_Aligned.sortedByCoord.out.bam"
        ),
        log1=temp("{sequencing_name}/alignment/{sample_name}_{species}_Log.final.out"),
        log2=temp("{sequencing_name}/alignment/{sample_name}_{species}_Log.out"),
        log3=temp("{sequencing_name}/alignment/{sample_name}_{species}_Log.progress.out"),
        SJ="{sequencing_name}/alignment/{sample_name}_{species}_SJ.out.tab",
        dir1=temp(directory("{sequencing_name}/alignment/{sample_name}_{species}__STARgenome/")),
        dir2=temp(directory("{sequencing_name}/alignment/{sample_name}_{species}__STARpass1/")),
        dir3=temp(directory("{sequencing_name}/alignment/{sample_name}_{species}__STARtmp/")),
    log:
        "logs/step3_alignment/star_align_{sequencing_name}_{sample_name}_{species}.log",
    params:
        rg="ID:{sequencing_name}_{sample_name} SM:{sample_name} PL:ILLUMINA",
        sampleName="{sequencing_name}/alignment/{sample_name}_",
        extra=config["settings"]["star"],
    threads: 8
    message:
        "Aligning reads with STAR."
    shell:
        """
        STAR {params.extra} --genomeDir {input.index} --runThreadN {threads} \
        --readFilesIn {input.R2} --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outFileNamePrefix {params.sampleName} \
        --outSAMattrRGline {params.rg} --outSAMattributes NH HI AS nM NM MD jM jI MC ch XS \
        --twopassMode Basic --twopass1readsN -1 \
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
    message:
        "Indexing BAM."
    shell:
        "sambamba index -t {threads} {input} {output} >& {log}"
