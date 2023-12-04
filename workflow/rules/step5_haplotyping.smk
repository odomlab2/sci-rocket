#   * Perform additional haplotype-specific analysis (optional). *
#
#   1. mgp_download:                Download the strain-specific mice SNPs from the Mouse Genomes Project (MGP).
#   2. mgp_chr_prefix:              Add chr-prefix to the MGP database.
#   3. generate_hybrid_vcf:         Retrieve informative heterozygous SNPs strain1 x strain2 from the MGP data.
#   4. normalize_hybrid_vcf:        Normalize the F1 VCF: left-align and at least two differing alleles.
#   5. download_repeatmasker:       Download Repeatmasker files (GRCm39)
#   6. filter_repeatmasker:         Retain SNPs and remove SNPs that are in known repetitive regions (simple repeats; e.g. AAAAAA / TTTTTT etc.).
#   7. run_haplotag:                Tag SNP-overlapping reads with haplotype origin.
#   8. haplotype_split_h1:          Split haplotype-tagged reads into haplotype 1.
#   9. haplotype_split_h2:          Split haplotype-tagged reads into haplotype 2.
#   10. haplotype_split_ua:         Split haplotype-tagged reads into unassigned reads.
#   11. count_haplotagged_reads:    Count haplotype-tagged reads.
#   12. join_counts:                Join haplotype-tagged read counts into a single file.
#############

rule mgp_download:
    output:
        mgp_snp=temp("resources/MGP/mgp_REL2021_snps.vcf.gz"),
        mgp_snp_idx=temp("resources/MGP/mgp_REL2021_snps.vcf.gz.csi"),
        mgp_indel=temp("resources/MGP/mgp_REL2021_indels.vcf.gz"),
        mgp_indel_idx=temp("resources/MGP/mgp_REL2021_indels.vcf.gz.csi"),
        mgp_combined=temp("resources/MGP/mgp_REL2021_snps_indels.vcf.gz"),
    log:
        "logs/haplotyping/prepare_mgp.log",
    threads: 10
    params:
        path_mgp=config["path_mgp"],
    conda:
        "envs/sci-haplotyping.yaml",
    message:
        "Downloading the MGP database (~30GB). This will take a moment."
    shell:
        """
        # If the MGP database (path_mgp) is already downloaded, symlink it.
        if [[ ! -z {params.path_mgp} ]]; then
            touch {output.mgp_snp}
            touch {output.mgp_indel}
            touch {output.mgp_snp_idx}
            touch {output.mgp_indel_idx}
            touch {output.mgp_combined}
        else
            wget -O {output.mgp_snp} https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz >& {log}
            wget -O {output.mgp_indel} https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_indels.vcf.gz >& {log}

            wget -O {output.mgp_snp_idx} https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz.csi >& {log}
            wget -O {output.mgp_indel_idx} https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_indels.vcf.gz.csi >& {log}

            # Combine SNP and InDels.
            bcftools concat --threads {threads} -a {output.mgp_snp} {output.mgp_indel} -O z -o {output.mgp_combined}
        fi
        """


rule mgp_chr_prefix:
    input:
        "resources/MGP/mgp_REL2021_snps_indels.vcf.gz",
    output:
        vcf="resources/MGP/mgp_REL2021_snps_indels_chr_prefix.vcf.gz",
        tbi="resources/MGP/mgp_REL2021_snps_indels_chr_prefix.vcf.gz.tbi",
    threads: 8
    conda:
        "envs/sci-haplotyping.yaml",
    params:
        path_mgp=config["path_mgp"],
    message:
        "Adding chr-prefix to the MGP database: {input}"
    shell:
        """
        if [[ ! -z {params.path_mgp} ]]; then
            ln -s {params.path_mgp} {output.vcf}
            ln -s {params.path_mgp}.tbi {output.tbi}
        else
            zcat {input} | awk '{{if($1 ~ "^#") {{gsub("contig=<ID=", "contig=<ID=chr"); gsub("contig=<ID=chrMT", "contig=<ID=chrM"); print}} else {{gsub("^MT", "M"); print "chr"$0}}}}' | bcftools view -O z -o {output.vcf}
            bcftools index -t {output.vcf}
        fi
        """


rule generate_hybrid_vcf:
    input:
        mgp="resources/MGP/mgp_REL2021_snps_indels_chr_prefix.vcf.gz",
        mgp_index="resources/MGP/mgp_REL2021_snps_indels_chr_prefix.vcf.gz.tbi",
    output:
        vcf=temp("resources/MGP/{strain1}_{strain2}_hybrid.vcf.gz"),
    log:
        "logs/haplotyping/generate_hybrid_vcf_{strain1}_{strain2}.log",
    threads: 2
    resources:
        mem_mb=1024 * 2,
    conda:
        "envs/sci-haplotyping.yaml",
    message:
        "Generating cross-hybrid VCF file: {wildcards.strain1} x {wildcards.strain2}"
    shell:
        """
        # If strain1 or strain2 is B6, only use the other strain.
        if [[ "{wildcards.strain1}" == "B6" ]]; then
            python3 {workflow.basedir}/rules/scripts/haplotyping/generate_crosshybrid.py \
                --haplotype {input.mgp} \
                --h1 {wildcards.strain2} \
                --out {output.vcf} \
                --highconfidence \
                >& {log}
        elif [[ "{wildcards.strain2}" == "B6" ]]; then
            python3 {workflow.basedir}/rules/scripts/haplotyping/generate_crosshybrid.py \
                --haplotype {input.mgp} \
                --h1 {wildcards.strain1} \
                --out {output.vcf} \
                --highconfidence \
                >& {log}
        else
            python3 {workflow.basedir}/rules/scripts/haplotyping/generate_crosshybrid.py \
                --haplotype {input.mgp} \
                --h1 {wildcards.strain1} \
                --h2 {wildcards.strain2} \
                --out {output.vcf} \
                --highconfidence \
                >& {log}
        fi
        """


rule normalize_hybrid_vcf:
    input:
        vcf="resources/MGP/{strain1}_{strain2}_hybrid.vcf.gz",
    output:
        vcf=temp("resources/MGP/{strain1}_{strain2}_hybrid_norm.vcf.gz"),
        idx=temp("resources/MGP/{strain1}_{strain2}_hybrid_norm.vcf.gz.tbi"),
    threads: 1
    params:
        fasta=lambda w: config["species"]["mouse"]["genome"],
    conda:
        "envs/sci-haplotyping.yaml",
    message:
        "Normalizing cross-hybrid VCF file: {wildcards.strain1} x {wildcards.strain2}"
    shell:
        """
        bcftools norm  -f {params.fasta} {input.vcf} -m -any | \
        bcftools view -I --trim-alt-alleles | \
        bcftools view --min-alleles 2 | \
        bcftools norm -m+ -o {output.vcf} -O z -f {params.fasta}

        bcftools index -t {output.vcf}
        """


rule download_repeatmasker:
    output:
        repeatmasker=temp("resources/MGP/rmsk.bed"),
    threads: 1
    resources:
        mem_mb=1024 * 2,
    params:
        url_repeatmasker=config["url_repeatmasker"],
    conda:
        "envs/sci-haplotyping.yaml",
    message: "Downloading Repeatmasker file (GRCm39)."
    shell:
        """
        wget -O resources/MGP/rmsk.txt.gz {params.url_repeatmasker}

        # Convert to BED format.
        zgrep -E "\(.\)n" resources/MGP/rmsk.txt.gz | awk '{{print $6"\t"$7"\t"$8"\t"$11"\t"$12"\t"$13}}' > {output.repeatmasker}

        # Remove temporary file.
        rm resources/MGP/rmsk.txt.gz
        """


rule filter_repeatmasker:
    input:
        vcf="resources/MGP/{strain1}_{strain2}_hybrid_norm.vcf.gz",
        idx="resources/MGP/{strain1}_{strain2}_hybrid_norm.vcf.gz.tbi",
        repeatmasker="resources/MGP/rmsk.bed",
    output:
        vcf="resources/MGP/{strain1}_{strain2}_hybrid_norm_SNPs_norepeats.vcf.gz",
        idx="resources/MGP/{strain1}_{strain2}_hybrid_norm_SNPs_norepeats.vcf.gz.tbi",
    threads: 1
    conda:
        "envs/sci-haplotyping.yaml",
    shell:
        """
        bcftools view --types snps -T ^{input.repeatmasker} -o {output.vcf} -O z {input.vcf}
        bcftools index -t {output.vcf}
        """

rule run_haplotag:
    input:
        bam="{sequencing_name}/alignment/{sample_name}_mouse_Aligned.sortedByCoord.out.bam",
        vcf="resources/MGP/{strain1}_{strain2}_hybrid_norm_SNPs_norepeats.vcf.gz",
        idx="resources/MGP/{strain1}_{strain2}_hybrid_norm_SNPs_norepeats.vcf.gz.tbi",
    output:
        bam="{sequencing_name}/alignment/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX.bam",
        bai="{sequencing_name}/alignment/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX.bam.bai",
    log:
        "logs/haplotyping/haplotag_{sequencing_name}_{sample_name}_{strain1}_{strain2}.log",
    threads: 10
    resources:
        mem_mb=1024 * 40,
    params:
        fasta=lambda w: config["species"]["mouse"]["genome"],
    conda:
        "envs/sci-haplotyping.yaml",
    message:
        "Tagging (HP tag) SNP-overlapping reads with haplotype-origin: {input.bam}."
    shell:
        """
        # Haplotype tagging.
        whatshap haplotag --regions chrX --sample F1 --ignore-read-groups -o {output.bam} --reference {params.fasta} --output-threads {threads} {input.vcf} {input.bam} >& {log}

        # Index BAM file.
        sambamba index -t {threads} {output.bam}
        """


rule haplotype_split_h1:
    input:
        bam="{sequencing_name}/alignment/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX.bam",
        bai="{sequencing_name}/alignment/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX.bam.bai",
    output:
        bam=temp("{sequencing_name}/haplotyping/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX_h1.bam"),
        bai=temp("{sequencing_name}/haplotyping/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX_h1.bam.bai"),
    threads: 10
    resources:
        mem_mb=1024 * 10,
    conda:
        "envs/sci-haplotyping.yaml",
    shell:
        """
        sambamba view -t {threads} -h -f bam -F "[HP]==1 and [CR]!=null and [GX]!=null" {input.bam} > {output.bam}
        sambamba index -t {threads} {output.bam}
        """


rule haplotype_split_h2:
    input:
        bam="{sequencing_name}/alignment/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX.bam",
        bai="{sequencing_name}/alignment/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX.bam.bai",
    output:
        bam=temp("{sequencing_name}/haplotyping/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX_h2.bam"),
        bai=temp("{sequencing_name}/haplotyping/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX_h2.bam.bai"),
    threads: 10
    resources:
        mem_mb=1024 * 10,
    conda:
        "envs/sci-haplotyping.yaml",
    shell:
        """
        sambamba view -t {threads} -h -f bam -F "[HP]==2 and [CR]!=null and [GX]!=null" {input.bam} > {output.bam}
        sambamba index -t {threads} {output.bam}
        """


rule haplotype_split_ua:
    input:
        bam="{sequencing_name}/alignment/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX.bam",
        bai="{sequencing_name}/alignment/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX.bam.bai",
    output:
        bam=temp("{sequencing_name}/haplotyping/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX_ua.bam"),
        bai=temp("{sequencing_name}/haplotyping/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX_ua.bam.bai"),
    threads: 10
    resources:
        mem_mb=1024 * 10,
    conda:
        "envs/sci-haplotyping.yaml",
    shell:
        """
        sambamba view -t {threads} -h -f bam -F "[HP]==null and [CR]!=null and [GX]!=null" {input.bam} > {output.bam}
        sambamba index -t {threads} {output.bam}
        """


rule count_haplotagged_reads:
    input:
        bam="{sequencing_name}/haplotyping/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX_{type}.bam",
        bai="{sequencing_name}/haplotyping/{sample_name}_mouse_Aligned.sortedByCoord.out.haplotagged_{strain1}_{strain2}.chrX_{type}.bam.bai"
    output:
        counts="{sequencing_name}/haplotyping/{sample_name}_{strain1}_{strain2}_haplotagged_readcounts_{type}.txt",
    log:
        "logs/haplotyping/count_haplotagged_reads_{sequencing_name}_{sample_name}_{strain1}_{strain2}_{type}.log",
    threads: 2
    resources:
        mem_mb=1024 * 10,
    conda:
        "envs/sci-haplotyping.yaml",
    message:
        "Counting haplotagged reads in {input.bam}."
    shell:
        """
        umi_tools count --per-gene --extract-umi-method=tag --gene-tag=GX --cell-tag=CR --umi-tag=UR --per-cell -I {input.bam} -S {output.counts} >& {log}
        """


rule join_counts:
    input:
        counts_h1="{sequencing_name}/haplotyping/{sample_name}_{strain1}_{strain2}_haplotagged_readcounts_h1.txt",
        counts_h2="{sequencing_name}/haplotyping/{sample_name}_{strain1}_{strain2}_haplotagged_readcounts_h2.txt",
        counts_ua="{sequencing_name}/haplotyping/{sample_name}_{strain1}_{strain2}_haplotagged_readcounts_ua.txt",
    output:
        counts="{sequencing_name}/haplotyping/{sample_name}_{strain1}_{strain2}_haplotagged_readcounts.txt",
    log:
        "logs/haplotyping/retrieve_counts_{sequencing_name}_{sample_name}_{strain1}_{strain2}.log",
    threads: 1
    resources:
        mem_mb=1024 * 20,
    conda:
        "envs/sci-haplotyping.yaml",
    message:
        "Retrieving counts from {input.counts_h1}, {input.counts_h2} and {input.counts_ua}."
    shell:
        """
        python3 {workflow.basedir}/rules/scripts/haplotyping/join_counts.py \
            --h1 {input.counts_h1} \
            --h2 {input.counts_h2} \
            --ua {input.counts_ua} \
            --output {output.counts}
        """
