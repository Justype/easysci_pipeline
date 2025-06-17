from os import path

configfile: "config.yaml"

rule star:
    params:
        output_prefix = lambda wildcards: path.join(config["output_dir"], wildcards.species, "bam", wildcards.prefix),
        count_log = lambda wildcards: path.join(config["log_dir"], wildcards.species, f"{wildcards.prefix}.log"),
    input:
        genome_dir = path.join(config["star_index_folder"], "{species}"),
        r1 = path.join(config["output_dir"], "{species}/fastqs_trimmed/{prefix}_R1_val_1.fq.gz"),
        r2 = path.join(config["output_dir"], "{species}/fastqs_trimmed/{prefix}_R2_val_2.fq.gz"),
    output:
        path.join(config["output_dir"], "{species}/bam/{prefix}.Aligned.out.bam"),
    log:
        path.join(config["log_dir_run"], "{species}/{prefix}/3_star.log")
    conda:
        "../envs/sci_rna.yaml"
    threads:
        config["threads"]["star"]
    resources:
        mem_mb = lambda wildcards, threads: 40960 + threads * 150,  # 40GB + 150MB per thread
        runtime = 720,  # 12 hours
    shell:
        """
        STAR \\
            --runThreadN {threads} \\
            --genomeDir {input.genome_dir} \\
            --readFilesCommand zcat \\
            --readFilesIn {input.r1} {input.r2} \\
            --outSAMtype BAM Unsorted \\
            --outSAMstrandField intronMotif \\
            --outFileNamePrefix {params.output_prefix}. &> {log}
        
        # remove _STARtmp directory if it exists
        if [ -d "{params.output_prefix}._STARtmp" ]; then
            rm -r {params.output_prefix}._STARtmp
        fi
        """

        # Record stats
        """
        # Remove previous count log if it exists
        sed -i '/reads_mapped/d' {params.count_log}
        sed -i '/alignments/d' {params.count_log}
        
        # Count number of mapped reads (64: first read, 260=256+4: primary map + not unmapped)
        n_mapped=$(samtools view -@ {threads} -c -f 64 -F 260 {output})
        n_alignments=$(samtools view -@ {threads} -c -f 64 -F 4 {output})
        printf "$n_mapped\\treads_mapped\\n" >> {params.count_log}
        printf "$n_alignments\\talignments\\n" >> {params.count_log}
        """
