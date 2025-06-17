from os import path

configfile: "config.yaml"

rule bam_filter_dedup:
    params:
        count_log = lambda wildcards: path.join(config["log_dir"], wildcards.species, f"{wildcards.prefix}.log")
    input:
        path.join(config["output_dir"], "{species}/bam/{prefix}.Aligned.out.bam"),
    output:
        path.join(config["output_dir"], "{species}/bam_filtered/{prefix}.filtered.dedup.bam"),
    log:
        path.join(config["log_dir_run"], "{species}/{prefix}/4_bam_filter_dedup.log")
    conda:
        "../envs/sci_rna.yaml"
    threads:
        config["threads"]["bam_filter_dedup"]
    resources:
        mem_mb = lambda wildcards, threads: threads * 850,  # samtools uses 768M per thread by default
        runtime = 360, # 6 hours
    shell:
        """
        samtools view -@ {threads} -h -q 30 -f 2 -F 780 -b {input} | \\
            samtools sort -@ {threads} -n | \\
            python workflow/scripts/bam_dedup_paired.py --log {params.count_log} \\
            > {output} 2> {log}
        """

