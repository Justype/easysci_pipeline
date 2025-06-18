from os import path

configfile: "config.yaml"

rule trim:
    params:
        output_folder = lambda wildcards: path.join(config["output_dir"], wildcards.species, "fastqs_trimmed"),
        count_log = lambda wildcards: path.join(config["log_dir"], wildcards.species, f"{wildcards.prefix}.log")
    input:
        r1 = path.join(config["output_dir"], "{species}/fastqs/{prefix}_R1.fastq.gz"),
        r2 = path.join(config["output_dir"], "{species}/fastqs/{prefix}_R2.fastq.gz"),
    output:
        r1 = path.join(config["output_dir"], "{species}/fastqs_trimmed/{prefix}_R1_val_1.fq.gz"),
        r2 = path.join(config["output_dir"], "{species}/fastqs_trimmed/{prefix}_R2_val_2.fq.gz"),
    log:
        path.join(config["log_dir_run"], "{species}/{prefix}/2_trim.log")
    conda:
        "../envs/sci_rna.yaml"
    threads:
        config["threads"]["trim"]
    resources:
        mem_mib = lambda wildcards, threads: 1024 + threads * 512,  # 1GB + 512MB per thread
        runtime = 360,  # 6 hours
    shell:
        """
        trim_galore \\
            --paired {input.r1} {input.r2} \\
            -a2 AAAAAAAA \\
            --stringency 3 \\
            --cores {threads} \\
            -o {params.output_folder} &> {log}
        """
        
        # Record stats
        """
        sed -i '/reads_trimmed/d' {params.count_log} # remove reads_trimmed line if it exists
        # n reads after trimming
        nafter=$(($(zcat {output.r1} | wc -l) / 4))
        printf "$nafter\\treads_trimmed\\n" >> {params.count_log}
        """
