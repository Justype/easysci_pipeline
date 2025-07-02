from os import path

configfile: "config.yaml"

PREFIXES = [prefix for prefix in open(config["input_prefix_file"]).read().splitlines() if prefix]

rule create_rmsk_gtf:
    # Download RepeatMasker GTF file from UCSC
    params:
        rmsk = lambda wildcards: path.join(config["rmsk_file"][wildcards.species]),
        is_url = lambda wildcards: config["rmsk_file"][wildcards.species].startswith("http"),
        is_gzipped = lambda wildcards: config["rmsk_file"][wildcards.species].endswith(".gz"),
    output:
        path.join(config["genome_folder"], "{species}_rmsk.gtf"),
    log:
        path.join(config["log_dir_run"], "{species}/0_create_rmsk.log")
    conda:
        "../envs/sci_rna.yaml"
    threads: 1
    resources:
        mem_mib = 512,  # 512MB
        runtime = 120,  # 2 hours
    shell:
        """
        # Check if input is url or local file and check if gzipped
        echo "IS URL: {params.is_url}, GZIPPED: {params.is_gzipped}" | tee -a {log}
        # Check if params.is_url is True or False
        if [ "{params.is_url}" = "True" ]; then
            wget -qO {output} "{params.rmsk}"
        else
            # if [ "{params.is_gzipped}" = "True" ]; then
            #     echo "Unzipping gzipped RMSK GTF file: {params.rmsk}" | tee -a {log}
            #     pigz -p {threads} -d < "{params.rmsk}" > {output}
            # else
                echo "Linking RMSK GTF file: {params.rmsk}" | tee -a {log}
                ln -s "{params.rmsk}" {output} 2>/dev/null || cp "{params.rmsk}" {output}
            # fi
        fi
        """

rule bam_add_bc:
    # Add cell barcode and UMI to BAM files
    params:
        prefix = lambda wildcards: wildcards.prefix,
    input:
        bam = path.join(config["output_dir"], "{species}/bam/{prefix}.Aligned.out.bam"),
        rt_barcode_tsv = lambda wildcards: config["barcodes"]["rt_barcode_tsv"][wildcards.species],
    output:
        path.join(config["output_dir"], "{species}/bam_bc/{prefix}.sorted.bam"),
    log:
        path.join(config["log_dir_run"], "{species}/{prefix}/5_sort_bam_add_bc.log")
    conda:
        "../envs/sci_rna.yaml"
    threads:
        config["threads"]["velocyto"]
    resources:
        mem_mib = lambda wildcards, threads: threads * 850,  # samtools uses 768M per thread by default
        runtime = 360, # 6 hours
    shell:
        """
        samtools sort -@ {threads} -o - {input.bam} | \\
            python workflow/scripts/bam_add_cb_ub.py \\
                --i7_prefix {params.prefix} \\
                --rt_barcode_tsv {input.rt_barcode_tsv} \\
            > {output} 2> {log}
        """

rule velocyto:
    params:
        output_folder = lambda wildcards: path.join(config["output_dir"], wildcards.species, "velocyto", wildcards.prefix),
        cell_sorted_bam = lambda wildcards: path.join(config["output_dir"], wildcards.species, "bam_bc", f"cellsorted_{wildcards.prefix}.sorted.bam"),
    input:
        bam = path.join(config["output_dir"], "{species}/bam_bc/{prefix}.sorted.bam"),
        gtf = path.join(config["genome_folder"], "{species}_annotation.gtf"),
        rmsk = path.join(config["genome_folder"], "{species}_rmsk.gtf"),
    output:
        path.join(config["output_dir"], "{species}/velocyto/{prefix}.loom"),
    log:
        path.join(config["log_dir_run"], "{species}/{prefix}/6_velocyto.log")
    conda:
        "../envs/sci_rna.yaml"
    threads:
        config["threads"]["velocyto"]
    resources:
        mem_mib = lambda wildcards, threads: threads * 850,  # velocyto uses 768M per thread by default
        runtime = 360, # 6 hours
    shell:
        """
        velocyto run -@ {threads} \\
            -o {params.output_folder} \\
            -m {input.rmsk} \\
            {input.bam} \\
            {input.gtf} &> {log}

        # Move the output loom file to the expected output location
        mv {params.output_folder}/*.loom {output} &>> {log}

        # Remove Intermediate output bam file
        rm {params.cell_sorted_bam} 2>/dev/null || true

        rmdir {params.output_folder} 2>/dev/null || true
        """

rule merge_velocyto:
    input:
        loom_files = lambda wildcards: expand(
            path.join(config["output_dir"], wildcards.species, "velocyto", "{prefix}.loom"),
            prefix=PREFIXES
        ),
    output:
        path.join(config["output_dir_final"], "{species}.loom"),
    log:
        path.join(config["log_dir_run"], "{species}/9_merge_velocyto.log")
    conda:
        "../envs/sci_rna.yaml"
    threads: 1
    resources:
        mem_mib = 2048,  # 2GB
        runtime = 120,  # 2 hours
    shell:
        """
        python workflow/scripts/merging_velocyto.py \\
            --input {input.loom_files} \\
            --output {output} &> {log}
        """
