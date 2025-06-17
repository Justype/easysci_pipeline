from os import path

configfile: "config.yaml"

rule barcoding:
    params:
        input_prefix = lambda wildcards: path.join(config["input_dir"], wildcards.prefix),
        output_prefix = lambda wildcards: path.join(config["output_dir"], wildcards.species, "fastqs", wildcards.prefix),
        count_log = lambda wildcards: path.join(config["log_dir"], wildcards.species, f"{wildcards.prefix}.log"),
        barcode_pattern = lambda wildcards: ",".join([
            str(config["read_barcode_format"]["ligation_cell_barcode_length"]),
            str(config["read_barcode_format"]["rt_umi_length"]),
            str(config["read_barcode_format"]["rt_cell_barcode_length"]),
            str(config["read_barcode_format"]["polyT_length"])
        ])
    input:
        r1 = path.join(config["input_dir"], "{prefix}_R1.fastq.gz"),
        r2 = path.join(config["input_dir"], "{prefix}_R2.fastq.gz"),
        r3 = path.join(config["input_dir"], "{prefix}_R3.fastq.gz"),
        ligation_barcode = config["barcodes"]["ligation_barcode"],
        rt_barcode = lambda wildcards: config["barcodes"]["rt_barcode"][wildcards.species],
        rt_barcode_tsv = lambda wildcards: config["barcodes"]["rt_barcode_tsv"][wildcards.species],
    output:
        r1 = path.join(config["output_dir"], "{species}/fastqs/{prefix}_R1.fastq.gz"),
        r2 = path.join(config["output_dir"], "{species}/fastqs/{prefix}_R2.fastq.gz"),
    log:
        path.join(config["log_dir_run"], "{species}/{prefix}/1_barcoding.log")
    conda:
        "../envs/sci_rna.yaml"
    threads: 1
    resources:
        mem_mb = 500, 
        runtime = 240,  # 4 hours
    shell:
        """
        python workflow/scripts/barcoding_paired.py \\
            --input_prefix {params.input_prefix} \\
            --output_prefix {params.output_prefix} \\
            --ligation_barcode {input.ligation_barcode} \\
            --rt_barcode {input.rt_barcode} \\
            --rt_barcode_tsv {input.rt_barcode_tsv} \\
            --barcodes_pattern "{params.barcode_pattern}" \\
            --min_length 20 \\
            --log {params.count_log} &> {log}
        """
