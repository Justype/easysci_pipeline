from os import path

configfile: "config.yaml"

PREFIXES = [prefix for prefix in open(config["input_prefix_file"]).read().splitlines() if prefix]

rule count_gene:
    params:
        input_folder = lambda wildcards: path.join(config["output_dir"], wildcards.species, "bam_filtered"),
        output_folder = lambda wildcards: path.join(config["output_dir"], wildcards.species, "gene"),
        i7_prefix_file = config["input_prefix_file"],
    input:
        bams = lambda wildcards: expand(
            path.join(config["output_dir"], wildcards.species, "bam_filtered", "{prefix}.filtered.dedup.bam"),
            prefix=PREFIXES
        ),
        gtf = path.join(config["genome_folder"], "{species}.gtf"),
        rt_barcode_tsv = lambda wildcards: config["barcodes"]["rt_barcode_tsv"][wildcards.species],
    output:
        # counts = lambda wildcards: expand(
        #     path.join(config["output_dir"], wildcards.species, "gene", "{prefix}.count.gz"),
        #     prefix=PREFIXES
        # ),
        # cell_ids = lambda wildcards: expand(
        #     path.join(config["output_dir"], wildcards.species, "gene", "{prefix}_cell_ids.csv.gz"),
        #     prefix=PREFIXES
        # ),
        gene_ids = path.join(config["output_dir"], "{species}", "gene", "gene_ids.csv.gz"),
    log:
        path.join(config["log_dir_run"], "{species}/5_count_gene.log")
    conda:
        "../envs/sci_rna.yaml"
    threads:
        config["threads"]["count_gene"]
    resources:
        mem_mb = lambda wildcards, threads: 2048 + threads * 256,  # 2GB + 256MB per thread
        runtime = "12:00:00"
    shell:
        """
        python workflow/scripts/counting_gene_paired_parallel.py \
            --threads {threads} --gzip \
            --input_folder {params.input_folder} \
            --output_folder {params.output_folder} \
            --i7_prefix_file {params.i7_prefix_file} \
            --rt_barcode_tsv {input.rt_barcode_tsv} \
            --gtf {input.gtf} &> {log}
        """

rule merge_gene:
    params:
        input_folder = lambda wildcards: path.join(config["output_dir"], wildcards.species, "gene"),
        output_folder = path.join(config["output_dir_final"], "{species}_gene"),
        i7_prefix_file = config["input_prefix_file"],
    input:
        gene_ids = path.join(config["output_dir"], "{species}", "gene", "gene_ids.csv.gz"),
        rt_barcode_tsv = lambda wildcards: config["barcodes"]["rt_barcode_tsv"][wildcards.species],
    output:
        mtx = path.join(config["output_dir_final"], "{species}_gene", config["output_matrix"]["mtx"]),
        cells = path.join(config["output_dir_final"], "{species}_gene", config["output_matrix"]["cells"]),
        features = path.join(config["output_dir_final"], "{species}_gene", config["output_matrix"]["features"]),
    log:
        path.join(config["log_dir_run"], "{species}/6_merge_gene.log")
    conda:
        "../envs/sci_rna.yaml"
    threads: 1
    resources:
        mem_mb = 100 * len(PREFIXES),  # 100MB per i7 prefix
        runtime = "4:00:00",
    shell:
        """
        python workflow/scripts/merging_count_gene.py \
            --gzip \
            --input_folder {params.input_folder} \
            --output_folder {params.output_folder} \
            --i7_prefix_file {params.i7_prefix_file} \
            --rt_barcode_tsv {input.rt_barcode_tsv} &> {log}
        """
