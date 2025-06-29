from os import path

configfile: "config.yaml"

rule generate_index:
    params:
        prebuilt_star_index = lambda wildcards: config["prebuilt_star_index"][wildcards.species],
        index_length = config["seq_length"] - 1,
    input:
        fasta = path.join(config["genome_folder"], "{species}.fa"),
        gtf = path.join(config["genome_folder"], "{species}_annotation.gtf"),
    output:
        directory(path.join(config["star_index_folder"], "{species}"))
    log:
        path.join(config["log_dir_run"], "{species}/0_generate_index.log")
    conda:
        "../envs/sci_rna.yaml"
    threads:
        config["threads"]["generate_index"]
    resources:
        mem_mib = 40960,  # 40GB
        runtime = 720,  # 12 hours
    shell:
        """
        if [ -d {params.prebuilt_star_index} ]; then
            echo "Using prebuilt STAR index: {params.prebuilt_star_index} for {wildcards.species}" | tee {log}
            ln -s {params.prebuilt_star_index} {output} 2>/dev/null || cp -r {params.prebuilt_star_index} {output}
        else
            echo "Generating STAR index for {wildcards.species}" | tee {log}
            STAR --runThreadN {threads} \
                 --runMode genomeGenerate \
                 --genomeDir {output} \
                 --genomeFastaFiles {input.fasta} \
                 --sjdbGTFfile {input.gtf} \
                 --sjdbOverhang {params.index_length} &>> {log}
        fi
        """

rule create_fasta:
    params:
        prebuilt_star_index = lambda wildcards: config["prebuilt_star_index"][wildcards.species],
        fasta = lambda wildcards: path.join(config["fasta_file"][wildcards.species]),
        is_url = lambda wildcards: config["fasta_file"][wildcards.species].startswith("http"),
        is_gzipped = lambda wildcards: config["fasta_file"][wildcards.species].endswith(".gz"),
    output:
        path.join(config["genome_folder"], "{species}.fa"),
    log:
        path.join(config["log_dir_run"], "{species}/0_create_fasta.log")
    conda:
        "../envs/sci_rna.yaml"
    threads:
        config["threads"]["generate_index"]
    resources:
        mem_mib = 1024,  # 1GB
        runtime = 120,  # 2 hours
    shell:
        """
        if [ -d {params.prebuilt_star_index} ]; then
            echo "STAR index already exists: {params.prebuilt_star_index} for {wildcards.species}" | tee {log}
            touch {output}  # Create an empty file to satisfy the output requirement
        else
            # Check if params.fasta is url or local file and check if gzipped
            echo "IS URL: {params.is_url}, GZIPPED: {params.is_gzipped}" | tee -a {log}
            if [ "{params.is_url}" = "True" ]; then
                echo "Downloading FASTA file from URL: {params.fasta}" | tee -a {log}
                if [ "{params.is_gzipped}" = "True" ]; then
                    wget -qO- "{params.fasta}" | pigz -p {threads} -d > {output}
                else
                    wget -qO {output} "{params.fasta}"
                fi
            else
                if [ "{params.is_gzipped}" = "True" ]; then
                    echo "Unzipping gzipped FASTA file: {params.fasta}" | tee -a {log}
                    pigz -p {threads} -d < {params.fasta} > {output}
                else
                    echo "Linking FASTA file: {params.fasta}" | tee -a {log}
                    ln -s {params.fasta} {output} 2>/dev/null || cp {params.fasta} {output}
                fi
            fi
        fi
        """

rule create_gtf: # GTF cannot be prebuilt, so we always create it
    params:
        gtf = lambda wildcards: path.join(config["gtf_file"][wildcards.species]),
        is_url = lambda wildcards: config["gtf_file"][wildcards.species].startswith("http"),
        is_gzipped = lambda wildcards: config["gtf_file"][wildcards.species].endswith(".gz"),
    output:
        path.join(config["genome_folder"], "{species}_annotation.gtf"),
    log:
        path.join(config["log_dir_run"], "{species}/0_create_gtf.log")
    conda:
        "../envs/sci_rna.yaml"
    threads:
        config["threads"]["generate_index"]
    resources:
        mem_mib = 1024,  # 1GB
        runtime = 120,  # 2 hours
    shell:
        """
        # Check if input is url or local file and check if gzipped
        echo "IS URL: {params.is_url}, GZIPPED: {params.is_gzipped}" | tee -a {log}
        # Check if params.is_url is True or False
        if [ "{params.is_url}" = "True" ]; then
            echo "Downloading GTF file from URL: {params.gtf}" | tee -a {log}
            if [ "{params.is_gzipped}" = "True" ]; then
                wget -qO- "{params.gtf}" | pigz -p {threads} -d > {output}
            else
                wget -qO {output} "{params.gtf}"
            fi
        else
            if [ "{params.is_gzipped}" = "True" ]; then
                echo "Unzipping gzipped GTF file: {params.gtf}" | tee -a {log}
                pigz -p {threads} -d < {params.gtf} > {output}
            else
                echo "Linking GTF file: {params.gtf}" | tee -a {log}
                ln -s {params.gtf} {output} 2>/dev/null || cp {params.gtf} {output}
            fi
        fi
        """
