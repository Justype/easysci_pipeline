from os import path

configfile: "config.yaml"

rule generate_py_barcode:
    params:
        base_dir = lambda wildcards, output: path.dirname(output.ligation_barcode)
    input:
        barcode_xlsx = config["barcodes_xlsx"]
    output:
        ligation_barcode = config["barcodes"]["ligation_barcode"],
        rt_barcode_human = config["barcodes"]["rt_barcode"]["human"],
        rt_barcode_mouse = config["barcodes"]["rt_barcode"]["mouse"],
        rt_barcode_tsv_human = config["barcodes"]["rt_barcode_tsv"]["human"],
        rt_barcode_tsv_mouse = config["barcodes"]["rt_barcode_tsv"]["mouse"],
    log:
        path.join(config["log_dir_run"], "generate_py_barcode.log")
    conda:
        "../envs/sci_rna.yaml"
    threads: 1
    resources:
        mem_mb = 500,
        runtime = "00:30:00"
    shell:
        """
        python workflow/scripts/generate_py_barcode.py {input.barcode_xlsx} {params.base_dir} &> {log}

        # Create empty files if they do not exist
        if [ ! -f {output.rt_barcode_human} ]; then
            touch {output.rt_barcode_human}
        fi
        if [ ! -f {output.rt_barcode_tsv_human} ]; then
            touch {output.rt_barcode_tsv_human}
        fi
        if [ ! -f {output.rt_barcode_mouse} ]; then
            touch {output.rt_barcode_mouse}
        fi
        if [ ! -f {output.rt_barcode_tsv_mouse} ]; then
            touch {output.rt_barcode_tsv_mouse}
        fi
        """
