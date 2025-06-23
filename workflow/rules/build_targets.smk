from os import path

# Based on the config file, this rule generates a list of targets for building final objects.
# human or mouse in the config file and the barcodes are not empty.

checkpoint generate_final_targets:
    input:
        config["barcodes"]["ligation_barcode"]
    params:
        build_objects = " ".join(config["output_species"])
    output:
        config["build_objects"]
    threads: 1
    resources:
        mem_mib = 100,  # 100MB should be enough for a small file
        runtime = 5, # 5 minutes
    shell:
        """
        echo {params.build_objects} > {output}
        """

def generate_final_build_objects(wildcards):
    """
    Read the build_objects file and return a list of targets.

    human/mouse gene/exon files are determined based on the config settings.
    """
    # Run the checkpoint to generate tsv before reading the build objects
    build_objects_path = checkpoints.generate_final_targets.get().output

    targets = []

    with open(str(build_objects_path), 'r') as f:
        for species in f.read().strip().split():
            if path.exists(config["barcodes"]["rt_barcode"][species]) and \
                    path.getsize(config["barcodes"]["rt_barcode"][species]) > 0:
                targets.append(species)

    build_targets = []

    for species in targets:
        if not species:
            continue  # Skip empty species
        if config["output_type"]["gene"]:
            for file_key in config["output_matrix"]:
                build_targets.append(path.join(config["output_dir_final"], f"{species}_gene", config["output_matrix"][file_key]))
        if config["output_type"]["exon"]:
            for file_key in config["output_matrix"]:
                build_targets.append(path.join(config["output_dir_final"], f"{species}_exon", config["output_matrix"][file_key]))

    return build_targets
