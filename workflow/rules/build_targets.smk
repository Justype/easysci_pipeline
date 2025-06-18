from os import path

def has_human_sample():
    """
    Check if config["barcodes"]["rt_barcode"]["human"] file exists and is not empty.
    Returns True if it exists, False otherwise.
    """
    return path.exists(config["barcodes"]["rt_barcode"]["human"]) and \
               path.getsize(config["barcodes"]["rt_barcode"]["human"]) > 0

def has_mouse_sample():
    """
    Check if config["barcodes"]["rt_barcode"]["mouse"] file exists and is not empty.
    Returns True if it exists, False otherwise.
    """
    return path.exists(config["barcodes"]["rt_barcode"]["mouse"]) and \
           path.getsize(config["barcodes"]["rt_barcode"]["mouse"]) > 0

def create_build_objects():
    """
    Create the build_objects file based on the presence of human and mouse samples.
    Returns the string of build objects to be used in the pipeline.

    If both samples are present, it returns "human mouse".
    """
    build_objects = ""
    if has_human_sample():
        # print("Human sample detected.") # Test
        build_objects += "human "
    if has_mouse_sample():
        # print("Mouse sample detected.") # Test
        build_objects += "mouse "

    return build_objects.strip()

checkpoint generate_final_targets:
    input:
        config["barcodes"]["ligation_barcode"]
    params:
        build_objects = create_build_objects() # in rules/generate_py_barcode.py
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

    with open(str(build_objects_path), 'r') as f:
        targets = f.read().strip().split()

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
