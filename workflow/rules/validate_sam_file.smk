def get_skip_mate_validation_flag(wildcards):
    # The contents of this file will already be formatted as a string for use in the shell command template
    skip_mate_validation_file = Path(f"{config.results_dir}/check_prevalidation/{wildcards.sample}.skip_mate_validation.txt")
    with skip_mate_validation_file.open("r") as f:
        skip_mate_validation = f.read().strip()
    return skip_mate_validation

rule validate_sam_file:
    input:
        bam = f"{config.results_dir}/convert_to_cram/{{sample}}.bam",
        reference_fasta = config.input_files.reference_fasta,
    output:
        validation_report = f"{config.results_dir}/validate_sam_file/{{sample}}.validation_report",
    params:
        skip_mate_validation = get_skip_mate_validation_flag,
    container: config.environments.default
    log: f"{config.log_dir}/validate_sam_file/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        java -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m -jar /usr/bin/picard.jar \
        ValidateSamFile \
            -INPUT {input_bam} \
            -OUTPUT {output.validation_report} \
            -REFERENCE_SEQUENCE {input.reference_fasta} \
            -MAX_OUTPUT 1000000000 \
            -IGNORE MISSING_TAG_NM \
            -MODE VERBOSE \
            {params.skip_mate_validation} \
            -IS_BISULFITE_SEQUENCED false
        """
