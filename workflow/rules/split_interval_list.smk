rule split_interval_list:
    input:
        interval_list = config.input_files.split_interval_list\
        ref_fasta = config.input_files.reference_fasta,
    output:
        outdir = directory(f"{config.results_dir}/split_interval_list"),
    params:
        scatter_count = SPLIT_INTERVAL_SCATTER_COUNT,
        scatter_mode = config.split_interval_list.scatter_mode,
    container: config.environments.gatk
    log: f"{config.log_dir}/split_interval_list.log"
    shell:
        """
        exec >> {log} 2>&1

        gatk --java-options "-Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m" SplitIntervals \
            -L {input.interval_list} -O {output.outdir} -scatter {params.scatter_count} -R {input.ref_fasta} \
            -mode {params.scatter_mode} --interval-merging-rule OVERLAPPING_ONLY
        """