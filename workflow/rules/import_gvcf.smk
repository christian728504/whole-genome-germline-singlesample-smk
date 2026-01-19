def get_split_interval(wildcards):
    scatter = int(wildcards.scatter)
    return f"{config.results_dir}/split_interval_list/{{split}}-scattered.interval_list "

rule import_gvcf:
    input:
        interval = get_split_interval,
        sample_map = rules.sample_map.output.sample_map,
    output:
        signal = f"{config.results_dir}/import_gvcf/{{split}}.continue",
        outdir = directory(f"{config.results_dir}/genomicsdb/"),
    params:
        batch_size = config.import_gvcf.batch_size,
        initial_heap_size = lambda resources: round(resources.mem_mb * (1/3), -3),
    container: config.environments.gatk
    log: f"{config.log_dir}/import_gvcf/{{split}}.log"
    shell:
        """
        exec >> {log} 2>&1

        gatk --java-options "-Xms{params.initial_heap_size}m -Xmx{resources.mem_mb}m" \
            GenomicsDBImport \
            --genomicsdb-workspace-path {output.outdir} \
            --batch-size {params.batch_size} \
            -L {input.interval} \
            --sample-name-map {input.sample_map} \
            --reader-threads {resources.threads} \
            --merge-input-intervals \
            --consolidate

        touch {output.signal}
        """