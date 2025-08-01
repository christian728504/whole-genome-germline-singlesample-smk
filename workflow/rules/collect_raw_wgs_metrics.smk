rule collect_raw_wgs_metrics:
    input:
        bam = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam",
        reference_fasta = config.input_files.reference_fasta,
        wgs_coverage_interval_list = config.input_files.wgs_coverage_interval_list,
    output:
        metrics = f"{config.results_dir}/collect_raw_wgs_metrics/{{sample}}.wgs_metrics",
    params:
        read_length = config.get("read_length", 250),
    container: config.environments.default
    log: f"{config.log_dir}/collect_raw_wgs_metrics/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        java -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m -jar /usr/bin/picard.jar \
        CollectRawWgsMetrics \
            -INPUT {input_bam} \
            -VALIDATION_STRINGENCY SILENT \
            -REFERENCE_SEQUENCE {input.reference_fasta} \
            -INCLUDE_BQ_HISTOGRAM true \
            -INTERVALS {input.wgs_coverage_interval_list} \
            -OUTPUT {input.metrics} \
            -USE_FAST_ALGORITHM true \
            -READ_LENGTH {params.read_length}
        """
