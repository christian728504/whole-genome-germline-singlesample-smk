rule collect_quality_yield_metrics:
    input:
        bam = f"{config.results_dir}/unmapped_bam/{{sample}}.bam"
    output:
        metrics = f"{config.results_dir}/collect_quality_yield_metrics/{{sample}}.metrics"
    container: config.environments.default
    log: f"{config.log_dir}/collect_quality_yield_metrics/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        java -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m -jar /usr/bin/picard.jar CollectQualityYieldMetrics \
            -INPUT {input.bam} \
            -OQ true \
            -OUTPUT {output.metrics}
        """
