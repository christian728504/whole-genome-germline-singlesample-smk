rule collect_unsorted_readgroup_bam_quality_metrics:
    input: 
        bam = f"{config.results_dir}/unmapped_bam_to_aligned/{{sample}}.bam"
    output:
        insert_size_metrics = f"{config.results_dir}/collect_unsorted_readgroup_bam_quality_metrics/{{sample}}.insert_size_metrics",
        insert_size_histogram = f"{config.results_dir}/collect_unsorted_readgroup_bam_quality_metrics/{{sample}}.insert_size_histogram.pdf"
    params:
        tmp_dir = config.tmp_dir,
        results_dir = config.results_dir
    container: config.environments.default
    log : f"{config.log_dir}/collect_unsorted_readgroup_bam_quality_metrics/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        java -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m -jar /usr/bin/picard.jar CollectMultipleMetrics \
            -INPUT {input.bam} \
            -OUTPUT {params.results_dir}/collect_unsorted_readgroup_bam_quality_metrics/{wildcards.sample} \
            -ASSUME_SORTED true \
            -PROGRAM null \
            -PROGRAM CollectBaseDistributionByCycle \
            -PROGRAM CollectInsertSizeMetrics \
            -PROGRAM MeanQualityByCycle \
            -PROGRAM QualityScoreDistribution \
            -METRIC_ACCUMULATION_LEVEL null \
            -METRIC_ACCUMULATION_LEVEL ALL_READS
        """

