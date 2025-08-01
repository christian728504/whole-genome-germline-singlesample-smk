rule collect_variant_calling_metrics:
    input:
        gvcf = f"{config.results_dir}/reblock/{{sample}}.rb.g.vcf.gz",
        dbsnp_vcf = config.input_files.dbsnp_vcf,
        reference_fasta_dict = config.input_files.reference_fasta_dict,
        evaluation_interval_list = config.input_files.evaluation_interval_list,
    output:
        variant_calling_detail_metrics = f"{config.results_dir}/collect_variant_calling_metrics/{{sample}}.variant_calling_detail_metrics",
        variant_calling_summary_metrics = f"{config.results_dir}/collect_variant_calling_metrics/{{sample}}.variant_calling_summary_metrics",
    params:
        ouput_prefix = f"{config.results_dir}/collect_variant_calling_metrics/{{sample}}",
    container:
    log: f"{config.log_dir}/collect_variant_calling_metrics/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        java -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m -jar /usr/bin/picard.jar \
        CollectVariantCallingMetrics \
            -INPUT {input.gvcf} \
            -OUTPUT {params.output_prefix} \
            -DBSNP {input.dbsnp_vcf} \
            -SEQUENCE_DICTIONARY {input.reference_fasta_dict} \
            -TARGET_INTERVALS {input.evaluation_interval_list} \
            -GVCF_INPUT true
        """
