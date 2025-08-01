rule calculate_readgroup_checksum::
    input:
        bam = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam",
    output:
        md5 = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam.md5",
    container: config.environments.default
    log: f"{config.log_dir}/calculate_readgroup_checksum/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        java -Xms{resources.mem_mb}m -Xmx{reources.mem_mb}m -jar /usr/bin/picard.jar \
        CalculateReadGroupChecksum \
            INPUT= {input.bam} \
            OUTPUT= {output.md5}
        """
