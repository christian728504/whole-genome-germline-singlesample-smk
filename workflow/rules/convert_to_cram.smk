rule convert_to_cram:
    input:
        bam = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam",
    output:
        cram_md5 = f"{config.results_dir}/convert_to_cram/{{sample}}.cram.md5",
        cram = f"{config.results_dir}/convert_to_cram/{{sample}}.cram"
    container: config.environments.default
    log: f"{config.log_dir}/convert_to_cram/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1

        REF_CACHE_TMP=$(mktemp -d -p {config.tmp_dir} ref_cache.XXXXXX)
        trap "rm -rf $REF_CACHE_TMP" EXIT

        samtools view -@ {resources.threads} -C -T {input.fasta} {input.bam} | \
        md5sum | awk '{print $1}' > {output.cram_md5}

        # Create REF_CACHE. Used when indexing a CRAM
        seq_cache_populate.pl -root $REF_CACHE_TMP {input.fasta}
        export REF_PATH=:
        export REF_CACHE=$REF_CACHE_TMP/%2s/%2s/%s

        samtools index -@ {resources.threads} {output.cram}
        """

