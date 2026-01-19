rule collect_raw_wgs_metrics:
    input:
        bam = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam",
        reference_fasta = config.input_files.reference_fasta,
        reference_fasta_dict = config.input_files.reference_fasta_dict,
    output:
        metrics = f"{config.results_dir}/collect_raw_wgs_metrics/{{sample}}.raw_wgs_metrics",
    params:
        read_length = config.get("read_length", 250),
        compression_level = config.get("compression_level", 2),
        picard_tmp_dir = config.get("picard_tmp_dir", ".picard"),
    container: config.environments.default
    log: f"{config.log_dir}/collect_raw_wgs_metrics/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        
        PICARD_TMP=$(mktemp -d -p {params.picard_tmp_dir} picard_tmp.XXXXXX)
        trap "rm -rf $PICARD_TMP" EXIT

        java -Dsamjdk.compression_level={params.compression_level} \
        -Dorg.xerial.snappy.tempdir=$PICARD_TMP -Djava.io.tmpdir=$PICARD_TMP \
        -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m -jar /usr/bin/picard.jar \
        CollectRawWgsMetrics \
            -INPUT {input.bam} \
            -VALIDATION_STRINGENCY SILENT \
            -REFERENCE_SEQUENCE {input.reference_fasta} \
            -INCLUDE_BQ_HISTOGRAM true \
            -INTERVALS {input.reference_fasta_dict} \
            -OUTPUT {output.metrics} \
            -USE_FAST_ALGORITHM true \
            -READ_LENGTH {params.read_length}
        """
