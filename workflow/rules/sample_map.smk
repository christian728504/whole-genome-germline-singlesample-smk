rule sample_map:
    input:
        gvcfs = expand(f"{config.results_dir}/reblock/{{sample}}.rb.g.vcf.gz", sample=SAMPLES),
    output:
        sample_map = f"{config.results_dir}/sample_map.smap",
    container: config.environments.gatk
    log: f"{config.log_dir}/sample_map.log"
    run:
        import os, glob

        gvcf_files_abs = [os.path.abspath(x) for x in input.gvcfs]
        records = list(zip(SAMPLES, gvcf_files_abs))
        assert len(records) == len(SAMPLES), "Number of samples does not match number of gVCF files."
        assert len(set(SAMPLES)) == len(SAMPLES), "Samples are not unique."
        assert all(os.path.exists(x) for x in gvcf_files_abs), "gVCF files do not exist."

        with open(output.sample_map, "w") as f:
            for record in records:
                f.write(f"{record[0]}\t{record[1]}\n")