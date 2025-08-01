def get_bamout_bams(wildcards):
    return expand(
        bam = f"{config.results_dir}/sort_bamout/{wildcards.sample}/{scatter}.bam",
        scatter=range(1, config.scatter_interval_list.scatter_count + 1)),
    )

def reformat_bamout_bams(wildcards):
    reformatted_bams = ""
    for scatter in range(1, config.scatter_interval_list.scatter_count + 1):
        reformatted_bams += f"{config.results_dir}/sort_bamout/{wildcards.sample}/{scatter}.bam "
    return reformatted_bams

rule merge_bamout:
    input:
        bams = get_bamout_bams
    output:
        bam = f"{config.results_dir}/merge_bamout/{{sample}}.bam",
        bam_index = f"{config.results_dir}/merge_bamout/{{sample}}.bam.bai",
    params:
        bams = reformat_bamout_bams
    container:
    log: f"{config.log_dir}/merge_bamout/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        samtools merge -@ {resources.threads} -m 2G -o {output.bam} {params.bams}
        samtools index -@ {resources.threads} {output.bam}
        """
