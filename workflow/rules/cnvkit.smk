__author__ = "Martin Rippin"
__copyright__ = "Copyright 2022, Martin Rippin"
__email__ = "martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_access:
    input:
        fasta=config["reference"]["fasta"],
    output:
        bed=temp("wgs_pon_grasnatter/cnvkit_access/access.bed"),
    params:
        extra=config.get("cnvkit_access", {}).get("extra", ""),
        min_gap_size=config.get("cnvkit_access", {}).get("min_gap_size", "10000"),
    log:
        "wgs_pon_grasnatter/cnvkit_access/access.bed.log",
    benchmark:
        repeat(
            "wgs_pon_grasnatter/cnvkit_access/access.bed.benchmark.tsv",
            config.get("cnvkit_access", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_access", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("cnvkit_access", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_access", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_access", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_access", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_access", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_access", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: generate cnvkit access file wgs_pon_grasnatter/cnvkit_access/access.bed"
    shell:
        "(cnvkit.py access "
        "{input.fasta} "
        "-s {params.min_gap_size} "
        "-o {output.bed} "
        "{params.extra}) &> {log}"


rule cnvkit_batch_pon:
    input:
        bam=expand(
            "parabricks/fq2bam/{sample}_N.bam",
            sample=samples.index,
        ),
        bai=expand(
            "parabricks/fq2bam/{sample}_N.bam.bai",
            sample=samples.index,
        ),
        fasta=config["reference"]["fasta"],
        bed="wgs_pon_grasnatter/cnvkit_access/access.bed",
    output:
        cnn="wgs_pon_grasnatter/cnvkit_batch_pon/pon.cnn",
    params:
        extra=config.get("cnvkit_batch_pon", {}).get("extra", ""),
    log:
        "wgs_pon_grasnatter/cnvkit_batch_pon/pon.cnn.log",
    benchmark:
        repeat(
            "wgs_pon_grasnatter/cnvkit_batch_pon/pon.cnn.benchmark.tsv",
            config.get("cnvkit_batch_pon", {}).get("benchmark_repeats", 1),
        )
    resources:
        threads=config.get("cnvkit_batch_pon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_batch_pon", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_batch_pon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_batch_pon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_batch_pon", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_batch_pon", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: generate pon for cnvkit wgs_pon_grasnatter/cnvkit_batch_pon/pon.cnn"
    shell:
        "(cnvkit.py batch "
        "-m wgs "
        "-n {input.bam} "
        "-f {input.fasta} "
        "-g {input.bed} "
        "-p {threads} "
        "-d wgs_pon_grasnatter/cnvkit_batch_pon/ "
        "--output-reference {output.cnn} "
        "{params.extra}) &> {log}"
