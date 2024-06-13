rule quality_control:
    input:
        files=expand(
            "steps/slim/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.Rds",
            seed=range(100, 110),
            n=[20, 50, 100],
            model=["m3.1", "m3.2", "m3.5"],
        ),
        script="scripts/01_local_adaptation_scenarios.R",
    output:
        "results/local_adaptation_scenarios/quality_control.csv.gz",
    log:
        "logs/local_adaptation_scenarios/m3_offsets.log",
    resources:
        mem_mb=20000,
        runtime=150,
    threads: 2,
    shell:
        "Rscript --vanilla {input.script} {input.files} {output} > {log} 2> {log}"

rule only_local_foreign:
    input:
        files=expand(
            "steps/slim/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.Rds",
            seed=range(100, 103),
            n=[20, 50, 100],
            model=["m3.5"],
        ),
        script="scripts/01_local_adaptation_scenarios.R",
    output:
        "foreign_local.csv.gz",
    log:
        "foreign_local.log",
    resources:
        mem_mb=16000,
        runtime=400,
    shell:
        "Rscript --vanilla {input.script} {input.files} {output} > {log} 2> {log}"

rule spurious_genomic_offset:
    input:
        script="scripts/11-spurious_genomic_offset.R",
    output:
        "results/local_adaptation_scenarios/spurious_genomic_offsets.csv.gz",
    log:
        "local_adaptation_scenarios/spurious_genomic_offsets.log",
    resources:
        mem_mb=8000,
        runtime=15,
    shell:
        "Rscript --vanilla {input.script} {output} > {log} 2> {log}"
