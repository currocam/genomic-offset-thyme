rule quality_control:
    input:
        files=expand(
            "steps/slim/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.Rds",
            seed=range(100, 110),
            n=[20, 50, 100],
            model=["m3.1", "m3.2"],
        ),
        script="scripts/02_local_adaptation_scenarios.R",
    output:
        "results/local_adaptation_scenarios/m3_offsets.csv",
    log:
        "logs/local_adaptation_scenarios/m3_offsets.log",
    resources:
        mem_mb=3000,
        runtime=10,
    shell:
        "Rscript --vanilla {input.script} > {log}"

rule check_ld_quality_control:
    input:
        files=expand(
            "steps/pairwise_ld/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.ld",
            seed=range(100, 110),
            n=[20, 50, 100],
            model=["m3.2"],
        ),