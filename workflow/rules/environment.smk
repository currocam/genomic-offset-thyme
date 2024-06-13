rule env_uncorrelated_env:
    input:
        files1=expand(
            "steps/slim/m4.1_s{seed}_nQTL1s{n}.Rds",
            seed=range(100, 110),
            n=[20, 50, 100],
        ),
        script="scripts/06_uncorrelated_env.R"
    output:
        "results/local_adaptation_scenarios/m4_offsets_uncorrelated.csv", 
    log:
        "logs/local_adaptation_scenarios/m4_offsets_uncorrelated.log",
    resources:
        mem_mb=10000,
        runtime=50,
    threads: 5
    shell:
        "Rscript --vanilla {input.script} {input.files1} {output}> {log} 2> {log}"

rule env_confounded_env:
    input:
        files=expand(
            "steps/slim/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.Rds",
            seed=range(100, 200),
            n=[100],
            model=["m4.3"],
        ),
        script="scripts/07_confounded_env.R"
    output:
        "results/local_adaptation_scenarios/m4_offsets_confounded.csv.gz", 
    log:
        "logs/local_adaptation_scenarios/m4_offsets_confounded.log",
    resources:
        mem_mb=10000,
        runtime=100,
    threads: 5
    shell:
        "Rscript --vanilla {input.script} {input.files} {output} > {log} 2> {log}"
