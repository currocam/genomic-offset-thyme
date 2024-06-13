rule quality_control_boots:
    input:
        files=expand(
             "steps/slim/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.Rds",
             seed=range(100, 110),
             n=[50],
             model=["m3.2"],
        ),
        script="scripts/02_local_adaptation_scenarios_bootstraps.jl",
    output:
        "results/local_adaptation_scenarios/quality_control_bootstrapping.csv.gz",
    log:
        "logs/local_adaptation_scenarios/quality_control_bootstrapping.log",
    threads: 10
    resources:
        mem_mb=16000,
        runtime=350,
    shell:
        "julia -t {threads} --project=. {input.script} {input.files} {output} > {log} 2> {log}"


rule boot_regressions:
    input:
        files=expand(
            "steps/slim/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.Rds",
            seed=range(102, 105),
            n=[100],
            model=["m3.2"],
        ),
        script="scripts/03_fitness_offset_regression_bootstraps.jl",
    output:
        "results/local_adaptation_scenarios/regression_bootstraps.csv.gz",
    log:
        "logs/local_adaptation_scenarios/regression_bootstraps.log",
    threads: 3
    resources:
        mem_mb=16000,
        runtime=400,
    shell:
        "julia -t {threads} --project=. {input.script} {input.files} {output} > {log} 2> {log}"
