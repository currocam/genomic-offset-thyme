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

rule quality_control_boots:
    input:
        files=expand(
             "steps/slim/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.Rds",
             seed=range(100, 110),
             n=[20, 100],
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
            model=["m3.1", "m3.2"],
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
