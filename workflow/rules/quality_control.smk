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
        mem_mb=16000,
        runtime=300,
    shell:
        "Rscript --vanilla {input.script} > {log} 2> {log}"

rule check_ld_quality_control:
    input:
        files=expand(
            "steps/pairwise_ld/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.ld",
            seed=range(100, 110),
            n=[20, 50, 100],
            model=["m3.2"],
        )

rule parse_slim_with_fst:
    input:
        "steps/slim/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.vcf",
        "steps/slim/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.txt",
    output:
        "steps/slim/{model}_s{seed}_nQTL1s{n}__nQTL2s{n}_fst.Rds"
    resources:
        mem_mb=1000,
        runtime=10,
    wildcard_constraints:
        seed=r"\d+",
        n=r"\d+",
    params:
        script="workflow/scripts/parse_slim_with_fst.jl",
    shell:
        """
        julia  --project=. {params.script} {input} {output}
        """

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
        mem_mb=3000,
        runtime=100,
    shell:
        "Rscript --vanilla {input.script} > {log} 2> {log}"

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
        mem_mb=3000,
        runtime=100,
    shell:
        "Rscript --vanilla {input.script} {input.files} {output} > {log} 2> {log}"

rule asymmetric_fitness_variance:
    input: 
        files = expand(
            "steps/slim/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.Rds",
            seed=range(100, 106),
            n=[20, 50, 100],
            model=["m3.3", "m3.4"],
        ),
        script="scripts/04_assymetric_traits.R"
    output:
        "results/local_adaptation_scenarios/m3_offsets_asymmetric.csv", 
    log:
        "logs/local_adaptation_scenarios/m3_offsets_asymmetric.log",
    resources:
        mem_mb=3000,
        runtime=10,
    shell:
        "Rscript --vanilla {input.script} > {log} 2> {log}"
