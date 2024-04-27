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
        "results/local_adaptation_scenarios/asymmetric_fitness_variance.csv.gz", 
    log:
        "logs/local_adaptation_scenarios/asymmetric_fitness_variance.log",
    resources:
        mem_mb=3000,
        runtime=10,
    shell:
        "Rscript --vanilla {input.script} {input.files} {output} > {log} 2> {log}"
