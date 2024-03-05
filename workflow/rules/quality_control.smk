rule sims:
    input:
        expand(
            "steps/slim/{model}_s{seed}_nQTL1s{n}_nQTL2s_{n}.Rds",
            seed=range(100, 103),
            n=[20, 50, 100], model = ["m3.1", "m3.2"]
        )