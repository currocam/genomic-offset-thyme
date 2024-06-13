rule thyme_msprime_simulation:
    input:
        slim = "src/slim/thyme_msprime.slim",
        Rscript = "src/R/landscape_ac.R",
        script = "scripts/10_thyme_simulation.py",
    output:
        "steps/msprime/thyme_{seed}_n{ticksAfter}.tree",
    log:
        "logs/thyme_simulations/msprime_{seed}_n{ticksAfter}.log",
    wildcard_constraints:
        seed=r"\d+",
        ticksAfter=r"\d+",
    params:
        sigmaFitness = 1.0,
    conda:
        "../envs/msprime.yaml"
    resources:
        mem_mb=10000,
        runtime=200,
    threads: 1
    shell:
        "python {input.script} {input.slim} {input.Rscript} {params.sigmaFitness} {wildcards.ticksAfter} {wildcards.seed} > {output} 2> {log}"

rule thyme_msprime_past:
    input:
        script = "scripts/10_thyme_past.py",
        tree = "steps/msprime/thyme_{seed}_n{ticksAfter}.tree",
    output:
        "steps/msprime/thyme_{seed}_n{ticksAfter}_past.csv"
    log:
        "logs/thyme_{seed}_n{ticksAfter}_past.log",
    conda:
        "../envs/msprime_analysis.yaml"
    resources:
        mem_mb=8000,
        runtime=45,
    threads: 1
    shell:
        "python {input.script} {input.tree} > {output} 2> {log}"

rule thyme_msprime_tuning:
    input:
        script = "scripts/10_thyme_tuning.py",
        tree = "steps/msprime/thyme_{seed}_n{ticksAfter}.tree",
    output:
        "steps/msprime/thyme_{seed}_n{ticksAfter}_tuning.csv"
    log:
        "logs/thyme_{seed}_n{ticksAfter}_tuning.log",
    conda:
        "../envs/msprime_analysis.yaml"
    resources:
        mem_mb=8000,
        runtime=20,
    threads: 1
    shell:
        "python {input.script} {input.tree} > {output} 2> {log}"

rule thyme_msprime_tuning_output:
    input:
        expand("steps/msprime/thyme_{seed}_n{ticksAfter}_tuning.csv", seed=range(500, 520), ticksAfter = [5, 20, 100])
    output:
        "results/msprime/thyme_tuning.csv"
    log:
        "logs/thyme_tuning.log",
    shell:
        # Print header first and then concatenate all files
        "head -n 1 {input[0]} > {output} && "
        "tail -n +2 -q {input} >> {output}"

rule thyme_msprime_past_several:
    input:
        expand("steps/msprime/thyme_{seed}_n{ticksAfter}_past.csv", seed=range(500, 520), ticksAfter = [5, 20, 100]),
    output:
        "results/msprime/thyme_past.csv.gz"
    shell:
        "echo 'file,fraction,category,dataset,type,statistic,pvalue,r2adj' > results/msprime/thyme_past.csv && "
        "tail -n +2 -q {input} >> results/msprime/thyme_past.csv && "
        "gzip results/msprime/thyme_past.csv"

rule thyme_generation_time:
    input:
        infiles = expand("steps/msprime/thyme_{seed}_n{ticksAfter}.tree", seed=range(500, 520), ticksAfter = [5, 20, 100]),
        script = "scripts/10_thyme_generation.py",
    output:
        "results/msprime/generation_time.npy"
    log:
        "logs/msprime_generation_time.log",
    conda:
        "../envs/msprime_analysis.yaml"
    resources:
        mem_mb=8000,
        runtime=20,
    shell:
        "python {input.script} {input.infiles} {output} > {log} 2> {log}"

rule thyme_bootstrap:
    input:
        infiles = "steps/msprime/thyme_500_n20.tree",
        script = "scripts/12_thyme_past_bootstrap.py",
    output:
        "results/msprime/thyme_500_n20_bootstrap.csv"
    log:
        "logs/thyme_500_n20_bootstrap.log",
    conda:
        "../envs/msprime_analysis.yaml"
    resources:
        mem_mb=8000,
        runtime=60,
    threads: 5
    shell:
        "python -X juliacall-handle-signals=yes -X juliacall-threads={threads} {input.script} {input.infiles} {output} > {log} 2> {log}"
