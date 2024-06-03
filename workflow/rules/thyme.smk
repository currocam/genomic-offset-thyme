rule gea:
  input:
    data = "steps/genomic_offset/genos27_03_2024_prepro.Rds",
    script = "scripts/07_gea_thyme.qmd"
  output:
    "results/GEA/genos27_03_2024.html"
  shell:
    """
    Rscript -e 'rmarkdown::render(\"{input.script}\", output_file = \"../{output}\")'
    """

rule genomic_offset_thyme:
    input:
        files=expand(
      "steps/slim/{model}_s{seed}.Rds",
      model=["thyme"],
      seed=range(100, 120)
      ),
        script="scripts/08-thyme_simulations.jl",
    output:
        "results/thyme_simulations/gradient_forest.csv",
    log:
        "logs/thyme_simulations/gradient_forest.log",
    resources:
        mem_mb=8000,
        runtime=420,
    threads: 1
    shell:
        "julia --project=. -t {threads} {input.script} {input.files} > {output} 2> {log}"

rule thyme_climate:
    input:
        files=expand(
      "steps/slim/{model}_s{seed}.Rds",
      model=["thymeclimate"],
      seed=range(100, 140)
      ),
        script="scripts/09_thyme_climate_change_gradient_forest.R",
    output:
        "results/thyme_simulations/climate_change_gradient_forest.csv",
    log:
        "logs/thyme_simulations/climate_change_gradient_forest.log",
    resources:
        mem_mb=10000,
        runtime=150,
    threads: 3
    shell:
        "Rscript --vanilla {input.script} {input.files} {output} > {log} 2> {log}"

rule causal_thyme:
    input:
        files=expand(
      "steps/slim/{model}_s{seed}.Rds",
      model=["thyme"],
      seed=range(100, 120)
      ),
        script="sandbox/causal_thyme.jl",
    output:
        "sandbox/causal_thyme.csv",
    log:
        "logs/thyme_simulations/causal_thyme.log",
    resources:
        mem_mb=8000,
        runtime=10,
    threads: 1
    shell:
        "julia --project=. {input.script} {input.files} > {output} 2> {log}"

rule empirical_thyme:
    input:
        files=expand(
      "steps/slim/{model}_s{seed}.Rds",
      model=["thyme"],
      seed=range(100, 120)
      ),
        script="sandbox/empirical_thyme.jl",
    output:
        "sandbox/empirical_thyme.csv",
    log:
        "logs/thyme_simulations/empirical_thyme.log",
    resources:
        mem_mb=8000,
        runtime=200,
    threads: 1
    shell:
        "julia --project=. {input.script} {input.files} > {output} 2> {log}"

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
        runtime=60,
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
        runtime=70,
    threads: 1
    shell:
        "python {input.script} {input.tree} > {output} 2> {log}"

rule thyme_msprime_analysis:
    input:
        script = "scripts/10_thyme_analysis.py",
        tree = "steps/msprime/thyme_{seed}_n{ticksAfter}.tree",
    output:
        "steps/msprime/thyme_{seed}_n{ticksAfter}_cors.csv"
    log:
        "logs/thyme_{seed}_n{ticksAfter}_cors.log",
    conda:
        "../envs/msprime_analysis.yaml"
    resources:
        mem_mb=8000,
        runtime=200,
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

rule thyme_msprime_several:
    input:
        expand("steps/msprime/thyme_{seed}_n{ticksAfter}_cors.csv", seed=range(500, 520), ticksAfter = [5, 20, 100]),
    output:
        "results/msprime/thyme_correlation.csv"
    shell:
        "echo 'file,method,climate_change_rate,generation_since,causal2empirical_mean,causal2empirical_std,causal2neglog_mean,causal2neglog_std,empirical2neglog_mean,empirical2neglog_std' > {output} && "
        "cat {input} >> {output}"


rule thyme_msprime_past_several:
    input:
        expand("steps/msprime/thyme_{seed}_n{ticksAfter}_past.csv", seed=range(500, 520), ticksAfter = [5, 20, 100]),
