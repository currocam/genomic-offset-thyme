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

rule geometric_offset_thyme:
    input:
        files=expand(
      "steps/slim/{model}_s{seed}.Rds",
      model=["thyme", "thymeclimate"],
      seed=range(100, 120)
      ),
        script="scripts/08-thyme_geometric_offset_simulations.R",
    output:
        "results/thyme_simulations/geometric_genomic_offset.csv.gz",
    log:
        "logs/thyme_simulations/geometric_genomic_offset.log",
    resources:
        mem_mb=20000,
        runtime=400,
    threads: 5
    shell:
        "Rscript --vanilla {input.script} {input.files} {output} > {log} 2> {log}"

rule thyme_simulation:
  input:
    expand(
      "steps/slim/{model}_s{seed}.Rds",
      model=["thyme", "thymeclimate"],
      seed=range(100, 120)
      )