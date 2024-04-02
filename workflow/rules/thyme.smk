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