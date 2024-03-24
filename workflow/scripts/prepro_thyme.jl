using Pkg
Pkg.activate(".")
using SnpArrays, CSV, DataFrames, RCall
R"library(tidyverse)"
function run(bedfile, annotation_tsv, sample_tsv, environment_rds, outfile)
  genome = convert(Matrix{Float64}, SnpArray(bedfile); impute=true)

  gene_ann = DataFrame(
    CSV.File(annotation_tsv; header=["Chromosome", "Position", "Id", "Gene"])
  )
  @assert(size(genome, 2) == size(gene_ann, 1))

  metadata = DataFrame(CSV.File(sample_tsv))
  @assert(size(genome, 1) == size(metadata, 1))

  environment = innerjoin(metadata, rcopy(R"read_rds($environment_rds)"), on=:site_id)
  select!(environment, Not(names(metadata)))
  @assert(size(genome, 1) == size(environment, 1))

  @rput genome
  @rput gene_ann
  @rput metadata
  @rput environment

  R"""
  res <- list(
    genome = genome,
    gene_ann = gene_ann,
    metadata = metadata,
    environment = environment
  )
  write_rds(res, $outfile)
  """
end

run(
  snakemake.input["bedfile"],
  snakemake.input["annotation_tsv"],
  snakemake.input["sample_tsv"],
  snakemake.input["environment_rds"],
  snakemake.output[1]
)