using Pkg
Pkg.activate(".")
using SnpArrays, CSV, DataFrames, RCall
R"library(tidyverse)"

function plink(bedfile)
  mat = convert(Matrix{Float64}, SnpArray(bedfile); impute=true)
  # Remove extension
  basename = replace(bedfile, r"\.bed$" => "");
  bim = DataFrame(CSV.File("$basename.bim"; header=false))
  fam = DataFrame(CSV.File("$basename.fam"; header=false))
  sample_names = fam[:, 2]
  snp_names = bim[:, 2]
  return mat, sample_names, snp_names
end

function run(prunned, unprunned, annotation_tsv, sample_tsv, environment_rds, outfile)
  prunned, prunned_samples, prunned_snps = plink(prunned)
  unprunned, unprunned_samples, unprunned_snps = plink(unprunned)
  gene_ann = DataFrame(
    CSV.File(annotation_tsv; header=["Chromosome", "Position", "Id", "Gene"])
  )
  metadata = DataFrame(CSV.File(sample_tsv))
  environment = leftjoin(metadata, rcopy(R"read_rds($environment_rds)"), on=:site_id)
  select!(environment, Not(names(metadata)))

  @rput prunned prunned_samples prunned_snps unprunned unprunned_samples unprunned_snps gene_ann metadata environment
  R"""
  rownames(prunned) <- prunned_samples
  colnames(prunned) <- prunned_snps
  rownames(unprunned) <- unprunned_samples
  colnames(unprunned) <- unprunned_snps
  res <- list(
    prunned = prunned,
    unprunned = unprunned,
    gene_ann = gene_ann,
    metadata = metadata,
    environment = environment,
    info = "prunned contains SNPs prunned by LD R2 > 0.1. Gene annotation is for the unprunned genome. Metadata is info for the samples, and environment is the environmental data from Thomas."
  )
  write_rds(res, $outfile)
  """
end

run(
  snakemake.input["prunned"],
  snakemake.input["unprunned"],
  snakemake.input["annotation_tsv"],
  snakemake.input["sample_tsv"],
  snakemake.input["environment_rds"],
  snakemake.output[1]
)