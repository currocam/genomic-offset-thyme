library(tidyverse)
library(LEA)
library(conflicted)
source("src/R/gain_offset.R")

# Parses non local and local adapted simulations
parse_simulation <- function(non_local_file, local_file) {
  non_localadapt <- read_rds(non_local_file)
  localadapt <- read_rds(local_file)
  causal_data  <- tibble(
    scenario = c("Local adaptation", "Non local adaptation"),
    inference = "causal",
    simulation = list(localadapt, non_localadapt),
  ) |>
    rowwise() |>
    mutate(
      shifted_fitness = list(simulation[["Future fitness"]]),
      loci = list(c(simulation[["Index QTLs 1"]], simulation[["Index QTLs 2"]])),
      Y = list(simulation$Genotype[,loci]),
      X = list(matrix(
        c(simulation[["Current env 1"]], simulation[["Current env 2"]]),
        ncol=2)),
      Xstar = list(matrix(
        c(simulation[["Future env 1"]], simulation[["Future env 2"]]),
        ncol=2))
    )
  empirical_data  <- tibble(
    scenario = c("Local adaptation", "Non local adaptation"),
    inference = "empirical",
    simulation = list(localadapt, non_localadapt),
  ) |>
    rowwise() |>
    mutate(
      shifted_fitness = list(simulation[["Future fitness"]]),
      Y = list(simulation$Genotype),
      X = list(matrix(c(
        simulation[["Current env 1"]], simulation[["Current env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)),
      Xstar = list(matrix(c(
        simulation[["Future env 1"]], simulation[["Future env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)),
    )

  incomplete_data  <- tibble(
    scenario = c("Local adaptation", "Non local adaptation"),
    inference = "incomplete",
    simulation = list(localadapt, non_localadapt),
  ) |>
    rowwise() |>
    mutate(
      shifted_fitness = list(simulation[["Future fitness"]]),
      loci = list(c(
        sample(simulation[["Index QTLs 1"]], 5), sample(simulation[["Index QTLs 2"]], 5),
        seq(ncol(simulation$Genotype)-100, ncol(simulation$Genotype))
        )),
      Y = list(simulation$Genotype[, loci]),
      X = list(matrix(c(
        simulation[["Current env 1"]], simulation[["Current env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)),
      Xstar = list(matrix(c(
        simulation[["Future env 1"]], simulation[["Future env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)),
    )

  only_qtl1_data  <- tibble(
    scenario = c("Local adaptation", "Non local adaptation"),
    inference = "only_QTL1",
    simulation = list(localadapt, non_localadapt),
  ) |>
    rowwise() |>
    mutate(
      shifted_fitness = list(simulation[["Future fitness"]]),
      loci = list(c(
        sample(simulation[["Index QTLs 1"]], 5),
        seq(ncol(simulation$Genotype)-105, ncol(simulation$Genotype))
        )),
      Y = list(simulation$Genotype[, loci]),
      X = list(matrix(c(
        simulation[["Current env 1"]], simulation[["Current env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)),
      Xstar = list(matrix(c(
        simulation[["Future env 1"]], simulation[["Future env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)),
    )

  only_qtl2_data  <- tibble(
    scenario = c("Local adaptation", "Non local adaptation"),
    inference = "only_QTL2",
    simulation = list(localadapt, non_localadapt),
  ) |>
    rowwise() |>
    mutate(
      shifted_fitness = list(simulation[["Future fitness"]]),
      loci = list(c(
        sample(simulation[["Index QTLs 2"]], 5),
        seq(ncol(simulation$Genotype)-105, ncol(simulation$Genotype))
        )),
      Y = list(simulation$Genotype[, loci]),
      X = list(matrix(c(
        simulation[["Current env 1"]], simulation[["Current env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)),
      Xstar = list(matrix(c(
        simulation[["Future env 1"]], simulation[["Future env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)),
    )

  bind_rows(causal_data, empirical_data, incomplete_data, only_qtl1_data, only_qtl2_data)
}

# Custom function that computes the offset
# Same interface as the rest from Gain et al. 
# Number of latent factors are computed using tracy.widom test
go_genetic_gap <- function(Y, X, X.pred, snps.set){
  Y <- Y[, snps.set]
  compute_k <- function(Y, threshold = 0.05) {
    infile <- tempfile(fileext = ".geno")
    pca.res <- pca(write.geno(Y, output.file = infile))
    m1 <- tracy.widom(pca.res)
    dirname(infile) |>
      list.files(sub('\\.geno$', '', basename(infile)), full.names = TRUE) |>
      unlink(recursive = TRUE)
    project <- sub('\\.geno$', '', basename(infile))
    remove.pcaProject(paste0(project, ".pcaProject"))
    sum(m1$pvalues < threshold)
  }
  genetic.gap(input = Y, env = X, pred.env = X.pred, K=compute_k(Y))$offset  
}

# Create a nicely formatted table for a given dataset and a offset function
handle_simulation_pair <- function(data, offset_function){
    data |>
      mutate(offset = list(offset_function(Y, X, Xstar, 1:ncol(Y)))) |>
      select(scenario, shifted_fitness, inference, offset) |>
      pivot_wider(names_from = "inference", values_from = "offset", names_prefix = "offset_") |>
      unnest_longer(c(
        shifted_fitness, offset_causal, offset_empirical,
        offset_incomplete, offset_only_QTL1, offset_only_QTL2
        ))
} 

offset_fns <- list(
  "Geometric GO" = go_genetic_gap,
  "RONA" = go_rona,
  "RDA GO" = go_rda
)

run <- function() {
  outfile <- "results/local_adaptation_scenarios/m3_offsets.csv"
  non_local_files <- list.files("steps/slim/", "m3.1.+.Rds", full.names = TRUE)
  local_files <- list.files("steps/slim/", "m3.2.+.Rds", full.names = TRUE)
  res <- map2(non_local_files, local_files, \(non_local_files, local_files) {
    data_sim <- parse_simulation(non_local_files, local_files)
    map(offset_fns, \(f) handle_simulation_pair(data_sim, f)) |> bind_rows(.id = "Method")
  })
  names(res) <- sub('\\.Rds$', '', basename(non_local_files))
  res <- bind_rows(res, .id = "File") |>
    separate_wider_delim(File, "_", names = c(NA, "seed", NA, NA, "QTLs")) |>
    mutate(QTLs = as.numeric(QTLs))
  if (!dir.exists(dirname(outfile))){
  dir.create(dirname(outfile))
  }
  write_csv(res, outfile)
}
run()