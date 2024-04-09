library(tidyverse)
library(LEA)
library(conflicted)
source("src/R/go_offset.R")

handle_simulation <- function(file){
  data <- read_rds(file)
  Y <- data$Genotype
  # Add 100 random snps at the end
  freqs <- runif(100, 0, 1)
  Yrand <- replicate(nrow(Y), rbinom(n = 100, size = 2, prob = freqs))
  Y <- cbind(Y, Yrand)
  # Prepare env data
  X <- matrix(c(
    data[["Current env 1"]], data[["Current env 2"]],
    rnorm(100), rnorm(100)
  ),ncol=4)
  X.pred <-  matrix(c(
    data[["Future env 1"]], data[["Future env 2"]],
    rnorm(100), rnorm(100)
  ),ncol=4)
  # Causal offset
  causal_loci <- c(data[["Index QTLs 1"]], data[["Index QTLs 2"]])
  causal <- go_genetic_gap(Y, X[,1:2], X.pred[,1:2], causal_loci)
  # All SNPs
  empirical <- go_genetic_gap(Y, X, X.pred, 1:ncol(Y))
  filtered <- go_genetic_gap_test(Y, X, X.pred)
  # Incomplete QTLs
  incomplete_loci <- c(
    sample(data[["Index QTLs 1"]], 5), sample(data[["Index QTLs 2"]], 5),
    # random snps
    seq(ncol(Y)-94, ncol(Y))
  )
  incomplete <- go_genetic_gap(Y, X, X.pred, incomplete_loci)
  # Missing one qtl
  only_qtl1 <- c(
    sample(data[["Index QTLs 1"]], 5),
    # random snps
    seq(ncol(Y)-99, ncol(Y))
  )
  only_one <- go_genetic_gap(Y, X, X.pred, only_qtl1)
  only_qtl2 <- c(
    sample(data[["Index QTLs 2"]], 5),
    # random snps
    seq(ncol(Y)-99, ncol(Y))
  )
  only_two <- go_genetic_gap(Y, X, X.pred, only_qtl2)
  # Only random
  random_snps <- seq(ncol(Y)-99, ncol(Y))
  random <- go_genetic_gap(Y, X, X.pred, random_snps)
  fn <- possibly(go_genetic_gap_test, NULL)
  random_putative <- fn(Y[,random_snps], X, X.pred)
  env_dist <- sqrt(rowSums((X-X.pred)^2))
  tibble(
    current_fitness = data[["Current fitness"]],
    future_fitness = data[["Future fitness"]],
    causal_offset = causal,
    empirical_offset = empirical,
    empirical_putative = filtered,
    incomplete_offset = incomplete,
    only_QTL1_offset = only_one,
    only_QTL2_offset = only_two,
    random_offset = random,
    random_putative = random_putative,
    env_dist = env_dist
  )
}

run <- function() {
  outfile <- "results/local_adaptation_scenarios/m3_offsets_asymmetric.csv"
  infiles <- c(
    list.files("steps/slim/", "m3.3.+.Rds", full.names = TRUE),
    list.files("steps/slim/", "m3.4.+.Rds", full.names = TRUE)
  )
  names(infiles) <- basename(infiles) |> str_remove(".Rds")
  set.seed(1000)
  res <- map(infiles, handle_simulation) |>
    bind_rows(.id = "File") |>
    separate_wider_delim(File, "_", names = c("model", "seed", NA, NA, "QTLs")) |>
    mutate(QTLs = as.numeric(QTLs))
  
  if (!dir.exists(dirname(outfile))){
    dir.create(dirname(outfile))
  }
  write_csv(res, outfile)
}

run()