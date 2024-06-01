library(tidyverse)
library(LEA)
library(gradientForest)
library(conflicted)
library(future)
library(furrr)

source("src/R/gain_offset.R")
source("src/R/go_offset.R")
# So GF doesn't complain
l2norm <- \(u, v) sqrt(sum((u-v)^2))

# Main function of the script
compute_genomic_offset <- function(file, fn) {
  simulation <- read_rds(file)
  shifted_fitness <- simulation[["Future fitness"]]
  causal_loci <- c(simulation[["Index QTLs 1"]], simulation[["Index QTLs 2"]])
  Y <- simulation$Genotype
  X <- matrix(c(
        simulation[["Current env 1"]], simulation[["Current env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)
  X.pred <- matrix(c(
        simulation[["Future env 1"]], simulation[["Future env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)

  k <- max(1, compute_k(Y))
  test <- lfmm2(Y, X, k) |> lfmm2.test(Y, X, genomic.control = TRUE, full = TRUE)
  qvalues <- qvalue::qvalue(test$pvalues)$qvalues
  candidates <- as.numeric(which(qvalues < 0.05))
  
  offset_causal <- fn(Y, X, X.pred, causal_loci)
  offset_empirical <- fn(Y, X, X.pred, candidates)
  res <- tibble(file, offset_causal, offset_empirical, shifted_fitness)
  colnames(X) <- c("CurrentCausalEnv1", "CurrentCausalEnv2", "CurrentNonCausalEnv1", "CurrentNonCausalEnv2")
  colnames(X.pred) <- c("FutureCausalEnv1", "FutureCausalEnv2", "FutureNonCausalEnv1", "FutureNonCausalEnv2")
  cbind(res, X, X.pred)
}

run <- function(infiles, outfile){
  print("Running RONA...")
  rona <- infiles |>
    future_map(.options = furrr_options(seed = TRUE),
    \(f) compute_genomic_offset(f, go_rona)
    ) |>
    bind_rows() |>
    mutate(method = "RONA")

  print("Running RDA...")
  rda <- infiles |>
    future_map(.options = furrr_options(seed = TRUE),
    \(f) compute_genomic_offset(f, go_rda)
    ) |>
    bind_rows() |>
    mutate(method = "GO RDA")

  print("Running GF...")
  gf <- infiles |>
    future_map(.progress = TRUE, .options = furrr_options(seed = TRUE), 
    \(f) compute_genomic_offset(
      f, \(Y, X, X_pred, causal_set) go_gf(Y, X, X_pred, causal_set)$go
      )) |>
    bind_rows() |>
    mutate(method = "GF GO")

  print("Running Geometric GO...")
  geometric <- infiles |>
    future_map(.options = furrr_options(seed = TRUE),
    \(f) compute_genomic_offset(f,  \(Y, X, X_pred, causal_set) go_genetic_gap(Y, X, X_pred, causal_set, X))
    ) |>
    bind_rows() |>
    mutate(method = "Geometric GO")

  print("Writing output...")
  bind_rows(rona, rda, gf, geometric) |>
    separate_wider_delim(file, "_", names = c("slim", "seed", NA, NA, "QTLs")) |>
    mutate(
      slim = basename(slim),
      QTLs = str_remove_all(QTLs, ".Rds") |> as.integer(),
      #m3.1 non local, m3.2 local, other NA
      scenario = case_when(
        grepl("m3.1", slim) ~ "Non local adaptation",
        grepl("m3.2", slim) ~ "Home-away & local-foreign",
        grepl("m3.5", slim) ~ "Local-foreign",
        TRUE ~ NA_character_
      )
    ) |>
    write_csv(outfile)
}

plan(cluster)
print(paste0(c("Available workers:", availableWorkers())))


args <- commandArgs(trailingOnly=TRUE)
print(args)
infiles <- args[1:length(args)-1]
outfile <- args[length(args)]
print(infiles)
run(infiles, outfile)