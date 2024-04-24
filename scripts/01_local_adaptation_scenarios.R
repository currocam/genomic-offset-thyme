library(tidyverse)
library(LEA)
library(gradientForest)
library(conflicted)
source("src/R/gain_offset.R")
source("src/R/go_offset.R")
# So GF doesn't complain
l2norm <- \(u, v) sqrt(sum((u-v)^2))

# Main function of the script
compute_genomic_offset <- function(file, causal_function, empirical_function) {
  simulation <- read_rds(file)
  shifted_fitness <- simulation[["Future fitness"]]
  causal_loci <- c(simulation[["Index QTLs 1"]], simulation[["Index QTLs 2"]])
  Y <- simulation$Genotype
  X <- matrix(c(
        simulation[["Current env 1"]], simulation[["Current env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)
  X.pred <- matrix(c(
        simulation[["Current env 1"]], simulation[["Current env 2"]],
        rnorm(100), rnorm(100)
      ),ncol=4)
  
  offset_causal <- causal_function(Y, X, X.pred, causal_loci)
  offset_empirical <- empirical_function(Y, X, X.pred)
  tibble(file, offset_causal, offset_empirical, shifted_fitness)
}

run <- function(infiles, outfile){
  print("Running RONA...")
  rona <- infiles |>
    map(\(f) compute_genomic_offset(
      f, go_rona, \(Y, X, X.pred) go_rona(Y, X, X.pred, 1:ncol(Y)))
      ) |>
    bind_rows() |>
    mutate(method = "rona")

  print("Running RDA...")
  rda <- infiles |>
    map(\(f) compute_genomic_offset(
      f, go_rda, \(Y, X, X.pred) go_rda(Y, X, X.pred, 1:ncol(Y)))
      ) |>
    bind_rows() |>
    mutate(method = "GO RDA")

  print("Running GF...")
  gf <- infiles |>
    map(\(f) compute_genomic_offset(
      f, \(Y, X, X_pred, causal_set) go_gf(Y, X, X_pred, causal_set)$go, \(Y, X, X_pred, causal_set) go_gf(Y, X, X_pred, 1:ncol(Y))$go )
      ) |>
    bind_rows() |>
    mutate(method = "GF GO")

  print("Running Geometric GO...")
  geometric <- infiles |>
    map(\(f) compute_genomic_offset(f, go_genetic_gap, go_genetic_gap_test)) |>
    bind_rows() |>
    mutate(method = "Geometric GO")

  print("Writing output...")
  bind_rows(rona, rda, gf, geometric) |>
    separate_wider_delim(file, "_", names = c("slim", "seed", NA, NA, "QTLs")) |>
    mutate(
      slim = basename(slim),
      #m3.1 non local, m3.2 local, other NA
      scenario = case_when(
        grepl("m3.1", slim) ~ "Non local adaptation",
        grepl("m3.2", slim) ~ "Local adaptation",
        TRUE ~ NA_character_
      )
    ) |>
    write_csv(outfile)
}

args <- commandArgs(trailingOnly=TRUE)
print(args)
infiles <- args[1:length(args)-1]
outfile <- args[length(args)]
print(infiles)
run(infiles, outfile)