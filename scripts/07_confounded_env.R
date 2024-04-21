#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(LEA)
library(cluster) 
library(popkin)
library(conflicted)
source("src/R/go_offset.R")

go_genetic_gap_test_possibly <- possibly(go_genetic_gap_test, NULL)

identify_subpopulations <- function(df){
  kmeans(df, centers = 11*11, nstart=100)$cluster
}

handle_simulation <- function(file){
  data <- read_rds(file)
  
  subpops <- c(data[["Location X"]], data[["Location Y"]]) |>
    matrix(ncol = 2) |>
    identify_subpopulations()
  
  kinship <- data[["Genotype"]] |>
    t() |>
    popkin(subpops)
  w <- weights_subpops(subpops)
  fst <- fst(kinship, w)
  
  Y <- data[["Genotype"]]
  variable_loci <- which(apply(Y, 2, sd) > 0)
  

  # Prepare env data
  n <- nrow(Y)
  X <- matrix(c(
    data[["Current env 1"]], data[["Current env 2"]]
  ),ncol=2) |>
    cbind(replicate(2, rnorm(n)))
  X.pred <-  matrix(c(
    data[["Future env 1"]], data[["Future env 2"]]
  ),ncol = 2) |>
    cbind(replicate(2, rnorm(n)))
  # Causal offset
  causal_loci <- c(data[["Index QTLs 1"]])
  causal <- go_genetic_gap(Y, X[,1], X.pred[,1], causal_loci)
  causal_one_confounded <- go_genetic_gap(Y, X[,c(1, 2)], X.pred[,c(1, 2)], causal_loci)
  # All SNPs
  empirical <- go_genetic_gap_test_possibly(Y[,variable_loci], X[,c(1, 3, 4)], X.pred[,c(1, 3, 4)])
  empirical_one_confounded <- go_genetic_gap_test_possibly(Y[,variable_loci], X[,1:4], X.pred[,1:4])
  tibble(
    fst = fst,
    mobility = data[["Later mobility"]],
    current_fitness = data[["Current fitness"]],
    future_fitness = data[["Future fitness"]],
    causal_offset = causal,
    causal_confounded_offset = causal_one_confounded,
    empirical_offset = empirical,
    empirical_confounded_offset = empirical_one_confounded,
  )
}

run <- function(args) {
  set.seed(1000)
  outfile <- tail(args, 1)
  infiles <- head(args, length(args)-1)
  names(infiles) <- basename(infiles) |> str_remove(".Rds")

  res <- map(infiles, handle_simulation) |>
    bind_rows(.id = "File") |>
    separate_wider_delim(File, "_", names = c("model", "seed",NA, NA, "QTLs")) |>
    mutate(QTLs = as.numeric(QTLs))

  mob <- res |> group_by(model, seed, QTLs) |> summarize (mobility = unique(mobility)) |> pull("mobility")
  fst <- res |> group_by(model, seed, QTLs) |> summarize (Fst = unique(fst)) |> pull("Fst")

  print("Correlation: ")  
  print(cor(fst, mob))
  stopifnot(cor(fst, mob) < -0.5)


  if (!dir.exists(dirname(outfile))){
  dir.create(dirname(outfile))
  }
  write_csv(res, outfile)
}

run(args)