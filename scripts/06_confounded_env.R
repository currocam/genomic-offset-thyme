library(tidyverse)
library(LEA)
library(conflicted)
source("src/R/go_offset.R")

go_genetic_gap_test_possibly <- possibly(go_genetic_gap_test, NULL)

handle_simulation <- function(file){
  data <- read_rds(file)
  Y <- data$Genotype
  # Prepare env data
  X <- matrix(c(
    data[["Current env 1"]], data[["Current env 2"]]
  ),ncol=2) |>
    cbind(replicate(100, rnorm(100)))
  X.pred <-  matrix(c(
    data[["Future env 1"]], data[["Future env 2"]]
  ),ncol = 2) |>
    cbind(replicate(100, rnorm(100)))
  # Causal offset
  causal_loci <- c(data[["Index QTLs 1"]])
  causal <- go_genetic_gap(Y, X[,1], X.pred[,1], causal_loci)
  causal_one_confounded <- go_genetic_gap_test_possibly(Y, X[,c(1, 2)], X.pred[,c(1, 2)])
  # All SNPs
  empirical <- go_genetic_gap_test_possibly(Y, X[,1], X.pred[,1])
  empirical_one_confounded <- go_genetic_gap_test_possibly(Y, X[,c(1, 2)], X.pred[,c(1, 2)])
  res <- tibble(
    current_fitness = data[["Current fitness"]],
    future_fitness = data[["Future fitness"]],
    causal_offset = causal,
    causal_confounded_offset = causal_one_confounded,
    empirical_offset = empirical,
    empirical_confounded_offset = empirical_one_confounded,
  )
  walk(seq(1, 70, by = 5), \(index) {
    res[paste0("causal_random_", index, "_offset")] <<- go_genetic_gap(Y, X[,c(1, 3:(2+index))], X.pred[,c(1, 3:(2+index))], causal_loci)
    res[paste0("causal_random_", index, "_offset")] <<- go_genetic_gap_test_possibly(Y, X[,c(1, 3:(2+index))], X.pred[,c(1, 3:(2+index))])
  })
  res
}

handle_simulation2 <- function(file){
  data <- read_rds(file)
  mobility <- data$Mobility
  Y <- data$Genotype
  # Prepare env data
  X <- matrix(c(
    data[["Current env 1"]], data[["Current env 2"]]
  ),ncol=2) |>
    cbind(replicate(100, rnorm(100)))
  X.pred <-  matrix(c(
    data[["Future env 1"]], data[["Future env 2"]]
  ),ncol = 2) |>
    cbind(replicate(100, rnorm(100)))
  # Causal offset
  causal_loci <- c(data[["Index QTLs 1"]])
  causal <- go_genetic_gap(Y, X[,1], X.pred[,1], causal_loci)
  causal_one_confounded <- go_genetic_gap_test_possibly(Y, X[,c(1, 2)], X.pred[,c(1, 2)])
  # All SNPs
  empirical <- go_genetic_gap_test_possibly(Y, X[,c(1, 3, 4)], X.pred[,c(1, 3, 4)])
  empirical_one_confounded <- go_genetic_gap_test_possibly(Y, X[,1:4], X.pred[,1:4])
  tibble(
    mobility = mobility,
    current_fitness = data[["Current fitness"]],
    future_fitness = data[["Future fitness"]],
    causal_offset = causal,
    causal_confounded_offset = causal_one_confounded,
    empirical_offset = empirical,
    empirical_confounded_offset = empirical_one_confounded,
  )
}

run <- function() {
  set.seed(1000)
  outfile1 <- "results/local_adaptation_scenarios/m4_offsets_uncorrelated.csv"
  infiles <- c(list.files("steps/slim/", "m4.1.+.Rds", full.names = TRUE))
  names(infiles) <- basename(infiles) |> str_remove(".Rds")
  res <- map(infiles, handle_simulation) |>
    bind_rows(.id = "File") |>
    separate_wider_delim(File, "_", names = c("model", "seed","QTLs")) |>
    mutate(QTLs = str_remove(QTLs, "nQTL1s") |> as.numeric())
    
  if (!dir.exists(dirname(outfile1))){
  dir.create(dirname(outfile1))
  }
  write_csv(res, outfile1)

  outfile2 <- "results/local_adaptation_scenarios/m4_offsets_confounded.csv"
  infiles <- c(list.files("steps/slim/", "m4.2.+.Rds", full.names = TRUE))
  names(infiles) <- basename(infiles) |> str_remove(".Rds")
  res <- map(infiles, handle_simulation2) |>
    bind_rows(.id = "File") |>
    separate_wider_delim(File, "_", names = c("model", "seed","QTLs")) |>
    mutate(QTLs = str_remove(QTLs, "nQTL1s") |> as.numeric())
    
  if (!dir.exists(dirname(outfile2))){
  dir.create(dirname(outfile2))
  }
  write_csv(res, outfile2)
}

run()