library(tidyverse)
library(LEA)
library(conflicted)
source("src/R/go_offset.R")
source("src/R/go_offset.R")

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
  causal_one_counfounded <- go_genetic_gap(Y, X[,c(1, 2)], X.pred[,c(1, 2)], causal_loci)
  # All SNPs
  empirical <- go_genetic_gap(Y, X[,1], X.pred[,1], 1:ncol(Y))
  empirical_one_counfounded <- go_genetic_gap(Y, X[,c(1, 2)], X.pred[,c(1, 2)], 1:ncol(Y))
  res <- tibble(
    current_fitness = data[["Current fitness"]],
    future_fitness = data[["Future fitness"]],
    causal_offset = causal,
    causal_counfounded_offset = causal_one_counfounded,
    empirical_offset = empirical,
    empirical_counfounded_offset = empirical_one_counfounded,
  )
  walk(seq(1, 70, by = 5), \(index) {
    res[paste0("causal_random_", index, "_offset")] <<- go_genetic_gap(Y, X[,c(1, 3:(2+index))], X.pred[,c(1, 3:(2+index))], causal_loci)
    res[paste0("causal_random_", index, "_offset")] <<- go_genetic_gap(Y, X[,c(1, 3:(2+index))], X.pred[,c(1, 3:(2+index))], 1:ncol(Y))
  })
  res
}

run <- function() {
  outfile <- "results/local_adaptation_scenarios/m4_offsets_counfounded.csv"
  infiles <- c(list.files("steps/slim/", "m4.1.+.Rds", full.names = TRUE))
  names(infiles) <- basename(infiles) |> str_remove(".Rds")
  res <- map(infiles, handle_simulation) |>
    bind_rows(.id = "File") |>
    separate_wider_delim(File, "_", names = c("model", "seed","QTLs")) |>
    mutate(QTLs = str_remove(QTLs, "nQTL1s") |> as.numeric())
    
  if (!dir.exists(dirname(outfile))){
  dir.create(dirname(outfile))
  }
  write_csv(res, outfile)
}

run()
