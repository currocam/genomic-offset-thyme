library(tidyverse)
library(LEA)
library(conflicted)
library(future)
library(furrr)

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
  causal <- go_genetic_gap(Y, X[,1,drop=FALSE], X.pred[,1,drop=FALSE], causal_loci, X[,1,drop=FALSE])
  causal_one_confounded <- go_genetic_gap_test_possibly(Y, X[,c(1, 2)], X.pred[,c(1, 2)], X[,c(1, 2)])
  # All SNPs
  empirical <- go_genetic_gap_test_possibly(Y, X[,1,drop=FALSE], X.pred[,1,drop=FALSE], X[,1,drop=FALSE])
  empirical_one_confounded <- go_genetic_gap_test_possibly(Y, X[,c(1, 2)], X.pred[,c(1, 2)], X[,c(1, 2)])
  res <- tibble(
    current_fitness = data[["Current fitness"]],
    future_fitness = data[["Future fitness"]],
    causal_offset = causal,
    causal_confounded_offset = causal_one_confounded,
    empirical_offset = empirical,
    empirical_confounded_offset = empirical_one_confounded,
  )
  walk(seq(1, 70, by = 5), \(index) {
    res[paste0("causal_random_", index, "_offset")] <<- go_genetic_gap(Y, X[,c(1, 3:(2+index))], X.pred[,c(1, 3:(2+index))], causal_loci, X[,c(1, 3:(2+index))])
    res[paste0("causal_random_", index, "_offset")] <<- go_genetic_gap_test_possibly(Y, X[,c(1, 3:(2+index))], X.pred[,c(1, 3:(2+index))], X[,c(1, 3:(2+index))])
  })
  res
}

run <- function(infiles, outfile1) {
  set.seed(1000)
  names(infiles) <- basename(infiles) |> str_remove(".Rds")
  res <- future_map(infiles, handle_simulation, .options = furrr_options(seed = TRUE)) |>
    bind_rows(.id = "File") |>
    separate_wider_delim(File, "_", names = c("model", "seed","QTLs")) |>
    mutate(QTLs = str_remove(QTLs, "nQTL1s") |> as.numeric())
    
  if (!dir.exists(dirname(outfile1))){
  dir.create(dirname(outfile1))
  }
  write_csv(res, outfile1)
}

plan(cluster)
print(paste0(c("Available workers:", availableWorkers())))
args <- commandArgs(trailingOnly=TRUE)
print(args)
infiles <- args[1:length(args)-1]
outfile <- args[length(args)]
print(infiles)
run(infiles, outfile)
