library(tidyverse)
library(LEA)
library(cluster) 
library(popkin)
library(conflicted)
source("src/R/go_offset.R")

go_genetic_gap_test_possibly <- possibly(go_genetic_gap_test, NULL)

identify_subpopulations <- function(df){
  # n_pops <- 2:200
  # argmin <- n_pops |>
  #   map(\(k){
  #      km <- kmeans(df, centers = k, nstart=25)
  #      ss <- silhouette(km$cluster, dist(df))
  #      mean(ss[, 3])
  #   }) |>
  #   as.numeric() |>
  #   which.min()
  # k <- n_pops[[argmin]]
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
    cbind(replicate(100, rnorm(n)))
  X.pred <-  matrix(c(
    data[["Future env 1"]], data[["Future env 2"]]
  ),ncol = 2) |>
    cbind(replicate(100, rnorm(n)))
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

run <- function() {
  set.seed(1000)
  outfile2 <- "results/local_adaptation_scenarios/m4_offsets_confounded.csv"
  infiles <- c(list.files("steps/slim/", "m4.2.+.Rds", full.names = TRUE))
  names(infiles) <- basename(infiles) |> str_remove(".Rds")

  res <- map(infiles, handle_simulation) |>
    bind_rows(.id = "File") |>
    separate_wider_delim(File, "_", names = c("model", "seed",NA, NA, "QTLs")) |>
    mutate(QTLs = as.numeric(QTLs))

  mob <- map(infiles, \(x) read_rds(x)[["Later mobility"]]) |> as.numeric()
  fst <- res |> group_by(seed) |> summarize (Fst = unique(fst)) |> pull("Fst")
    
  stopifnot(cor(fst, mob) < -0.90)

  if (!dir.exists(dirname(outfile2))){
  dir.create(dirname(outfile2))
  }
  write_csv(res, outfile2)
}

run()