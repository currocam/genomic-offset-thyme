library(tidyverse)
library(LEA)
library(conflicted)
library(future)
library(future.apply)
source("src/R/go_offset.R")

sample_individuals <- function(data, n){
    if (n > 0.5 * nrow(data$Genotype)) {
        return(sample(1:nrow(data$Genotype), n))
    }
    breaks <- c(0, 4, 8, 12)
    tibble(
      x = data[[paste0("LocationX_T", data$SamplingTimepoint)]],
      y = data[[paste0("LocationY_T", data$SamplingTimepoint)]],
      ) |>
      mutate(x = cut(x, breaks = breaks), y = cut(y, breaks = breaks)) |>
      mutate(row = row_number()) |>
      group_by(x, y) |>
      sample_n(ceiling(n/9), replace = TRUE) |>
      ungroup() |>
      sample_n(n) |>
      pull("row")
}

find_centers_geo <- function(data, K){
  loc <- tibble(x = data$LocationX_T0, y = data$LocationY_T0)
  kmeans(loc, centers = K)$centers
}

assign_population <- function(loc, centers){
  inner <- Vectorize(\(x, y){
    dists <- sqrt((centers[,1] - x)^2 + (centers[,2] - y)^2)
    which.min(dists)
  })
  inner(loc$x, loc$y) |> as.numeric()
}

handle_row <- function(file, sample_size, seed){

  message(
    paste0(
      "Processing file ", file, " with sample size ", sample_size, " and seed ", seed
    )
  )

  set.seed(seed)
  data <- read_rds(file)
  sampling_timepoint <- data$SamplingTimepoint
  last_timepoint <- data$LastTimepoint 
  Y <- data$Genotype
  population_size <- nrow(Y)
  colnames(Y) <- paste0("Locus", 1:ncol(Y))
  stopifnot(mean(data$Fitness_T0) > 0.9)

  n <- min(sample_size, nrow(Y))

  sample_index <- sample_individuals(data, n)
  Y <- Y[sample_index,]
  causal_loci <- paste0("Locus", data$`Index QTLs 1`)

  zero_variance_loci <- which(apply(Y, 2, var) == 0)
  if (length(zero_variance_loci) > 0) {
      Y <- Y[, -zero_variance_loci]
  }
  causal_loci_index <- which(colnames(Y) %in% causal_loci)

  X <- matrix(
    c(
      data[[paste0("Env1_T", sampling_timepoint)]],
      data[[paste0("Env2_T", sampling_timepoint)]],
      rnorm(population_size*2)
      ), ncol = 4
  )

  X.pred <- matrix(
    c(
      data[[paste0("PredEnv1_T", sampling_timepoint)]],
      data[[paste0("PredEnv2_T", sampling_timepoint)]],
      rnorm(population_size*2)
      ), ncol = 4
  )

  # Individual predicted 
  causal_ind_pred <- go_genetic_gap(Y, X[sample_index,1:2], X.pred[,1:2], causal_loci_index, X[,1:2])
  empirical_ind_pred <- go_genetic_gap_test(Y, X[sample_index,], X.pred, X)
  neglog_ind_pred <- -log(data[[paste0("PredFitness_T", sampling_timepoint)]])
  # Population predicted
  # First, we have to identify populations
  K <- 6*6
  pop_centers <- find_centers_geo(data, K)
  populations <- tibble(
    x = data[[paste0("LocationX_T", sampling_timepoint)]],
    y = data[[paste0("LocationY_T", sampling_timepoint)]]
  ) |>
  assign_population(pop_centers)

  causal_meanpop_pred <- 1:K |>
    map(\(pop) mean(causal_ind_pred[which(populations == pop)]) ) |>
    as.numeric()
  empirical_meanpop_pred <- 1:K |>
    map(\(pop) mean(empirical_ind_pred[which(populations == pop)]) ) |>
    as.numeric()
  neglog_meanpop_pred <- 1:K |>
    map(\(pop) mean(neglog_ind_pred[which(populations == pop)]) ) |>
    as.numeric() |>
    replace_na(Inf)
  # Actual population fitness
  future_populations <- tibble(
    x = data[[paste0("LocationX_T", last_timepoint)]],
    y = data[[paste0("LocationY_T", last_timepoint)]]
  ) |>
  assign_population(pop_centers)
  neglog_ind_actual <- -log(data[[paste0("PredFitness_T", last_timepoint)]])
  neglog_meanpop_actual <- 1:K |>
    map(\(pop) mean(neglog_ind_actual[which(future_populations == pop)]) ) |>
    as.numeric() |>
    replace_na(Inf)

  # Output
  tibble(
    file = file,
    sample_size = n,
    sampling_seed = seed,
    causal_ind_pred_cor = cor(causal_ind_pred, neglog_ind_pred),
    empirical_ind_pred_cor = cor(empirical_ind_pred, neglog_ind_pred),
    causal_meanpop_pred_cor = cor(causal_meanpop_pred, neglog_meanpop_pred),
    empirical_meanpop_pred_cor = cor(empirical_meanpop_pred, neglog_meanpop_pred),
    causal_meanpop_actual_cor = cor(causal_meanpop_pred, neglog_meanpop_actual),
    empirical_meanpop_actual_cor = cor(empirical_meanpop_pred, neglog_meanpop_actual),
    climate_change_rate = data$ClimateChangeRate
  )
}

handle_row_and_ignore_errors <- possibly(handle_row, NULL)

handle_file <- function(file, sample_sizes, seeds) {
  print("here")
  expand_grid(sample_size = sample_sizes, seed = seeds) |>
    pmap_dfr(\(sample_size, seed) handle_row_and_ignore_errors(file, sample_size, seed))
}

run <- function(infiles, outfile) {
    sample_size <- c(100, 300, 500)
    seeds <- c(1000, 2000, 3000, 4000)
    future_lapply(infiles, handle_file, sample_size, seeds, future.seed=TRUE, future.conditions = "message") |> 
        bind_rows() |>
        write_csv(outfile)
}

plan(cluster)
print(paste0(c("Available workers:", availableWorkers())))

args <- commandArgs(trailingOnly=TRUE)
infiles <- args[1:length(args)-1]
outfile <- args[length(args)]
run(infiles, outfile)