library(tidyverse)
library(LEA)
library(conflicted)
library(future)
library(future.apply)
source("src/R/go_offset.R")

sample_individuals <- function(data, n){
    if (n > 0.5 * nrow(data$Genotype)) {
        sample(1:nrow(data$Genotype), n)
    }
    with(data, {
        breaks <- c(0, 4, 8, 12)
        tibble(x = LocationX_T0, y = LocationY_T0, row = 1:length(x)) |>
        # Bin the locations
        mutate(x = cut(x, breaks = breaks), y = cut(y, breaks = breaks)) |>
        group_by(x, y) |>
        sample_n(ceiling(n/9), replace = TRUE) |>
        # Shuffle the rows
        ungroup() |>
        sample_n(n) |>
        pull(row)
    })
}
handle_row <- function(file, sample_size, seed){
    set.seed(seed)
    data <- read_rds(file)
    Y <- data$Genotype
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

    X <- with(data, {
        matrix(c(Env1_T0, Env2_T0, rnorm(length(Env2_T0)), rnorm(length(Env2_T0))), ncol=4)
    })

    X.pred <- with(data, {
        matrix(c(PredEnv1_T0, PredEnv2_T0, rnorm(length(PredEnv2_T0)), rnorm(length(PredEnv2_T0))), ncol=4)
    })

    neglogfitness <- -log(data$PredFitness_T0)
    causal <- go_genetic_gap(Y, X[sample_index,1:2], X.pred[,1:2], causal_loci_index, X[,1:2])
    
    empirical <- go_genetic_gap_test(Y, X[sample_index,], X.pred, X)
    
    tibble_row(
        file = file,
        causal_correlation = cor(neglogfitness, causal),
        empirical_correlation = cor(neglogfitness, empirical),
        seed = seed,
        sample_size = sample_size
    ) 
}

handle_row_and_ignore_errors <- possibly(handle_row, NULL)

handle_file <- function(file, sample_sizes, seeds) {
    expand_grid(sample_size = sample_sizes, seed = seeds) |>
    pmap_dfr(\(sample_size, seed) handle_row_and_ignore_errors(file, sample_size, seed))
}

run <- function(infiles, outfile) {
    sample_size <- c(100, 200, 300, 400, 500, Inf)
    seeds <- c(1000, 2000, 3000, 4000, 5000)

    future_lapply(infiles, handle_file, sample_size, seeds, future.seed=TRUE) |> 
        bind_rows() |>
        write_csv(outfile)
}

plan(cluster)
print(paste0(c("Available workers:", availableWorkers())))

args <- commandArgs(trailingOnly=TRUE)
infiles <- args[1:length(args)-1]
outfile <- args[length(args)]
run(infiles, outfile)