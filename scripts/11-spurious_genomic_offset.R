library(JuliaConnectoR)
library(tidyverse)
stopifnot(juliaSetupOk())
# Activate julia env for reproducibility
juliaEval('using Pkg; Pkg.activate(".")')
# Finally, load the required Julia package
GO <- juliaImport("GenomicOffsets")

adjusted_r2 <- function(x, y){
  model <- lm(y ~ x)
  summary(model)$adj.r.squared
}

compute_row <- function(seed, n_uncorrelated, L){
  set.seed(seed)
  x <- matrix(runif(n = 100)) # Causal env variable
  optimum_phenotype <- 1*x
  phenotypes <- optimum_phenotype + rnorm(n = 100, 0, 0.1) # nearly-optimal
  current_fitness <- exp(-(optimum_phenotype - phenotypes)^2 / 2)
  x.modified <- matrix(x + runif(n = 100, -1, 1))
  modified_optimum_phenotype <- 1*x.modified
  modified_fitness <- exp(-(modified_optimum_phenotype - phenotypes)^2 / 2)
  # No genetic compoent
  allele_frequencies <- runif(L)
  Y <- sapply(allele_frequencies, \(f) rbinom(100, 2, f))
  print(dim(Y))
  
  X <- x
  Xpred <- x.modified
  if (n_uncorrelated>0){
    X <- cbind(X, replicate(n= n_uncorrelated, rnorm(100)))
    Xpred <- cbind(Xpred, replicate(n= n_uncorrelated, rnorm(100)))
  }
  neglog <- -log(modified_fitness)
  model <- GO$fit(GO$GeometricGO, Y, X, tw_threshold = 1e-3)
  # You can check the chosen K
  print(model$K)
  pvalues <- GO$LFMM_Ftest(model$model, Y, X, genomic_control = TRUE, center = TRUE)
  qvalues <- qvalue::qvalue(pvalues)$qvalues
  candidates <- as.integer(which(qvalues < 0.05))
  with_test <- rep(0, length(neglog))
  if (length(candidates)> 0){
    with_test <- as.numeric( GO$genomic_offset(model, X, Xpred, as.list(candidates)))
  }
  # Select 10% of the loci at random
  Y <- Y[, sample(1:L, L/10)]
  model <- GO$fit(GO$GeometricGO, Y, X, tw_threshold = 1e-3)
  print(length(candidates))
  tibble(
    with_testing = adjusted_r2(neglog, with_test),
    without_testing = adjusted_r2(neglog, GO$genomic_offset(model, X, Xpred)),
    seed = seed, 
    L = L, 
    n_uncorrelated = n_uncorrelated
  )
}
compute_row_safe <- possibly(compute_row, otherwise = tibble(with_testing = NA_real_, without_testing = NA_real_, seed = NA_integer_, L = NA_integer_, n_uncorrelated = NA_integer_))

args <- commandArgs(trailingOnly=TRUE)
outfile <- args[[1]]
expand_grid(seed = seq(1, 15), n_uncorrelated = c(0, 1, 2, 5, 10, 20), L = c(50, 500, 1000)) |>
  pmap(compute_row_safe, .progress=TRUE) |>
  bind_rows() |>
  write_csv(outfile)