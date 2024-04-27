compute_k <- function(Y, threshold = 0.05) {
  infile <- tempfile(fileext = ".geno")
  pca.res <- pca(write.geno(Y, output.file = infile))
  m1 <- tracy.widom(pca.res)
  dirname(infile) |>
    list.files(sub('\\.geno$', '', basename(infile)), full.names = TRUE) |>
    unlink(recursive = TRUE)
  project <- sub('\\.geno$', '', basename(infile))
  remove.pcaProject(paste0(project, ".pcaProject"))
  sum(m1$pvalues < threshold)
}

go_genetic_gap <- function(Y, X, X.pred, snps.set){
  m.x <- apply(X, 2, mean)
  sd.x <- apply(X, 2, sd)
  X <- t(t(X) - m.x) %*% diag(1/sd.x)
  X.pred <- t(t(X.pred) - m.x) %*% diag(1/sd.x)
  Y <- Y[, snps.set]
  k <- max(1, compute_k(Y))
  genetic.gap(input = Y, env = X, pred.env = X.pred, K=k)$offset  
}

go_genetic_gap_test <- function(Y, X, X.pred){
  m.x <- apply(X, 2, mean)
  sd.x <- apply(X, 2, sd)
  X <- t(t(X) - m.x) %*% diag(1/sd.x)
  X.pred <- t(t(X.pred) - m.x) %*% diag(1/sd.x)
  k <- max(1, compute_k(Y))
  print(k)
  test <- lfmm2(Y, X, k) |> lfmm2.test(Y, X, genomic.control = TRUE, full = TRUE)
  candidate.loci <- which(test$pvalues < 0.05 / ncol(Y))
  if (length(candidate.loci) < 1) {
    return(rep(0, nrow(Y)))
  }
  genetic.gap(input = Y, env = X, pred.env = X.pred, K=k, candidate.loci = candidate.loci)$offset
}