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
  Y <- Y[, snps.set]
  genetic.gap(input = Y, env = X, pred.env = X.pred, K=compute_k(Y))$offset  
}

go_genetic_gap_test <- function(Y, X, X.pred){
  test <- lfmm2(Y, X, compute_k(Y)) |> lfmm2.test(Y, X, genomic.control = TRUE, full = TRUE)
  snps.set <- which(test$pvalues < 0.05 / ncol(Y))
  if (length(snps.set) < 1) {
    return(rep(0, nrow(Y)))
    
  }
  go_genetic_gap(Y, X, X.pred, snps.set)
}