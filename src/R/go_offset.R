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