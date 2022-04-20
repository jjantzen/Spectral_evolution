mvgls.pca_ed <- function (object, ...) 
{
  args <- list(...)
  if (is.null(args[["axes"]])) 
    axes <- c(1, 2)
  else axes <- args$axes
  if (is.null(args[["las"]])) 
    las <- 1
  else las <- args$las
  if (is.null(args[["main"]])) {
    if (object$method == "LL") 
      main <- "Phylogenetic PCA"
    else main <- "Regularized Phylogenetic PCA"
  }
  else {
    main <- args$main
  }
  if (is.null(args[["mode"]])) 
    mode <- "cov"
  else mode <- args$mode
  if (!inherits(object, "mvgls")) 
    stop("only works with \"mvgls\" class objects. See ?mvgls")
  covR <- object$sigma$Pinv
  if (mode == "corr") 
    covR <- cov2cor(covR)
  eig <- eigen(covR)
  values <- eig$values
  U <- eig$vectors
  resids <- object$residuals
  S <- resids %*% U
  res <- list(scores = S, values = values, vectors = U, rank = qr(covR)$rank)
  class(res) <- "mvgls.pca"
  invisible(res)
}
