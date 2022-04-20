#' Pick PLS ncomp from `caret`
#'
#' @param x fit object returned from `caret::train`
#' @param sigma number of standard deviations
#' @param SE use Standard Errors? Defaults to FALSE
#'
#' @return the best number of components
#' @export
#'
#' @author Jose Eduardo Meireles
pick_pls_ncomp_caret = function(x, sigma = 1, SE = FALSE){

  idx_max = which.max(x$results[ , x$metric])
  mu_max  = x$results[ idx_max , x$metric]
  sd_max  = x$results[ idx_max , paste(x$metric, "SD", sep = "") ]

  if(SE){
    sd_max = sd_max / sqrt(nrow(x$resample))
  }

  low = mu_max - sd_max * sigma

  pick = min(which(x$results[ , x$metric] >= low))

  x$results[ pick, "ncomp"]
}
