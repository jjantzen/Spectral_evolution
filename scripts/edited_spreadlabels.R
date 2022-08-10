spreadlabels <- function (tree, x, fsize = 1, cost = c(1, 1), range = NULL, 
          label.pos = NULL) 
{
  if (!is.null(label.pos)) 
    return(label.pos[tree$tip.label])
  else {
    if (is.null(range)) 
      range <- range(x)
    yy <- x[1:Ntip(tree)]
    zz <- setNames((rank(yy, ties.method = "random") - 1)/(length(yy) - 
                                                             1) * diff(range(yy)) + range(yy)[1], names(yy))
    mm <- max(fsize * strheight(tree$tip.label))
    ff <- function(zz, yy, cost, mo = 1, ms = 1) {
      ZZ <- cbind(zz - mm/2, zz + mm/2)
      ZZ <- ZZ[order(zz), ]
      oo <- 0
      for (i in 2:nrow(ZZ)) oo <- if (ZZ[i - 1, 2] > ZZ[i, 
                                                        1]) 
        oo <- oo + ZZ[i - 1, 2] - ZZ[i, 1]
      else oo <- oo
      pp <- sum((zz - yy)^2)
      oo <- if (oo < (1e-06 * diff(par()$usr[3:4]))) 
        0
      else oo
      pp <- if (pp < (1e-06 * diff(par()$usr[3:4]))) 
        0
      else pp
      oo/mo * cost[1] + pp/ms * cost[2]
    }
    mo <- ff(yy, zz, cost = c(1, 0))
    ms <- ff(yy, zz, cost = c(0, 1))
    if (mo == 0 && ms == 0) 
      return(yy)
    else {
      rr <- optim(zz, ff, yy = yy, mo = mo, ms = ms, cost = cost, 
                  method = "L-BFGS-B", lower = rep(range[1], length(yy)), 
                  upper = rep(range[2], length(yy)))
      return(rr$par)
    }
  }
}
