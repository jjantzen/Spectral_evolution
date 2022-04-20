library(ggrepel)
#plotting function - figure out ellipses
plot_pca_phylo <- function(pca_object, model, axes, groups,  cols, labels) {
  #getting values for plotting
  U <- pca_object$vectors
  resids <- model$residuals
  S <- resids %*% U
  #make as dataframe for plotting and name columns
  df_S <- as.data.frame(S[,axes[1:2]])
  names(df_S) <- c("xvar", "yvar")
  #assign groups based on separate trait for colouring
  #groups <- levels(factor(col))
  df_S$groups <- groups
  #get labels for axes
  tot <- sum(pca_object$values)
  valX <- round(pca_object$values[axes[1]] * 100/tot, digits = 2)
  valY <- round(pca_object$values[axes[2]] * 100/tot, digits = 2)
  xlabel <- paste("PC", axes[1], " (", valX, " %)", sep = "")
  ylabel <- paste("PC", axes[2], " (", valY, " %)", sep = "")
  #get parameters for ellipses
  #ellipse.prob <- 0.68
  #theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  #circle <- cbind(cos(theta), sin(theta))
  #ell <- ddply(df_S, "groups", function(x) {
  #  sigma <- var(cbind(df_S$xvar, df_S$yvar))
  #  mu <- c(mean(df_S$xvar), mean(df_S$yvar))
  #  ed <- sqrt(qchisq(ellipse.prob, df = 2))
  #  data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = "+"), groups = df_S$groups)
  #})
  #names(ell)[1:2] <- c("xvar", "yvar")
  #str(ell)
  #labels for points
  df_S$labels <- labels
  ggplot(data = df_S, aes(x = xvar, y = yvar))+
    geom_point(aes(colour = groups), shape = 19, size = 2)+
    geom_text_repel(label = labels, max.overlaps = 3, size = 3)+
    #geom_path(data = ell[which(ell$groups == 0),])#, aes(colour = groups, group = groups))+ #, inherit.aes = FALSE
    geom_hline(aes(yintercept = 0))+
    geom_vline(aes(xintercept = 0))+
    labs(y = ylabel, x = xlabel)+
    scale_colour_manual(name = "Mycorrhizal association", labels = c("AM", "EM"), values = rev(levels(factor(cols))))+
    theme_bw()+
    theme(axis.title=element_text(size=14))
}

