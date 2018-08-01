#'  GGMprojpred
#' @description Estimate Gaussian graphical models with projection predictive selection
#' @param X n by p data matrix
#' @param n_cl number of clusters for parallel
#' @param type regularized (using the horeshoe prior distribution)
#' @param iter number of saved posterior samples
#' @param threshold if TRUE a threhold is applied after fitting (see Details)
#' @return pcor_mat estimated partial correlation matrix
#' @return inv_cov  estimated inverse covariance matrix
#' @return adj_mat  adjacenty matrix
#'
#' @details Thresholds are often used after a graphical structure has been estiamted, where small values are set to zero.
#' From our perspective, this is an ad-hoc "fix" to inrease specificity in which the decision theroretic jusfication for
#' predictive covariance selection is compromised. Nonetheless, for those wishing to lower the false positive rate (and also decrease power),
#' the threshold argement can be set to TRUE. The current threhold is based the Fisher-z tranformed standard error sqrt(1/ (n - s - 3)), with s as
#' the number of variables controlled for (p - 1). Values less than 1.96 * SE are set to zero.
#'
#' @examples
#' library(BDgraph)
#'
#' # generate AR-2 graph
#' main <- bdgraph.sim(p = 100, n = 50, graph = "AR2")
#' X <- main$data
#'
#' fit <- GGMprojpred(X, n_cl = detectCores() - 1, type = "hs", iter = 1000)
#' plot_graph(fit$pcor_mat, layout=layout_in_circle, vertex.color = "white")
#' compare(fit$adj_mat, main$G)
#' @export
#' @references
#' Piironen, J., & Vehtari, A. (2017). Comparison of Bayesian predictive methods for model selection. Statistics and Computing, 27(3), 711-735.
#' \href{https://link.springer.com/article/10.1007/s11222-016-9649-y}{https://link.springer.com/article/10.1007/s11222-016-9649-y}
#'
#' Rubin, D. B. (1981). The bayesian bootstrap. The annals of statistics, 130-134. \href{https://www.jstor.org/stable/2240875}{https://www.jstor.org/stable/2240875}
#'
#' Williams, D. R., Piironen, J., Vehtari, A., & Rast, P. (2018). Bayesian Estimation of Gaussian Graphical Models with Projection Predictive Selection. arXiv preprint arXiv:1801.05725.
#' \href{https://arxiv.org/abs/1801.05725}{https://arxiv.org/abs/1801.05725}

GGMprojpred <- function(X, n_cl, iter, threshold = FALSE){
    n <- nrow(X)
    p <- ncol(X)
    s <- p - 1
    if(threshold == TRUE & p >= n){
      stop("Thresholds are currently only implemented for low-dimensional data (n > p)")
    }
  # intiial fitting
    fit <- hs_proj(X, n_cl, iter)



   # reproject
   mats <- re_project(fit$beta_mat, fit$fit_cv)

  # selected variables using the "or-rule"
  or_select <- ifelse(mats$pcor_or == 0, 0, 1)

  # acheive symmetry
  beta_mat <- beta_symmetric(fit$beta_mat, mats$mat_temp)  * or_select

  # compute the inverse covariance matrix
  inv_cov <- beta_to_inv(beta_mat = beta_mat, or_select, X)

  if(threshold == TRUE){
    threshold_value <- sqrt(1 / (n - s - 3)) * 1.96
    threshold_mat <-  ifelse(abs(mats$pcor_or) > threshold_value, 1, 0)
    results <- list(pcor_mat = mats$pcor_or, inv_cov = inv_cov$mat_inv, adj_mat = or_select, threshold_mat = threshold_mat)
    } else{
     results <- list(pcor_mat = mats$pcor_or, inv_cov = inv_cov$mat_inv, adj_mat = or_select)
    }

  return(results)

  }
