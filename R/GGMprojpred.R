#'  GGMprojpred
#' @description Estimate Gaussian graphical models with projection predictive selection
#' @param X n by p data matrix
#' @param n_cl number of clusters for parallel
#' @param type regularized (using the horeshoe prior distribution) or non-regularized (Bayesian bootstrap with least squares)
#' @param iter number of saved posterior samples
#'
#' @return pcor_mat estimated partial correlation matrix
#' @return inv_cov  estimated inverse covariance matrix
#' @return adj_mat  adjacenty matrix
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

GGMprojpred <- function(X, n_cl, iter){

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

  list(pcor_mat = mats$pcor_or, inv_cov = inv_cov$mat_inv, adj_mat = or_select)
}
