#'  GGMprojpred
#' @description Estimate Gaussian graphical models with projection predictive selection
#' @param X n by p data matrix
#' @param n_cl number of clusters for parallel
#' @param type regularized (using the horeshoe prior distribution) or non-regularized (Bayesian bootstrap with OLS)
#' @param iter number of saved posterior samples
#'
#' @return
#'
#' @examples
#' GGMprojpred()
#' @export

GGMprojpred <- function(X, n_cl, type, iter){

  # intiial fitting
  fit <- int_proj(X, n_cl)

  # reproject
  mats <- re_project(fit$beta_mat, fit$fit_cv)

  # selected variables using the "or-rule"
  or_select <- ifelse(mats$pcor_or == 0, 0, 1)

  # acheive symmetry
  beta_mat <- beta_symmetric(fit$beta_mat, mats$mat_temp)  * or_select

  # compute the inverse covariance matrix
  inv_cov <- beta_to_inv(beta_mat = beta_mat, or_select, X)

  list(pcor_mat = mats$pcor_or, inv_cov = inv_cov$mat_inv,
       or_pd = inv_cov$det_check, selected_mat = or_select)
}
