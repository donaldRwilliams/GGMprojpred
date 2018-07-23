# intial fit (re projected for OR rule)
hs_proj <- function(X, n_cl, iter){
  cl <- parallel::makeCluster(n_cl)
  doSNOW::registerDoSNOW(cl)
  colnames(X) <- 1:ncol(X)
  X <- scale(X, scale = F)
  par_fit <- foreach(i = 1:ncol(X), .export = c("node_fit", "project_node_fit",
                                                "bb_fit", "post_func"),
                     .packages = c("horseshoe", "projpred")) %dopar% {
                       X_new <- X[,-i]
                       y <- X[,i]

                       fit <- node_fit(X = X_new, y = y, iter = iter)

                       selected <- project_node_fit(X_new, y, beta = as.matrix(fit$beta), sigma = fit$sigma)

                     }
  parallel::stopCluster(cl)

  mat_select <- mat_all <- matrix(0, ncol = ncol(X), nrow = ncol(X))
  colnames(mat_select) <- 1:ncol(X)

  for(i in 1:ncol(X)){
    if(!is.null(par_fit[[i]]$mu_sel)){
      mat_select[i, ] <- rep(0,ncol(X))
    }
    if(!is.null(par_fit[[i]]$mu_sel)){
      mat_select[i, names(par_fit[[i]]$mu_sel)] <- par_fit[[i]]$mu_sel
    }

  }
  fit_cv <- list()
  for(i in 1:ncol(X)){
    fit_cv[[i]] <- par_fit[[i]]$fit_cv

  }
  fit_sigma <- list()
  for(i in 1:ncol(X)){
    fit_sigma[[i]] <- par_fit[[i]]$sigma
  }
  list(mat_sel = parcor::Beta2parcor(mat_select), beta_mat = mat_select, fit_cv = fit_cv, fit_sigma = fit_sigma)
}


bb_proj <- function(X, n_cl, iter){
  cl <- parallel::makeCluster(n_cl)
  doSNOW::registerDoSNOW(cl)
  colnames(X) <- 1:ncol(X)
  X <- scale(X, scale = F)
  par_fit <- foreach(i = 1:ncol(X), .export = c("node_fit", "project_node_fit",
                                                "bb_fit", "post_func"),
                     .packages = c("horseshoe", "projpred")) %dopar% {
                       X_new <- X[,-i]
                       y <- X[,i]




                       fit <- bb_fit(X = X_new, y = y, iter)

                       selected <- project_node_fit(X_new, y, beta = as.matrix(fit$beta), sigma = fit$sigma)

                     }
  parallel::stopCluster(cl)

  mat_select <- mat_all <- matrix(0, ncol = ncol(X), nrow = ncol(X))
  colnames(mat_select) <- 1:ncol(X)

  for(i in 1:ncol(X)){
    if(!is.null(par_fit[[i]]$mu_sel)){
      mat_select[i, ] <- rep(0,ncol(X))
    }
    if(!is.null(par_fit[[i]]$mu_sel)){
      mat_select[i, names(par_fit[[i]]$mu_sel)] <- par_fit[[i]]$mu_sel
    }

  }
  fit_cv <- list()
  for(i in 1:ncol(X)){
    fit_cv[[i]] <- par_fit[[i]]$fit_cv

  }
  fit_sigma <- list()
  for(i in 1:ncol(X)){
    fit_sigma[[i]] <- par_fit[[i]]$sigma
  }
  list(mat_sel = parcor::Beta2parcor(mat_select), beta_mat = mat_select, fit_cv = fit_cv, fit_sigma = fit_sigma)
}



node_fit <- function(X, y, burn = 1000, iter){
  fit <- horseshoe(y = y, X = X, method.tau = "halfCauchy", burn = burn,
                   method.sigma = "Jeffreys", nmc = iter)
  list(beta = t(fit$BetaSamples), sigma = sqrt(fit$Sigma2Samples))
}





project_node_fit <- function(X, y, beta, sigma){
  predfun <- function(xt) t(as.matrix(beta) %*% t(xt))
  reference <- init_refmodel(x = as.matrix(X), y = y, family = gaussian(), predfun = predfun,
                             intercept = F, dis = sigma)
  vs <- cv_varsel(reference, intercept = F, nloo = 25)
  nm <-  names(vs$varsel$pctch[vs$varsel$ssize,-1])
  if(vs$varsel$ssize == 0 |is.na(vs$varsel$ssize) ) {
    m <- rep(0, ncol(X))
  }
  if(is.numeric(vs$varsel$ssize)){
    p =  project(vs, vind =  vs$varsel$vind[1:vs$varsel$ssize])
    m <- colMeans(t(p$beta))
  }

  list(mu_sel = m, fit_cv = vs, sigma = sigma)
}



bb_fit <- function(X, y, iter){

  estimates <- do.call(rbind.data.frame, replicate(iter, post_func(X, y), simplify = F))
  list(beta = estimates[,1:ncol(X)], sigma = estimates[,ncol(estimates)])

}



post_func <- function(X, y){
  wts <- rexp(length(y), 1)
  wts <- wts / sum(wts)
  betas  <- solve(t(X) %*% diag(wts) %*% X) %*% t(X) %*% diag(wts) %*% y
  sig <- sd(y -  X %*% betas)
  estimates <-  t(as.matrix(c(t(as.matrix(betas)), sig)))
}


re_project <- function(x, fit_cv){
  beta_or <- list()
  select_mat <- ifelse(abs(x) == 0, 0, 1)
  mlt  <-  reshape2::melt(select_mat)
  sub_one <- subset(mlt, mlt[order(mlt$Var1), 3] != mlt[order(mlt$Var2),3])

  if(nrow(sub_one) == 0) {
    stop(list(pcor_and = parcor::Beta2parcor(x), pcor_and = parc))
  }
  for(i in 1:nrow(sub_one)){


    # i_th row corresponds to each regression model
    temp <- as.numeric(sub_one[i,])

    # temp[2] is the Beta weight
    temp_chr <- as.numeric(temp[2])

    # retrieve the regression model for re-projection
    temp_fit <- fit_cv[[temp[1]]]

    # check if condition is met, if not, i + 1
    # is.na denoting nothing was selected
    if(is.na(temp_fit$varsel$ssize)) {
      nm <- names(temp_fit$varsel$vind)[1]
    }#next(i)



    if(is.numeric(temp_fit$varsel$ssize)){
      nm <- names(temp_fit$varsel$vind)[1:temp_fit$varsel$ssize]
    }

    ####################################
    ### ensure only max is assessed ####
    ####################################
    if(is.na(match(temp_chr, colnames(temp_fit$varsel$pctch)[-1]))){
      selected <- nm
    }

    if(!is.na(match(temp_chr, colnames(temp_fit$varsel$pctch)[-1]))){
      selected <- unique(c(nm, temp_chr))
    }


    nv_max <- min(ncol(select_mat), floor(0.4*nrow(temp_fit$x)))
    if(nv_max > (ncol(select_mat) - 1 )) {nv_max <- ncol(select_mat) - 1}
    #if(nv_max > 20)
    # project
    beta_temp <- as.data.frame(t(project(temp_fit, nv = nv_max)$beta)[, as.character(selected)])

    # name projections
    colnames(beta_temp) <- selected

    # storage
    beta_or[[i]] <- as.matrix(beta_temp)

  }

  names(beta_or) <-  as.character(sub_one[,1])

  # temporary storage for the new projected matrix
  mat_temp <- matrix(0, nrow = ncol(x), ncol = ncol(x))
  colnames(mat_temp) <- 1:ncol(x)
  row.names(mat_temp)<- 1:ncol(x)
  for(i in 1:length(beta_or)){
    if(is.null(beta_or[[i]])) next(i)
    mat_temp[as.numeric(names(beta_or)[[i]]), colnames(beta_or[[i]])] <-
      as.numeric(colMeans(beta_or[[i]]))
  }

  # NAs as zeroes
  mat_temp[is.na(mat_temp)] <- 0

  # partial correlation matrix
  temp_or  <- parcor::Beta2parcor(mat_temp)
  pcor_and  <- parcor::Beta2parcor(x)


  or_upper <- as.numeric(temp_or[upper.tri(temp_or)])
  and_lower  <- as.numeric(pcor_and[upper.tri(pcor_and)])

  all <- ifelse(abs(or_upper) > 0, or_upper, and_lower)

  # matrix for "or" results
  mat_or <- matrix(nrow = ncol(x), ncol = ncol(x))
  mat_or[upper.tri(mat_or)] <- all
  mat_or <- as.matrix(forceSymmetric(mat_or))
  diag(mat_or) <- 1

  list(pcor_and = pcor_and, pcor_or = mat_or, mat_temp = mat_temp)
}

beta_symmetric <- function(initial, re_project){



  mat <- matrix(0, ncol(initial), ncol = ncol(initial))

  low <-  ifelse(abs(re_project[lower.tri(re_project)]) > 0 & initial[lower.tri(initial)] == 0,
                 re_project[lower.tri(re_project)], 0)
  up <-  ifelse(abs(re_project[upper.tri(re_project)]) > 0 & initial[upper.tri(initial)] == 0,
                re_project[upper.tri(re_project)], 0)
  mat[lower.tri(mat)] <-  low + initial[lower.tri(initial)]
  mat[upper.tri(mat)] <-  up + initial[upper.tri(initial)]

  mat
}


beta_to_inv <- function(beta_mat, or_select, X){
  # dimension
  beta_mat <- beta_mat * or_select
  p <- ncol(beta_mat)
  mat <- mat_inv <- matrix(0, p, p)
  # storage for residual sigma
  sig <- NA

  for(i in 1:p){
    # projected fitted values
    fitted <- beta_mat[i,-i] %*% t(X[,-i])
    # biased estimate of sigma
    sig[i] <- sd(fitted - X[,i])
    # beta to precision
    mat[i,-i] <- (beta_mat[i, -i] * - 1) / sig[i]^2
  }
  # take average
  mat_inv[upper.tri(mat_inv)] <-  (mat[upper.tri(mat)] + t(mat)[upper.tri(mat)]) / 2
  mat_inv[lower.tri(mat_inv)] <-  t(mat_inv)[lower.tri(mat_inv)]

  # inverse of variances
  diag(mat_inv) <- 1 / sig^2

  # structure of selected graph
  #mat_inv <- mat_inv * or_select

  # check if posiitive definite
  det_check <- det(mat_inv) > 0

  if(det_check == FALSE){
    for(i in 1:10000){
      print(i)
      d <- det(mat_inv)
      mat_inv  <- as.matrix(Matrix::nearPD(mat_inv, keepDiag = T, maxit = 10000)$mat) * or_select
      d <- det(mat_inv)
      if(d > 0) break
      if(i == 10000) mat_inv <- mat_inv + (abs(min(eigen(mat_inv)$values)) + 0.01)*diag(p)
    }}

  list(mat_inv = mat_inv, det_check = det_check)
}


