dwald <- function(t, alpha, tau, kappa, give_log=FALSE){
  dwald_gamma(t, alpha, tau, kappa, give_log)
}

dIG <- function(t, alpha, lambda, nu, give_log=FALSE){
  d <- vector("double", length=length(t))
  for (i in 1:length(t)) {
      d[i] <- .Call(dIG_c, t[i], alpha, lambda, nu, give_log)
  }
  return(d)
}

dwald_gamma <- function(t, alpha, tau, kappa, give_log=FALSE){
  d <- vector("double", length=length(t))
  for (i in 1:length(t)) {
    d[i] <- .Call(dwald_gamma_c, t[i], alpha, tau, kappa, give_log)
  }
  return(d)
}

dwald_trunc <- function(t, lambda, alpha, v, delta, give_log=FALSE){
  d <- vector("double", length=length(t))
  for (i in 1:length(t)) {
    d[i] <- .Call(dwald_trunc_c, t[i], lambda, alpha, v, delta, give_log)
  }
  return(d)
} 

dswald <- function(t, alpha, nu, theta, give_log=FALSE){
  d <- vector("double", length=length(t))
  for (i in 1:length(t)) {
    d[i] <- .Call(dswald_c, t[i], alpha, nu, theta, give_log)
  }
  return(d)
}

pswald <- function(t, alpha, nu, theta, lower.tail=TRUE, log.p=FALSE){
  d <- vector("double", length=length(t))
  for (i in 1:length(t)) {
    d[i] <- .Call(pswald_c, t[i], alpha, nu, theta, lower.tail, log.p)
  }
  return(d)
}

laguerrel <- function(n, a, x) {
  d <- .Call(laguerrel_c, n, a, x)
  return(d)
}
