#functions.R


#'
#' @export
#'
#' @examples
#' mu <- rnorm(10)
#' Sigma <- rWishart(1,11,diag(1,10))[,,1]
#' condparam_W(Sigma)
#'
condparam_W <- function(Sigma){
  library(magrittr)
  param <- 2:dim(Sigma)[1] %>%
    purrr::map(function(i){
      coef <- solve(Sigma[1:(i-1),1:(i-1)],Sigma[1:(i-1),i])
      var <-  Sigma[i,i]-sum(coef * Sigma[1:(i-1),i])
      list(var=var,coef=coef)
    }
    )
  c(list(list(var=Sigma[1,1],coef=0)),param)
}


#'
#' @export
#'
#' @examples
#' mu <- rnorm(10)
#' Sigma <- rWishart(1,11,diag(1,10))[,,1]
#' yval <- sample(0:1,10,replace=T)
#' lbvec  <- -1/yval +1
#' ubvec <- (1/!yval)-1
#' condparam <- condparam_W(Sigma)
#' gen_W(mu,condparam,lbvec,ubvec)
#'
#'
gen_W <- function(mu,condparam,lbvec,ubvec){
  W <- numeric(length(lbvec))
  W[1] <- extraDistr::rtnorm(1,mean=mu[1],sd=sqrt(condparam[[1]]$var),
                             a=lbvec[1],b=ubvec[1])
  for(i in 2:length(lbvec)){
    W[i] <- extraDistr::rtnorm(1,mean=mu[i] +sum(condparam[[i]]$coef *( W[1:(i-1)]-mu[1:(i-1)])),
                               sd=sqrt(condparam[[i]]$var),a=lbvec[i],b=ubvec[i])

  }
  W
#  mytruncNormal(mu,Sigma,lbvec,ubvec)
}

#'
#' @export
#'
gen_Wmat <- function(mumat,condparam,lbmat,ubmat){
  1:dim(lbmat)[1] %>%
    purrr::map(~gen_W(mu=mumat[.x,],condparam=condparam,
                      lbvec=lbmat[.x,],ubvec=ubmat[.x,]))%>%
    do.call("rbind",.)
}


#'
#' @export
#'
#' @examples
#' mu <- rnorm(10)
#' Sigma <- diag(runif(10)+1,10)
#' yval <- sample(0:1,10,replace=T)
#' lbvec  <- -1/yval +1
#' ubvec <- (1/!yval)-1
#'
#'
mytruncNormal <- function(mu,Sigma,lbvec,ubvec,Winit=NULL){
  1:length(mu) %>%
    purrr::map(~extraDistr::rtnorm(1,mean = mu[.x],
                                   sd=sqrt(Sigma[.x,.x]),
                                   a = lbvec[.x], b = ubvec[.x])) %>%
    unlist
}


sqrtmatinv <- function(mat) {
  eig <- eigen(mat)
  eig$vec %*% diag(1/sqrt(eig$values), nrow = nrow(mat)) %*% t(eig$vec)

}

sqrtmat <- function(mat) {
  eig <- eigen(mat)
  eig$vec %*%
    diag(sqrt(eig$values), nrow = nrow(mat)) %*%
    t(eig$vec)

}

find_gammas_from_A <- function(A) {
  dims <- dim(A)
  u <- dims[2]
  r <- sum(dims)
  CA <- matrix(0, nrow = r, ncol = u)
  DA <- matrix(0, nrow = r, ncol = r-u)
  CA[(u+1):r, ] <- A
  CA[1:u, 1:u] <- diag(1, u)
  DA[1:u, ] <- -t(A)
  DA[-(1:u), ] <- diag(1, r-u)
  CAtCA <- crossprod(CA)
  DAtDA <- crossprod(DA)
  gamma <- CA %*% sqrtmatinv(CAtCA)
  gamma0 <- DA %*% sqrtmatinv(DAtDA)

  list(gamma = gamma, gamma0 = gamma0,
       CA = CA, CAtCA = CAtCA,
       DA = DA, DAtDA = DAtDA)
}

find_A_from_gamma <- function(gamma) {
  m <- ncol(gamma)
  G1 <- as.matrix(gamma[1:m, ])
  # check if G1 is invertible - else reorganize the predictors
  if (abs(det(G1)) < 1e-7) {
    gamma.t <- t(gamma)
    X.order <- qr(gamma.t, tol = 1e-7)$pivot
    #X <- X[, X.order]
    gamma <- gamma[X.order, ]
  }

  G2 <- gamma[-(1:m), ]
  G2 %*% solve(G1)
}


scaletoCor <- function(Beta,Sigma){
  diag(diag(Sigma^(-1/2))) %*% Beta
  #sqrtmatinv(Sigma) %*% Beta
}

