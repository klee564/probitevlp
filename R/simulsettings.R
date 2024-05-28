

#'
#'
#' @export
#'
#' @examples
#'
#' r <- 10
#' p <- 10
#' u <- u.tru <- 5
#' all_pars <- generate_par(r, p, u)
#'
#'
generate_par <- function(r, p, u,leta=0,ueta=1,isuni=TRUE,iscor=TRUE) {
  mu.tru <- runif(r, -0, 0)
  if(isuni){
    eta.tru <- matrix(runif(u * p, min =leta, max = ueta), nrow = u, ncol = p) #2/sqrt(u)
  } else{
    eta.tru <- matrix(rnorm(u * p,  sd = 2), nrow = u, ncol = p) #2/sqrt(u)
  }
  A <- A.tru <- matrix(runif(u * (r - u), min = -1, max = 1), nrow = r-u, ncol = u)
  gamma_gamma0 <- find_gammas_from_A(A)
  gamma <- gamma_gamma0$gamma
  gamma0 <- gamma_gamma0$gamma0
  CA <- gamma_gamma0$CA

  # CA <- rbind(diag(1, u), A)
  # DA <- rbind(-t(A), diag(1, r-u))
  #
  # calculate beta
  beta.tru <- CA %*% eta.tru
  #beta.tru <- gamma %*% eta.tru
  #beta.tru <- rstiefel::rmf.matrix(diag(1,r))[,1:u] %*% eta.tru

  V.tru <- diag(sort(runif(u, 0.5,1), decreasing = TRUE),
                ncol = u, nrow = u)
  V0.tru <- diag(sort(runif(r-u, 9, 10), decreasing = TRUE),
                 ncol = r-u, nrow = r-u)

  if(iscor){
    R.tru <- make_corr(u)
    R0.tru <- make_corr(r-u)

    Omega.tru <- sqrt(V.tru) %*% R.tru %*% sqrt(V.tru)
    Omega0.tru <- sqrt(V0.tru) %*% R0.tru %*% sqrt(V0.tru)

  }else{
    Omega.tru <- V.tru
    Omega0.tru <- V0.tru
  }


  Sigma1 <- gamma %*% Omega.tru %*% t(gamma)
  Sigma2 <- gamma0 %*% Omega0.tru %*% t(gamma0)

  Sigma.tru <- Sigma1 + Sigma2

  #scaled_mubeta <- scaletoCor(cbind(mu.tru,beta.tru),Sigma.tru)
  #mu.tru <- scaled_mubeta[,1]
  #beta.tru <- scaled_mubeta[,-1]
  #Sigma.tru <- cov2cor(Sigma.tru)


  list(mu.tru = mu.tru,
       eta.tru=eta.tru,
       beta.tru = beta.tru,
       Omega.tru = Omega.tru,
       Omega0.tru = Omega0.tru,
       #A.tru = A.tru,
       gamma.tru = gamma,
       #gamma0.tru = gamma0,
       Sigma.tru = Sigma.tru,
       mux = 0,
       sigmax = 1)
}


#'
#'
#'
#'
make_corr <- function(d){
  if(d==1) return(1) else{
    Sigma=diag(1,d)
    Sigma <- 0.5^(abs(row(Sigma)-col(Sigma)))
    #return(cov2cor(rWishart(1, 10*d, Sigma/(10*d))[,,1]))
    return(cov2cor(Sigma))
  }

}




#'
#' @export
#'
#' @examples
#'
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' simuldata <- generate_data(param,50)
#'
#'
generate_data <- function(param,n,  ...)
{
  mux <- param$mux
  sigmax <- param$sigmax
  mu.tru <- param$mu.tru
  Sigma.tru <- param$Sigma.tru
  beta.tru <- param$beta.tru
  p <- dim(beta.tru)[2]
  r <- dim(beta.tru)[1]



  X <- matrix(rnorm(n*p, mean = mux, sd = sigmax),
              nrow = n, ncol = p)
  eps <- matrix(rnorm(n*r), n, r) %*% sqrtmat(Sigma.tru)
  Z <- tcrossprod(rep(1, n), mu.tru) + X %*% t(beta.tru) + eps
  Y <- apply(Z>0,2,as.integer)

  list(X = X, Y = Y,Z=Z)
}

