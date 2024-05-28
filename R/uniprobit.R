

#'
#' @export
#'
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,200)
#' mcmc.uni <- wrap_uniprobit(inputdata)
#'
wrap_uniprobit <- function(inputdata,mcmc.num=100,hyper=NULL,
                           burnin=as.integer(mcmc.num/2)){

  Y <- inputdata$Y
  r <- dim(Y)[2]

  Betasample <- abind::abind(
    purrr::map(1:r,~uniprobit(Y[,.x],inputdata$X,mcmc.num=mcmc.num,hyper=hyper)),
    along = 3) %>% aperm(3:1)
  list(Betasample=Betasample,
       Sigsample = purrr::map(1:dim(Betasample)[3],~diag(1,r)) %>% abind::abind(along=3))

}

#'
#' @export
#'
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,200)
#' mcmc.uni <- uniprobit(y=inputdata$Y[,1],X=inputdata$X,mcmc.num=2000)
#'
uniprobit <- function(y,X,mcmc.num=100,hyper=NULL,
                      burnin=as.integer(mcmc.num/2)){
  library(magrittr)

  postlist <- list()

  p <- dim(X)[2]
  n <- dim(X)[1]
  lbvec <- -1/y +1
  ubvec <- (1/!y)-1

  if(is.null(hyper)){
    M <- diag(10^(-6),p+1)
    M.half <- sqrtmat(M)
  }

  beta <- rep(0,p+1)

  V <- solve(crossprod(cbind(1,X)))

  for (iter in 1:mcmc.num) {
    # Update Mean of z
    mu_z <- cbind(1,X) %*% beta
    # Draw latent variable z from its full conditional: z | \theta, y, X
    z <- truncnorm::rtruncnorm(1, mean = mu_z, sd = 1, a = lbvec, b = ubvec)

    # Compute posterior mean of theta
    M <- V %*% crossprod(cbind(1,X), z)
    # Draw variable \theta from its full conditional: \theta | z, X
    beta <- c(mvnfast::rmvn(1, M, V))
    postlist[[iter]] <- beta

  }

  return(do.call("rbind",postlist[-(1:burnin)]))
}
