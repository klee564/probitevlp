
#'
#' @export
#'
#' @examples
#' r <- 6
#' p <- 10
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,200)
#' param$A.tru
#' res <- MCEM_envprobit(inputdata,u)
#' find_gammas_from_A(find_A_from_gamma(res$m$Gamma))$CA
#' m <- res$m
#' mcmc.out <- res$mcmc.out
#' Sigmean <- apply(mcmc.out$Sigsample,c(1,2),mean)
#' sqrt(mean((scaletoCor(m$beta,m$Sigma)-scaletoCor(param$beta.tru,param$Sigma.tru))^2))
#' sqrt(mean((scaletoCor(apply(mcmc.out$Betasample,c(1,2),mean),Sigmean)[,-1]-scaletoCor(param$beta.tru,param$Sigma.tru))^2))
#'
MCEM_envprobit <- function(inputdata,u,niter=25){
  Y <- inputdata$Y
  X <- inputdata$X
  lbmat <- -1/Y +1
  ubmat <- (1/!Y)-1

  m <- Renvlp::env(X, Y+runif(length(Y),-0.01,0.01), u)
  for(i in 1:niter){

    meanmat <- t(as.vector(m$mu) + m$beta %*% t(X))
    condparam <- condparam_W(m$Sigma)

    Wmat <- purrr::rerun(50,gen_Wmat(meanmat,condparam,lbmat,ubmat)) %>% do.call("rbind",.)
    #Wmat <- purrr::rerun(50,gen_Wmat2(meanmat,m$Sigma,lbmat,ubmat)) %>% do.call("rbind",.)
    Xmat <- purrr::rerun(50,X) %>% do.call("rbind",.)
    m <- Renvlp::env(Xmat, Wmat, u)
    #

    #print(i)
    #print(norm(scaletoCor(m$beta,m$Sigma)-scaletoCor(param$beta.tru,param$Sigma.tru),type="2"))
    #print(find_gammas_from_A(find_A_from_gamma(m$Gamma))$CA)
  }
  #find_A_from_init(m)
  m

#  mcmc.out <- mvprobit(inputdata)
#  Sigmean <- apply(mcmc.out$Sigsample,c(1,2),mean)
#  norm(scaletoCor(apply(mcmc.out$Betasample,c(1,2),mean),Sigmean)[,-1]-scaletoCor(param$beta.tru,param$Sigma.tru),type="2")
}


#'
#'
#'
find_A_from_init <- function(m){
  gg0 <- find_gammas_from_A(find_A_from_gamma(m$Gamma))
  u <- dim(m$Gamma)[2]

  list(mu = m$mu,eta = m$Gamma[1:u,] %*% m$eta, A= gg0$CA[-(1:u),],
       Omega =  t(gg0$gamma) %*% m$Sigma %*% gg0$gamma,
       Omega0 = t(gg0$gamma0) %*% m$Sigma %*% gg0$gamma0)
}
