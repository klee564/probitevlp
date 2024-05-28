
#'
#' @export
#'
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,50)
#' mcmc.env <- envprobit(inputdata,u,mcmc.num=2000)
#'
envprobit <- function(inputdata,u,mcmc.num=100,init=NULL,hyper=NULL,
                      burnin=as.integer(mcmc.num/2)){
  library(magrittr)

  #pre-define
  postlist <- list()
  Y <- inputdata$Y
  X <- inputdata$X
  #Wsam <- inputdata$Z
  r <- dim(inputdata$Y)[2]
  p <- dim(inputdata$X)[2]
  n <- dim(X)[1]

  lbmat <- -1/Y +1
  ubmat <- (1/!Y)-1

  barX <- colMeans(X)
  Xc <- scale(X,scale=FALSE)
  XctXc <- crossprod(Xc)
  barXbarXt <- crossprod(t(barX))

  #set hyperparam
  if(is.null(hyper)){
    Psi <- diag(10^(-0),u)
    nu <- u
    Psi0 <- diag(10^(-0),r-u)
    nu0 <- r-u
    A0 <- matrix(0,r-u,u)
    #K <- diag(10^4, r-u)
    #L <- diag(10^4, u)
    K <- diag(10^6, r-u)
    L <- diag(10^6, u)
    M <- diag(10^(-6),p)
    kappa0 <- 10^(-6)
    K.half.inv <- sqrtmatinv(K)
    L.half.inv <- sqrtmatinv(L)
    M.half <- sqrtmat(M)
  } else{
    #hyper <- list(Psi = diag(10^(-8),u),nu = u, Psi0 = diag(10^(-8),r-u),
    #              nu0 = r-u, A0 = matrix(0,r-u,u), K = diag(10^8, r-u),
    #              L = diag(10^8, u), M = diag(10^(-8),p), kappa0 = 10^(-8))
    Psi <- diag(hyper$Psi,u)
    nu <- hyper$nu
    Psi0 <- diag(hyper$Psi0,r-u)
    nu0 <- hyper$nu0
    A0 <- hyper$A0
    K <- diag(hyper$K, r-u)
    L <- diag(hyper$L, u)
    M <- diag(hyper$M,p)
    kappa0 <- hyper$kappa0
    K.half.inv <- sqrtmatinv(K)
    L.half.inv <- sqrtmatinv(L)
    M.half <- sqrtmat(M)

  }


  #set init
  if(is.null(init)){
    init <- MCEM_envprobit(inputdata,u) %>% find_A_from_init
    print("initial value computing is completed.")
  }

  Omega <- init$Omega
  Omega0 <- init$Omega0
  tmu <- init$mu
  teta <- init$eta
  A <- matrix(init$A,nrow=r-u,ncol=u)

  #A <- matrix(param$A.tru,nrow=r-u,ncol=u)

  #A <- matrix(runif(u * (r - u), min = -1, max = 1), nrow = r-u, ncol = u)

  #A로 부터 유도되는 값들
  gamma_gamma0 <- find_gammas_from_A(A)
  gamma <- gamma_gamma0$gamma
  gamma0 <- gamma_gamma0$gamma0

  CA <- gamma_gamma0$CA
  CAtCA <- gamma_gamma0$CAtCA
  CA_CAtCAinv <- CA %*% solve(CAtCA)

  Sigma <- gamma %*% Omega %*% t(gamma) + gamma0 %*% Omega0 %*% t(gamma0) ##위치 변경

  #barW <- colMeans(Y)
  #Wsam <- scale(Y,scale=FALSE)

  rw_var <- rep(1,dim(A)[2])
  for(iter in 1:mcmc.num){


    ## W sample

    #Dinv <- diag(1/sqrt(diag(Sigma)))
    #meanmat <- t(as.vector(tmu) + CA %*% teta %*% t(X)) %*% Dinv
    #R <- LaplacesDemon::as.positive.definite(Dinv %*% Sigma %*% Dinv)
    #condparam <- condparam_W(R)
    meanmat <- t(as.vector(tmu) + CA %*% teta %*% t(X))
    condparam <- condparam_W(Sigma)

    Wsam <- gen_Wmat(meanmat,condparam,lbmat,ubmat)

    #mumat <- t(as.vector(tmu) + CA %*% teta %*% t(X))
    #Wsam <- gen_Wmat2(mumat,LaplacesDemon::as.positive.definite(Sigma),lbmat,ubmat)



    barW <- colMeans(Wsam)
    Wc <- scale(Wsam,scale=FALSE)

    ## Omega,Omega0 sample
    C1 <- Wc %*% CA_CAtCAinv
    C2 <- t(barW)%*% CA_CAtCAinv
    tbareta <- solve(XctXc + n*kappa0/(n+kappa0)* barXbarXt + M) %*%
      (t(Xc)%*%C1 + n*kappa0/(n+kappa0)*barX %*% C2) ##solve(A)%*%B 개선

    C3 <- crossprod(Xc %*%tbareta-C1) + t(tbareta) %*% M %*% tbareta +
      n*kappa0/(n+kappa0)*crossprod(t(barX) %*% tbareta-C2) ##tbareta 개선

    Omega <-
      CholWishart::rInvWishart(1,n+nu,Psi+ sqrtmat(CAtCA) %*%
                                                C3 %*% sqrtmat(CAtCA))[,,1]

    Omega0 <-
      CholWishart::rInvWishart(1,n+nu0,Psi0 + t(gamma0) %*% (crossprod(Wc) +
                  n*kappa0/(n+kappa0)*crossprod(t(barW)) ) %*%gamma0)[,,1]

    ## eta sample
    teta <-LaplacesDemon::rmatrixnorm(t(tbareta),
                                      LaplacesDemon::as.positive.definite(sqrtmatinv(CAtCA)  %*% Omega %*% sqrtmatinv(CAtCA)),
                                      LaplacesDemon::as.positive.definite(solve(XctXc + n*kappa0/(n+kappa0)*barXbarXt +M))) #개선여지
    #teta <-LaplacesDemon::rmatrixnorm(t(tbareta),
    #                                  LaplacesDemon::as.positive.definite(as.matrix(Matrix::nearPD(sqrtmatinv(CAtCA)  %*% Omega %*% sqrtmatinv(CAtCA),posd.tol = 0.01)$mat)),
    #                                  LaplacesDemon::as.positive.definite(as.matrix(solve(Matrix::nearPD(XctXc + n*kappa0/(n+kappa0)*barXbarXt +M,posd.tol = 0.01)$mat))))

    ## mu sample
    Sigma <- gamma %*% Omega %*% t(gamma) + gamma0 %*% Omega0 %*% t(gamma0)
    tmu <- mvnfast::rmvn(1,mu=n/(n+kappa0)*(barW-CA %*% teta %*% barX) ,
                         sigma= Sigma/(n+kappa0))



    ## A sample and gamma 유도
    #A sample
    Omega0.inv <- solve(Omega0)
    Omega.inv <- solve(Omega)
    lpd_val <- lpd_A_pred(A, Wsam, X,teta,tmu,Omega0.inv,Omega.inv,
                          K.half.inv, L.half.inv, A0,M.half)
    MHsample <- rwmh_colwise(lpd_val, rw_var,A,Wsam, X,teta,tmu,Omega0.inv,
                             Omega.inv,K.half.inv, L.half.inv, A0,M.half)



    #tune_int <- iter-as.integer(mcmc.num/4)
    tune_int <- iter
    if(tune_int>0) rw_var <- rw_var * exp(tune_int^(-0.7)*(MHsample$alphas-0.44))
    A <- MHsample$A

    gamma_gamma0 <- find_gammas_from_A(A)
    gamma <- gamma_gamma0$gamma
    gamma0 <- gamma_gamma0$gamma0

    CA <- gamma_gamma0$CA
    CAtCA <- gamma_gamma0$CAtCA
    CA_CAtCAinv <- CA %*% solve(CAtCA)


    vtilde <- sum(diag(Omega)) +sum(diag(Omega0))
    Omega <- Omega/vtilde
    Omega0 <- Omega0/vtilde
    teta <- teta/sqrt(vtilde)
    tmu <- tmu/sqrt(vtilde)

    Sigma <- gamma %*% Omega %*% t(gamma) + gamma0 %*% Omega0 %*% t(gamma0) ##위치 변경



    postlist[[iter]] <- list(Sigma=Sigma,Omega=Omega,Omega0=Omega0,teta=teta,
                             tmu=tmu,A=A,alphas=MHsample$alphas)

    if(iter %% 100==0) cat(iter,"\n")


  }

  return(envprobit_after(postlist,burnin))
}

#'
#'
#'
envprobit_after <- function(mcmc.env,burnin){
  #postnum <- length(mcmc.env)
  list(Betasample = mcmc.env[-(1:burnin)] %>%
         purrr::map(~cbind(as.vector(.x$tmu),
                           rbind(diag(1,dim(.x$A)[2]),.x$A) %*%.x$teta)) %>%
         abind::abind(along=3),
       Sigsample = mcmc.env[-(1:burnin)] %>%
         purrr::map(~.x$Sigma) %>% abind::abind(along=3),
       Asample = mcmc.env[-(1:burnin)] %>%
         purrr::map(~.x$A) %>% abind::abind(along=3),
       tetasample = mcmc.env[-(1:burnin)] %>%
         purrr::map(~.x$teta) %>% abind::abind(along=3)
       )
}
