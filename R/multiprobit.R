
#'
#' @export
#'
gen_W2 <- function(muvec,covmat,lbvec,ubvec){
  #d <- length(Zvec)
  #for(i in 1:d){
  #  Zvec[i] <- condTruncMVN::rcmvtruncnorm(1, mean = muvec,sigma = covmat,
                                           #lower = lbvec,upper = ubvec,dependent.ind = i,
                                           #given.ind = (1:d)[-i], X.given = Zvec[-i])
  #}
  tmvtnorm::rtmvnorm(1,mean=muvec,sigma=covmat,lower=lbvec,upper =ubvec,algorithm = "gibbs")
}


#'
#' @export
#'
gen_Wmat2 <- function(mumat,covmat,lbmat,ubmat){
  1:dim(lbmat)[1] %>%
    purrr::map(~gen_W2(muvec=mumat[.x,],
                       covmat=LaplacesDemon::as.positive.definite(covmat),
                       lbvec=lbmat[.x,],ubvec=ubmat[.x,]))%>%
    do.call("rbind",.)
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
#' mcmc.env <- envprobit(inputdata,u,mcmc.num=2000)
#'
multiprobit <- function(inputdata,mcmc.num=1000,init=NULL,hyper=NULL,
                        burnin=as.integer(mcmc.num/2)){
  library(magrittr)

  #pre-define
  postlist <- list()
  Y <- inputdata$Y
  r <- dim(inputdata$Y)[2]

  X <- cbind(1,inputdata$X)
  #Wsam <- inputdata$Z
  p <- dim(inputdata$X)[2]
  n <- dim(X)[1]

  lbmat <- -1/Y +1
  ubmat <- (1/!Y)-1

  #set hyperparam
  Psi_inv <- diag(10^(-6),p+1)



  init <- Renvlp::env(X[,-1], Y, r)
  Sigma <- init$Sigma
  gamma <- cbind(init$mu,init$beta)
  Dinv <- diag(1/sqrt(diag(Sigma)))
  R <- Dinv %*% Sigma %*% Dinv
  beta <- Dinv %*% gamma
  Zmat <- scale(Y,scale=FALSE)

  for(iter in 1:mcmc.num){


    ## W sample
    #mumat <- X %*% t(beta)
    #Zmat <- gen_Wmat2(Zmat,mumat,R,lbmat,ubmat)
    meanmat <- X %*% t(beta)
    condparam <- condparam_W(R)

    Zmat <- gen_Wmat(meanmat,condparam,lbmat,ubmat)


    ##gen di^2
    D <- diag(solve(R)) %>%
      purrr::map(~1/sqrt(rgamma(1,(r+1)/2 , scale=2/.x))) %>% unlist %>%
      diag

    Wmat <- Zmat %*% D
    InvTheta <- crossprod(X) + Psi_inv
    M <- solve(InvTheta, crossprod(X,Wmat))
    Sigma <-
      CholWishart::rInvWishart(1,n+2,
                               LaplacesDemon::as.positive.definite(diag(1,r)+crossprod(Wmat) - t(M) %*% InvTheta %*% M ))[,,1]
    gamma <- LaplacesDemon::rmatrixnorm(t(M),LaplacesDemon::as.positive.definite(Sigma),
                                        LaplacesDemon::as.positive.definite(  solve(InvTheta) )
    ) #개선여지

    Dinv <- diag(1/sqrt(diag(Sigma)))
    beta <- Dinv %*% gamma
    R <- LaplacesDemon::as.positive.definite(Dinv %*% Sigma %*% Dinv)

    postlist[[iter]] <- list(Sigma=R, beta=beta)

    if(iter %% 100==0) cat(iter,"\n")


  }

  return(mulprobit_after(postlist,burnin))
}



#'
#'
#'
mulprobit_after <- function(mcmc.env,burnin){
  #postnum <- length(mcmc.env)
  list(Betasample = mcmc.env[-(1:burnin)] %>%
         purrr::map(~.x$beta) %>%
         abind::abind(along=3),
       Sigsample = mcmc.env[-(1:burnin)] %>%
         purrr::map(~.x$Sigma) %>% abind::abind(along=3))
}

