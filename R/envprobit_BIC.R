#'
#' @export
#'
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,100)
#' mcmc.env <- envprobit_BIC(inputdata,1:2,mcmc.num=20,ncore=2)
#'
envprobit_BIC <- function(inputdata,uvec,ncore=1,...){
  r <- dim(inputdata$Y)[2]
  p <- dim(inputdata$X)[2]
  n <- dim(inputdata$Y)[1]

  BIC_u <- function(u){
    mcmc.out <- envprobit(inputdata,u,...)
    postlik <- 1:dim(mcmc.out$Betasample)[3] %>%
      purrr::map(function(ind){
        loglik_probit(mcmc.out$Betasample[,,ind],
                      mcmc.out$Sigsample[,,ind],
                      inputdata)}) %>% unlist

    BICMCMC <- -2*max(postlik) + (r*(r + 1)/2 + r + p*u-1)*log(n)

    Dbar <- -2 * mean(postlik)
    Dhat <- -2 * loglik_probit(apply(mcmc.out$Betasample,c(1,2),mean),
                               apply(mcmc.out$Sigsample,c(1,2),mean),
                               inputdata)
    DIC <- 2*Dbar - Dhat
    data.frame(u=u,BIC=BICMCMC,DIC=DIC,IC=3*Dbar - Dhat)
  }
  if(ncore>1){
    library(furrr)
    plan(multisession, workers = ncore)
    BICs <- uvec %>% furrr::future_map(~BIC_u(.x),
                        .options = furrr_options(seed = TRUE))

  }else BICs <- uvec %>% purrr::map(~BIC_u(.x))
  do.call("rbind",BICs)

}



#'
#' @export
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,100)
#' BICres <- calc_BIC(inputdata,1:2)
#'
calc_BIC <- function(inputdata,uvec){
  p <- dim(inputdata$X)[2]
  r <- dim(inputdata$Y)[2]
  n <- dim(inputdata$X)[1]

  MLEres <- uvec %>% purrr::map(~MCEM_envprobit(inputdata,.x))
  BICvec <-
    purrr::map2(uvec,MLEres,
                ~ -2*loglik_probit(cbind(.y$mu,.y$beta),.y$Sigma,inputdata) +
                  (r*(r + 1)/2 + r + p*.x-1)*log(n)) %>% unlist
  data.frame(u = uvec, BIC = BICvec)

}


#'
#' @export
#'
#' @examples
#'
#' r <- 3
#' p <- 5
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,50)
#' loglik_probit(cbind(param$mu.tru,param$beta.tru),param$Sigma.tru,inputdata)
#'
loglik_probit <- function(mubeta,Sigma,inputdata){
  p <- dim(inputdata$X)[2]
  r <- dim(inputdata$Y)[2]

  colnames(inputdata$X) <- paste0("x",1:p)
  xMat <- cbind(const=1,inputdata$X)

  yMat <- inputdata$Y
  colnames(yMat) <- paste0("y",1:r)
  formula_str <- paste0("cbind(",paste(colnames(yMat),collapse = ","),")",
                        "~",paste(paste0("x",1:p),collapse = "+"))

  Sigma <- LaplacesDemon::as.positive.definite(Sigma)
  sum(mvProbit::mvProbitLogLik(as.formula(formula_str),data=as.data.frame(cbind(xMat,yMat)),
                               coef=c(t(scaletoCor(mubeta,Sigma))),sigma=cov2cor(Sigma)))
}

