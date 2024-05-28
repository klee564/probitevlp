#performance.R

#'
#' @export
#'
#' @examples
#' r <- 5
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,100)
#' mcmc.out <- mvprobit(inputdata)
#' newdata <- generate_data(param,100)
#' perf <- perf_pred(newdata,mcmc.out)
#'
perf_pred2 <- function(newdata,mcmc.out){
  library(pROC)
  testn <- dim(newdata$X)[1]
  #probpred <-
  #  1:dim(mcmc.out$Betasample)[3] %>%
  #  purrr::map(~pnorm(0,cbind(1,as.matrix(newdata$X)) %*%
  #                      t(scaletoCor(mcmc.out$Betasample[,,.x],
  #                                   mcmc.out$Sigsample[,,.x])) ,
  #                   lower.tail = FALSE)) %>%
  #  purrr::reduce(`+`)
  probpred <-
    1:dim(mcmc.out$Betasample)[3] %>%
    purrr::map(~pnorm(0,cbind(1,as.matrix(newdata$X)) %*% t(mcmc.out$Betasample[,,.x]) ,
                      sd = rep(1,testn) %*% t(sqrt(diag(mcmc.out$Sigsample[,,.x]))),
                      lower.tail = FALSE)) %>%
    purrr::reduce(`+`)
  probpred <- probpred/dim(mcmc.out$Betasample)[3]



  postmean <- 1:dim(mcmc.out$Betasample)[3] %>%
    purrr::map(~scaletoCor(mcmc.out$Betasample[,,.x],
                           mcmc.out$Sigsample[,,.x])) %>%
    abind::abind(along=3) %>% apply(c(1,2),mean)

  probpred2 <- pnorm(0,cbind(1,as.matrix(newdata$X)) %*%
                       t(postmean),lower.tail = FALSE)


  list(probpred=probpred,probpred2=probpred2,auc=pROC::auc(as.integer(newdata$Y),as.numeric(probpred)))
}

#'
#' @export
#'
#' @examples
#' r <- 5
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,100)
#' mcmc.out <- mvprobit(inputdata)
#' newdata <- generate_data(param,100)
#' perf <- perf_pred(newdata,mcmc.out)
#'
perf_pred <- function(newdata,mcmc.out){
  library(pROC)
  testn <- dim(newdata$X)[1]
  probpred <-
    1:dim(mcmc.out$Betasample)[3] %>%
    purrr::map(~pnorm(0,cbind(1,as.matrix(newdata$X)) %*% t(mcmc.out$Betasample[,,.x]) ,
                      sd = rep(1,testn) %*% t(sqrt(diag(mcmc.out$Sigsample[,,.x]))),
                      lower.tail = FALSE)) %>%
    purrr::reduce(`+`)
  probpred <- probpred/dim(mcmc.out$Betasample)[3]

  pROC::roc(as.integer(newdata$Y),as.numeric(probpred))
}

#'
#' @export
#'
#' @examples
#' r <- 5
#' p <- 10
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,100)
#' mcmc.out <- mvprobit(inputdata)
#' perf <- perf_coef(param,mcmc.out)
#' mean(perf$coverage); mean(mean(perf$pointwise))
#' mcmc.out2 <- envprobit(inputdata,2,mcmc.num=2000)
#' perf <- perf_coef(param,mcmc.out2)
#' mean(perf$coverage); mean(mean(perf$pointwise))
#'
#'
perf_coef <- function(param.tru,mcmc.out){
  Beta.tru <- scaletoCor(cbind(param.tru$mu.tru,param.tru$beta.tru),
                         param.tru$Sigma.tru)

  postmean <- 1:dim(mcmc.out$Betasample)[3] %>%
    purrr::map(~scaletoCor(mcmc.out$Betasample[,,.x],
                           mcmc.out$Sigsample[,,.x])) %>%
    abind::abind(along=3) %>% apply(c(1,2),mean)

  postarray <- 1:dim(mcmc.out$Betasample)[3] %>%
    purrr::map(~scaletoCor(mcmc.out$Betasample[,,.x],
                           mcmc.out$Sigsample[,,.x])) %>%
    abind::abind(along=3)

  postlower <- postarray %>% apply(c(1,2),function(x){quantile(x,prob=c(0.025))})
  postupper <- postarray %>% apply(c(1,2),function(x){quantile(x,prob=c(0.975))})



  list(pointwise = (postmean - Beta.tru)^2,
       spnorm = norm(postmean - Beta.tru,type="2"),
       postmean = postmean,
       coverage = postlower <= Beta.tru & postupper >= Beta.tru)
}


#'
#' @export
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,50)
#' mcmc.out <- mvprobit(inputdata,postnum=50)
#' newdata <- generate_data(param,50)
#' perf <- perf_LPPD(newdata,mcmc.out)
#'
perf_LPPD <- function(newdata,mcmc.out){
  postlik <- 1:dim(mcmc.out$Betasample)[3] %>%
    purrr::map(function(ind){
      loglik_probit(mcmc.out$Betasample[,,ind],
                    mcmc.out$Sigsample[,,ind],
                    newdata)}) %>% unlist

  matrixStats::logSumExp(postlik) - log(length(postlik))
}

