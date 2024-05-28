#mvprobit.R


#'
#' @export
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,2000)
#' mcmc.out <- mvprobit(inputdata)
#'
mvprobit <- function(inputdata,postnum=5000){

  library(magrittr)
  r <- dim(inputdata$Y)[2]
  p <- dim(inputdata$X)[2]+1

  Data1 <- list(p=r , X=cbind(1,inputdata$X) %x% diag(1,r),y=as.vector(t(inputdata$Y)))
  Mcmc1 = list(R=postnum, keep=1)
  ###
  k <- ncol(Data1$X)
  #prior <- list(A=diag(10^(-6),k),V = diag(10^(-6),r))
  prior <- list(A=diag(10^(-6),k),V = diag(10^(0),r))
  ##

  out = bayesm::rmvpGibbs(Data=Data1, Mcmc=Mcmc1,Prior= prior)


  Sigsample <- 1:dim(out$sigmadraw)[1] %>%
    purrr::map(~matrix(out$sigmadraw[.x,],r,r)) %>%
    abind::abind(along=3)

  Betasample <- 1:dim(out$sigmadraw)[1] %>%
    purrr::map(~matrix(out$betadraw[.x,],r,p)) %>%
    abind::abind(along=3)

  list(Betasample= Betasample[,,-(1:as.integer(postnum/2))],
       Sigsample=Sigsample[,,-(1:as.integer(postnum/2))])
}



#mvprobit.R


#'
#' @export
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,200)
#' mcmc.out <- mvprobit_noX(inputdata)
#'
mvprobit_noX <- function(inputdata,postnum=5000){

  library(magrittr)
  r <- dim(inputdata$Y)[2]
  p <- 1

  Data1 <- list(p=r , X=rep(1,dim(inputdata$X)[1]) %x% diag(1,r),y=as.vector(t(inputdata$Y)))
  Mcmc1 = list(R=postnum, keep=1)

  out = bayesm::rmvpGibbs(Data=Data1, Mcmc=Mcmc1)


  Sigsample <- 1:dim(out$sigmadraw)[1] %>%
    purrr::map(~matrix(out$sigmadraw[.x,],r,r)) %>%
    abind::abind(along=3)

  Betasample <- 1:dim(out$sigmadraw)[1] %>%
    purrr::map(~matrix(out$betadraw[.x,],r,p)) %>%
    abind::abind(along=3)

  list(Betasample= Betasample[,,-(1:as.integer(postnum/2))],
       Sigsample=Sigsample[,,-(1:as.integer(postnum/2))])
}




#'
#' @export
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,5000)
#' mcmc.out <- mvprobit(inputdata)
#' newdata <- generate_data(param,100)
#' perf <- perf_pred(newdata,mcmc.out)
#'
perf_pred_noX <- function(newdata,mcmc.out){
  library(pROC)
  mcmc.num <- dim(mcmc.out$Sigsample)[3]
  n <- dim(newdata$X)[1]
  r <- dim(newdata$Y)[2]
  probpred <-
    1:dim(mcmc.out$Sigsample)[3] %>%
    purrr::map(~pnorm(0,rep(1,dim(newdata$X)[1]) %*% t(mcmc.out$Betasample[,.x]) ,
                      sd = rep(1,n) %*% t(sqrt(diag(mcmc.out$Sigsample[,,.x]))),
                      lower.tail = FALSE)) %>%
    purrr::reduce(`+`)
  probpred <- probpred/dim(mcmc.out$Sigsample)[3]

  mymldr <- mldr::mldr_from_dataframe(data.frame(cbind(newdata$Y,newdata$X)),
                                      labelIndices = 1:r,
                                      name = "testMLDR")
  predperf <- mldr::mldr_evaluate(mymldr,probpred)
  predperf$auc <- predperf$roc$auc
  return(predperf)
}



