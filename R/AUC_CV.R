#'
#' @export
#'
#' @examples
#' library(magrittr)
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' traindata <- generate_data(param,200)
#' testdata <- generate_data(param,20)
#' init <- MCEM_envprobit(traindata,u) %>% find_A_from_init
#' rocres <- microAUC(traindata,testdata,init,u)
#'
microAUC <- function(traindata,testdata,u,init=NULL,mcmc.num=2000){
  testn <- dim(testdata$X)[1]

  #train과 init으로 fitting
  mcmc.env <- envprobit(traindata,u,init = init, mcmc.num=mcmc.num)
  perf_pred(testdata,mcmc.env)
  #probpred <-
  #  1:dim(mcmc.env$Betasample)[3] %>%
  #  purrr::map(~pnorm(0,cbind(1,testdata$X) %*% t(mcmc.env$Betasample[,,.x]) ,
  #                    sd = rep(1,testn) %*% t(sqrt(diag(mcmc.env$Sigsample[,,.x]))),
  #                    lower.tail = FALSE)) %>%
  #  purrr::reduce(`+`)
  #probpred <- probpred/dim(mcmc.env$Betasample)[3]

  #pROC::roc(as.integer(testdata$Y),as.numeric(probpred))

}


#'
#' @export
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,200)
#' auccv <- AUC_CV(u,inputdata,ncore=10)
#'
#'
AUC_CV <- function(u,inputdata,ncore=1,CVindex=NULL,mcmc.num=2000){
  library(magrittr)
  #init계산
  #init <- MCEM_envprobit(inputdata,u) %>% find_A_from_init
  init = NULL
  #train test 나누기
  if(is.null(CVindex)){
    CVindex <- caret::createMultiFolds(1:dim(inputdata$X)[1], k = 10, times = 5)
  }

  aucfun <- function(index){
    traindata <- list(X=inputdata$X[index,],Y=inputdata$Y[index,])
    testdata <- list(X=inputdata$X[-index,],Y=inputdata$Y[-index,])
    microAUC(traindata,testdata,init=init,u=u,mcmc.num=mcmc.num) %>%
      pROC::auc()
  }
  #병렬적으로 auc 계산 후 평균 return
  if(ncore>1){
    library(furrr)
    plan(multisession, workers = ncore)
    aucvec <- CVindex %>%
      furrr::future_map(~aucfun(.x),
                        .options = furrr_options(seed = TRUE)) %>% unlist
  } else{
    aucvec <- CVindex %>%
      purrr::map(~aucfun(.x)) %>% unlist
  }

  return(aucvec)
}


#'
#' @export
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,200)
#' auccv <- wrapAUC_CV(1:2,inputdata,ncore=10)
#'
#'
wrapAUC_CV <- function(uvec,inputdata,ncore=1,times=1,mcmc.num=2000){
  CVindex <- caret::createMultiFolds(1:dim(inputdata$X)[1], k = 10, times = times)
  aucvec <- uvec %>%
    purrr::map(~AUC_CV(.x,inputdata,ncore,CVindex,mcmc.num=mcmc.num)) %>%
    purrr::map(~mean(.x)) %>% unlist

  data.frame(u=uvec,auc=aucvec)
}
