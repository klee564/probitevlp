#'
#' @export
#' @examples
#' genedataEG <- readRDS("EGFR.rds")
#' procdata <- preproc_gene(genedata)
#'
#'
preproc_gene <- function(genedata,method="AG",PCscale=FALSE){
  rnadata <- genedata$rnadata
  mutdata <- genedata$mutdata
  responsedata <- genedata$responsedata

  NC <- dim(rnadata)[2]
  spca <- ClassDiscovery::SamplePCA(t(scale(rnadata)))
  if(method=="BS"){
    lambda <- spca@variances[1:(NC-1)]
    num_PC <- PCDimension::bsDimension(lambda)
  } else if(method=="AG"){
    ag.obj <- PCDimension::AuerGervini(spca)
    num_PC <- PCDimension::agDimension(ag.obj)
  }

  prres <- prcomp(rnadata,center = T,scale. = T)
  fullX <- data.frame(mutation= mutdata,prres$x[,1:num_PC])
  fullY <- responsedata %>% purrr::map(~as.integer(.x< -2)) %>%
    as.data.frame %>% as.matrix
  scaledX <- scale(fullX)
  if(PCscale) fullX <- scaledX
  ##scale(rnadata) %*% prres$rotation[,1:num_PC] %*% diag(1/apply(fullX,2,sd)[-1]) - scaledX[,-1]
  ##scale(rnadata) %*% VL %*% diag(1/PCscale[-1]) %*% t(beta)
  list(fullZ = responsedata,fullX=fullX,fullY=fullY,VL=prres$rotation[,1:num_PC],PCscale=apply(fullX,2,sd))


}
