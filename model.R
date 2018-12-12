init_model <- function(D,K,sigma_t,sigma_0,nb_epochs,nb_MB,nb_sampleVI){
  require(MASS)
  
  model <- list()
  model$sigma_t <- sigma_t
  model$sigma_0 <- sigma_0
  model$K <- K
  model$D <- D
  
  model$nb_MB <- nb_MB
  model$nb_epochs <- nb_epochs
  model$nb_sampleVI <- nb_sampleVI
  
  model$U <- list()
  model$V <- list()
  
  model$vp_m1 <- draw_infoless_VP(D,K)  
  
  model$timePosition <- 1
  return(model)
}

draw_infoless_VP <- function(D,K){
  vp <- list()
  vp$u_mu <- matrix(0,D,K)
  vp$v_mu <- matrix(0,D,K)
  vp$u_sigma <- matrix(1,D,K)
  vp$v_sigma <- matrix(1,D,K)
  return(vp)
}


likelyhood_t <- function(X_t,model){
  
}
likelyhood <- function(data,model){
  
}

write <- function(model){
  # TODO
}