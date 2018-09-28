init_model <- function(D,K,sigma_t,sigma_0,nb_epochs,nb_MB){
  require(MASS)
  
  model <- list()
  model$sigma_t <- sigma_t
  model$sigma_0 <- sigma_0
  model$K <- K
  model$D <- D
  
  model$nb_MB <- nb_MB
  model$nb_epochs <- nb_epochs
  
  model$U <- list()
  model$V <- list()
  
  model$timePosition <- 1
  return(model)
}


likelyhood_t <- function(X_t,model){
  
}
likelyhood <- function(data,model){
  
}

write <- function(model){
  # TODO
}