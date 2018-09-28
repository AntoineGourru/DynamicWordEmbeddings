train <- function(X,model){
  
  time_periods <- length(X)
  
  for (t in 1:time_periods) {
    A <- X[[t]]
    model <- optim(A,model)
  }
  return(model)
}

optim <- function(X_t,model){
  t <- model$timePosition
  
  if(t == 1){
    model$U[[1]] <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
    model$V[[1]] <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  }else{
    model$U[[t]] <- model$U[[t-1]]
    model$V[[t]] <- model$V[[t-1]]
  }
  
  for (i in 1:model$nb_epochs) {
    
    model$U[[t]] <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
    model$V[[t]] <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
    likelyhood_t(X_t,model)
  }
  
  model$timePosition <- model$timePosition + 1 
  return(model)
}