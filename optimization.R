train <- function(X,model){
  
  time_periods <- length(X)
  
  for (t in 1:time_periods) {
    A <- X[[t]]
    model <- optim(A,model)
  }
  return(model)
}

optim <- function(X_t,model){
  learning_rate <- 0.02
  t <- model$timePosition
  
  if (t == 1) {
    model$U[[1]] <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
    model$V[[1]] <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  }else{
    model$U[[t]] <- model$U[[t - 1]]
    model$V[[t]] <- model$V[[t - 1]]
  }
  
  for (i in 1:model$nb_epochs) {
    for (i in model$D) {
      
      
      vp <- draw_VP(model)
      
      #Learning, on peut mettre ici une condition d'arret
      
      #BLACKBOX VI pour U
      Us <- mvrnorm(n = 1,vp$u_mu,vp$u_sigma)
      grad_mu_ui <- 1
      grad_sigma_ui <- 1
      for (s in 2:nb_sampleVI) {
        Us <- mvrnorm(n = 1,vp$u_mu,vp$u_sigma)
        temp_grad <- 1
        grad_mu_ui <- grad_mu_ui + temp_grad
        grad_sigma_ui <- grad_sigma_ui + temp_grad
      }
      grad_mu_ui <- grad_mu_ui / nb_sampleVI
      grad_sigma_ui <- grad_sigma_ui / nb_sampleVI
      
      vp$u_mu <- vp$u_mu + learning_rate * grad_mu_ui
      vp$u_sigma <- diag(vp$u_sigma) + learning_rate * grad_sigma_ui
      
      #BLACKBOX VI pour V
      Vs <- mvrnorm(n = 1,vp$v_mu,vp$v_sigma)
      grad_mu_vi <- 1
      grad_sigma_vi <- 1
      for (s in 2:nb_sampleVI) {
        Vs <- mvrnorm(n = 1,vp$v_mu,vp$v_sigma)
        temp_grad <- 1
        grad_mu_vi <- grad_mu_vi + temp_grad
        grad_sigma_vi <- grad_sigma_vi + temp_grad
      }
      grad_mu_vi <- grad_mu_vi / nb_sampleVI
      grad_sigma_vi <- grad_sigma_vi / nb_sampleVI
      
      vp$v_mu <- vp$v_mu + learning_rate * grad_mu_vi
      vp$v_sigma <- diag(vp$v_sigma) + learning_rate * grad_sigma_vi
      
      model$U[[t]][i,] <- mvrnorm(n = 1,vp$u_mu,vp$u_sigma)
      model$V[[t]][i,] <- mvrnorm(n = 1,vp$v_mu,vp$v_sigma)
    }
    
    likelyhood_t(X_t,model)
  }
  model$timePosition <- model$timePosition + 1 
  return(model)
}

draw_VP <- function(model){
  vp$u_mu <- mvrnorm(n = 1, rep(0,model$K),diag(rep(1,model$K)))
  vp$v_mu <- mvrnorm(n = 1, rep(0,model$K),diag(rep(1,model$K)))
  s <- rnorm(0,1)
  vp$u_sigma <- diag(rep(s,model$K))
  vp$v_sigma <- diag(rep(s,model$K))
  return(vp)
}

ProdNorm_q <- function(the_i,model){
  sum <- 1
  t <- model$timePosition
  for (i in 1:D) {
    if (i != the_i) {
      un <- pnorm(U[[t]][i,],model$VP_u_mu[[t]][i,],diag(model$VP_u_sigma[[t]][i,]))
      vn <- pnorm(V[[t]][i,],model$VP_v_mu[[t]][i,],diag(model$VP_v_sigma[[t]][i,]))
      sum <- sum + un + vn
    }
  }
}