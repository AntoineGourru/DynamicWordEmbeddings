
library(mixtools)

train <- function(X,model){
  
  time_periods <- length(X)
  
  for (t in 1:time_periods) {
    A <- X[[t]]
    model <- optim(A,model)
  }
  return(model)
}

optim <- function(X_t,model,color){
  logLout <- c()
  learning_rate <- seq(5,2,length.out = model$nb_epochs)
  # print(learning_rate)
  t <- model$timePosition
  
  vp <- draw_VP(model) 
  
  for (epo in 1:model$nb_epochs) {
    # epo <- 1
    grad_mu_u <- matrix(0,model$D,model$K)
    grad_mu_v <- matrix(0,model$D,model$K)
    grad_sigma_u <- matrix(0,model$D,model$K)
    grad_sigma_v <- matrix(0,model$D,model$K)
    
    
    for (s in 1:model$nb_sampleVI) {
      # Us <- mvrnorm(n = 1,vp$u_mu[1,],diag(vp$u_sigma[1,]))
      # Vs <- mvrnorm(n = 1,vp$v_mu[1,],diag(vp$v_sigma[1,]))
      # for (i in 2:model$D) {
      #   Us <- rbind(Us,mvrnorm(n = 1,vp$u_mu[i,],diag(vp$u_sigma[i,])))
      #   Vs <-  rbind(Vs,mvrnorm(n = 1,vp$v_mu[i,],diag(vp$v_sigma[i,])))
      # }
      Us <- mvrnorm(n = 1,vp$u_mu[1,],diag(rep(1,model$K)))
      Vs <- mvrnorm(n = 1,vp$v_mu[1,],diag(rep(1,model$K)))
      for (i in 2:model$D) {
        Us <- rbind(Us,mvrnorm(n = 1,vp$u_mu[i,],diag(rep(1,model$K))))
        Vs <-  rbind(Vs,mvrnorm(n = 1,vp$v_mu[i,],diag(rep(1,model$K))))
      }
      # PQ <- probaQ(Us,Vs,vp,model$D)
      LL <- likely(X_t,Us,Vs,model$D)
      lPC <- LprobaCond(Us,Vs,model)
      lPQ <- logprobaQ(Us,Vs,vp,model$D)
      
      for (i in 1:model$D) {    
        
        #Learning, on peut mettre ici une condition d'arret
        
        #BLACKBOX VI pour U
        
        # #To verify
        # temp_gradU <- (2*(vp$u_m[i,] - Us[i,])/vp$u_sigma[i,]) * (LL + PC - 2 * lPQ)
        # 
        # grad_mu_u[i,] <-  grad_mu_u[i,] + temp_gradU
        # 
        # temp_gradV <- (2*(vp$v_m[i,] - Vs[i,])/vp$v_sigma[i,]) * (LL + PC - 2 * lPQ)
        # 
        # grad_mu_v[i,] <-  grad_mu_v[i,] + temp_gradV
        
        #To verify
        temp_gradU <- (LL + lPC - lPQ)*(-(vp$u_m[i,] - Us[i,])/vp$u_sigma[i,])
        
        grad_mu_u[i,] <-  grad_mu_u[i,] + temp_gradU
        
        temp_gradV <- (LL + lPC - lPQ)*(-(vp$v_m[i,] - Vs[i,])/vp$v_sigma[i,])
        
        grad_mu_v[i,] <-  grad_mu_v[i,] + temp_gradV
      }
      # grad_sigma <- grad_sigma_ui / model$nb_sampleVI
    }
    grad_mu_u <- grad_mu_u / model$nb_sampleVI
    grad_mu_v <- grad_mu_v / model$nb_sampleVI
    
    grad_mu_u <- grad_mu_u / sqrt(sum(grad_mu_u * grad_mu_u))
    grad_mu_v <- grad_mu_v / sqrt(sum(grad_mu_v * grad_mu_v))
    
    vp$u_mu <- vp$u_mu + learning_rate[epo] * grad_mu_u
    vp$v_mu <- vp$v_mu + learning_rate[epo] * grad_mu_v
    
    # print(vp$u_mu)
    # vp$u_sigma <- vp$u_sigma + learning_rate * grad_sigma_ui
    
    # Ut <- mvrnorm(n = 1,vp$u_mu[1,],diag(vp$u_sigma[1,]))
    # Vt <- mvrnorm(n = 1,vp$v_mu[1,],diag(vp$v_sigma[1,]))
    # for (i in 2:model$D) {
    #   Ut <- rbind(Ut,mvrnorm(n = 1,vp$u_mu[i,],diag(vp$u_sigma[i,])))
    #   Vt <- rbind(Vt,mvrnorm(n = 1,vp$v_mu[i,],diag(vp$v_sigma[i,])))
    # }
    Ut <- vp$u_mu
    Vt <- vp$v_mu
    model$U[[t]] <- Ut
    model$V[[t]] <- Vt
    # debug(likely)
    logL <- likely(X_t,model$U[[t]],model$V[[t]],model$D)
    print(logL)
    logLout <- c(logLout,logL)
    plot(Ut[,1:2],col = color)
  } # epo
  
  model$timePosition <- model$timePosition + 1 
  # return(model,logLout)
  return(logLout)
}

draw_VP <- function(model){
  vp <- list()
  vp$u_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  vp$v_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  vp$u_sigma <- abs(mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K))))
  vp$v_sigma <- abs(mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K))))
  return(vp)
}

likely <- function(X_t,U,V,D){
  ll <- 0
  for (i in 1:model$D) {
    for (j in 1:model$D) {
      u <- U[i,]
      v <- V[j,]
      x <- t(u) %*% v /(sqrt(sum(u * u)) * sqrt(sum(v * v)))
      
      if(x < -13){
        sig <- x 
      }else{
        if(x > 36){
          sig <- 0
        }else{
          sig <- log(1/(1 + exp(-x)))
        }      
      }
      gauche <- X_t$P[i,j] * sig
      
      # if(x > 13){
      #   sig <- -x 
      # }else{
      #   if(x < -36){
      #     sig <- 0
      #   }else{
      #     sig <- log(1/(1 + exp(x)))
      #   }      
      # }
      # droite <- X_t$N[i,j] * sig
      # ll <- ll + gauche + droite
      
      
      ll <- ll + gauche 
    
    }
  }
  return(ll)
}

logprobaQ <- function(U,V,vp,D){
  ll <- 0
  for (i in 1:D) {
    u <- U[i,]
    v <- V[i,]
    
    mu_u <- vp$u_mu[i,]
    mu_v <- vp$v_mu[i,]
    
    sigma_u <- diag(vp$u_sigma[i,])
    sigma_v <- diag(vp$v_sigma[i,])
    
    gauche <- log(dmvnorm(u,mu_u,sigma_u))
    droite <- log(dmvnorm(v,mu_v,sigma_v))
    
    ll <- ll + gauche + droite
  }
  if (is.infinite(ll)) {
    ll <- 0
  }
  return(ll)
}

probaCond <- function(U,V,model){
  ll <- 0
  for (i in 1:model$D) {
    u <- U[i,]
    v <- V[i,]
    
    sigma_u <- solve(solve(diag(model$vp_m1$u_sigma[i,]) + diag(rep(model$sigma_t,model$K))) + (diag(rep(1/model$sigma_0),model$K)))
    sigma_v <- solve(solve(diag(model$vp_m1$v_sigma[i,]) + diag(rep(model$sigma_t,model$K))) + (diag(rep(1/model$sigma_0,model$K))))
    mu_u <- (sigma_u * solve(diag(model$vp_m1$u_sigma[i,]))) %*% model$vp_m1$u_mu[i,]
    mu_v <- (sigma_u * solve(diag(model$vp_m1$v_sigma[i,]))) %*% model$vp_m1$v_mu[i,]
    
    gauche <- dmvnorm(u,mu_u,sigma_u)
    droite <- dmvnorm(v,mu_v,sigma_v)
    
    ll <- ll * gauche * droite
  }
  if (is.infinite(ll)) {
    ll <- 0
  }
  return(ll)
}

LprobaCond <- function(U,V,model){
  ll <- 0
  for (i in 1:model$D) {
    u <- U[i,]
    v <- V[i,]
    
    sigma_u <- solve(solve(diag(model$vp_m1$u_sigma[i,]) + diag(rep(model$sigma_t,model$K))) + (diag(rep(1/model$sigma_0),model$K)))
    sigma_v <- solve(solve(diag(model$vp_m1$v_sigma[i,]) + diag(rep(model$sigma_t,model$K))) + (diag(rep(1/model$sigma_0,model$K))))
    mu_u <- (sigma_u * solve(diag(model$vp_m1$u_sigma[i,]))) %*% model$vp_m1$u_mu[i,]
    mu_v <- (sigma_u * solve(diag(model$vp_m1$v_sigma[i,]))) %*% model$vp_m1$v_mu[i,]
    
    gauche <- dmvnorm(u,mu_u,sigma_u)
    droite <- dmvnorm(v,mu_v,sigma_v)
    
    ll <- ll + log(gauche) + log(droite)
  }
  if (is.infinite(ll)) {
    ll <- 0
  }
  return(ll)
}

source("model.R")
source("utils.R")
# (D,K,sigma_t,sigma_0,nb_epochs,nb_MB,nb_sampleVI)
model <- init_model(10,5,1,1,20,3,500)
model$timePosition
X_t <- list()

# GenerateData
n <- 10
K = 3
alpha = rep(1/K, K)
intra <- 0.9
inter <- (1 - intra)/(K-1)
pi = matrix(0.8, K, K)
diag(pi) = intra
e <- 1000
rSBM <- randomSBM(n,e,K,alpha,pi)

rSBM$Adj
rSBM$cluster


X_t$P <- as.matrix(rSBM$Adj)
X_t$P

X_t$N <- 5 * round((rowSums(X_t$P) %*% t(colSums(X_t$P))) / sum(X_t$P))
X_t$N

debug(optim)

blob <- optim(X_t,model,rSBM$cluster)
plot(blob)

