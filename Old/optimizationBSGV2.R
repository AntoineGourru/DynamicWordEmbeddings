
library(mixtools)

optim <- function(X_t,model,color){
  logLout <- c()
  learning_rate <- seq(0.00005,0.00001,length.out = model$nb_epochs)
  # print(learning_rate)
  
  vp <- draw_VP_withoutSigma(model) 
  # vp <- draw_VP(model) 
  Ut <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  Vt <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  
  for (epo in 1:model$nb_epochs) {
    print(paste("eopch :",epo))
    # epo <- 1
    grad_mu_u <- matrix(0,model$D,model$K)
    grad_mu_v <- matrix(0,model$D,model$K)
    
    # grad_sigma_u <- matrix(0,model$D,model$K)
    # grad_sigma_v <- matrix(0,model$D,model$K)
    
    for (i in 1:model$D) {  
      
      for (s in 1:model$nb_sampleVI) {
        
        # print(paste("sample :",s))
        Esu <- mvrnorm(n = 1,rep(0,model$K),diag(rep(1,model$K)))
        u <- vp$u_mu[i,] + Esu * vp$u_sigma[i,]
        
        temp_grad <- 0
        for (j in 1:model$D) {    
          
          v <- vp$v_mu[j,]
          # v <- Vt[j,]
          x <- t(u) %*% v
          # x <- t(u) %*% v /(sqrt(sum(u * u)) * sqrt(sum(v * v)))
          # print(x)
          if(x < -13) {
            sig <- x 
          }else{
            if(x > 36){
              sig <- 0
            }else{
              sig <- log(1/(1 + exp(-x)))
            }      
          }
          gauche <- X_t$N[i,j] * sig
          
          if(x > 13){
            sig <- -x
          }else{
            if(x < -36){
              sig <- 0
            }else{
              sig <- log(1/(1 + exp(x)))
            }
          }
          droite <- X_t$P[i,j] * sig
          
          temp_grad <- temp_grad + gauche + droite
          
        }
        grad_mu_u[i,] <-  grad_mu_u[i,] +  (u * temp_grad)
      }
      grad_mu_u[i,] <- grad_mu_u[i,] / model$nb_sampleVI - (vp$u_mu[i,])
      
      vp$u_mu[i,] <- vp$u_mu[i,] + learning_rate[epo] * grad_mu_u[i,]
    }
    
    for (j in 1:model$D) {    
      
      for (s in 1:model$nb_sampleVI) {
        
        # print(paste("sample :",s))
        Esu <- mvrnorm(n = 1,rep(0,model$K),diag(rep(1,model$K)))
        v <- vp$v_mu[j,] + Esu * vp$v_sigma[j,]
        
        temp_grad <- 0
        for (i in 1:model$D) {    
          
          u <- vp$u_mu[i,]
          # u <- Ut[i,]
          x <- t(u) %*% v
          # x <- t(u) %*% v /(sqrt(sum(u * u)) * sqrt(sum(v * v)))
          
          if(x < -13) {
            sig <- x 
          }else{
            if(x > 36){
              sig <- 0
            }else{
              sig <- log(1/(1 + exp(-x)))
            }      
          }
          gauche <- X_t$N[i,j] * sig
          
          if(x > 13){
            sig <- -x
          }else{
            if(x < -36){
              sig <- 0
            }else{
              sig <- log(1/(1 + exp(x)))
            }
          }
          droite <- X_t$P[i,j] * sig
          
          temp_grad <- temp_grad + gauche + droite
        }
        grad_mu_v[i,] <-  grad_mu_v[i,] +  (v * temp_grad)
      }
      grad_mu_v[i,] <- grad_mu_v[i,] / model$nb_sampleVI - (vp$v_mu[i,])
      
      vp$v_mu[i,] <- vp$v_mu[i,] + learning_rate[epo] * grad_mu_v[i,]
      
    }
    Ut <- vp$u_mu
    Vt <- vp$v_mu
    logL <- likely(X_t,Ut,Vt,model$D)
    print(logL)
    logLout <- c(logLout,logL)
    plot(Ut,col = color)
  } # epo
  
  model$timePosition <- model$timePosition + 1 
  # return(model,logLout)
  return(list(embedding = Ut,ll=logLout))
}

draw_VP_withoutSigma <- function(model){
  vp <- list()
  vp$u_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  vp$v_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  vp$u_sigma <- matrix(1,model$D,model$K)
  vp$v_sigma <- matrix(1,model$D,model$K)
  return(vp)
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
      x <- t(u) %*% v
      # x <- t(u) %*% v /(sqrt(sum(u * u)) * sqrt(sum(v * v)))
      
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
      
      if(x > 13){
        sig <- -x
      }else{
        if(x < -36){
          sig <- 0
        }else{
          sig <- log(1/(1 + exp(x)))
        }
      }
      droite <- X_t$N[i,j] * sig
      # droite <- 0
      ll <- ll + gauche + droite
      
      
      ll <- ll + gauche 
      
    }
    ll <- ll +  log(dmvnorm(U[i,],rep(0,model$K),diag(rep(1,model$K))))
    ll <- ll +  log(dmvnorm(V[i,],rep(0,model$K),diag(rep(1,model$K))))
  }
  return(ll)
}

source("model.R")
source("utils.R")
# (D,K,sigma_t,sigma_0,nb_epochs,nb_MB,nb_sampleVI)
model <- init_model(100,2,1,1,10,3,700)
model$timePosition
X_t <- list()

# GenerateData
n <- 100
K = 3
alpha = rep(1/K, K)
intra <- 0.9
inter <- (1 - intra)/(K-1)
pi = matrix(inter, K, K)
diag(pi) = intra
e <- 50000
rSBM <- randomSBM(n,e,K,alpha,pi)

# rSBM$Adj
# rSBM$cluster


X_t$P <- as.matrix(rSBM$Adj)
sum(X_t$P)
# X_t$P
X_t$N <- 5 * round((rowSums(X_t$P) %*% t(colSums(X_t$P))) / sum(X_t$P))
sum(X_t$N)
# X_t$N
e
# debug(optim)

blob <- optim(X_t,model,rSBM$cluster)
plot(blob$embedding,col = rSBM$cluster)
plot(blob$ll,type= "b")


