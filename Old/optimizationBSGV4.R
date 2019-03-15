
library(mixtools)

optim <- function(X_t,model,color){
  logLout <- c()
  # print(learning_rate)
  
  vp <- draw_VP_withoutSigma(model) 

  Ut <- vp$u_mu
  plot(Ut,col = color)
  logL <- likely(X_t,Ut,Ut,model$D)
  print(logL)
  # Vt <- vp$v_mu
  
  for (epo in 1:model$nb_epochs) {
    print(paste("eopch :",epo))
    # epo <- 1
    grad_mu_u <- matrix(0,model$D,model$K)
    # grad_mu_v <- matrix(0,model$D,model$K)
    
    # grad_sigma_u <- matrix(0,model$D,model$K)
    # grad_sigma_v <- matrix(0,model$D,model$K)
    
    for (i in 1:model$D) {  
      
      for(samp in 1:model$nb_sampleVI){
        
        # Esv <- mvrnorm(n = model$D,rep(0,model$K),diag(rep(1,model$K)))
        
        u <- vp$u_mu[i,] + mvrnorm(n = 1,rep(0,model$K),diag(rep(1,model$K))) * vp$u_sigma[i,]

        temp_grad <- 0
        
        for (j in 1:model$D) {
          
          # v <- vp$u_mu[j,] + Esv[j,] * vp$u_sigma[j,]
          v <- vp$u_mu[j,] 
          
          x <- t(u) %*% v
          
          sig <- 1/(1 + exp(x))
          
          gauche <- X_t$P[i,j] * sig
          
          sig <- 1/(1 + exp(-x))

          droite <- X_t$N[i,j] * sig
          
          # grad_mu_u[i,] <-  grad_mu_u[i,] + (v * (gauche - droite))
          temp_grad <- temp_grad + gauche + droite
          grad_mu_u[i,] <-  grad_mu_u[i,] + (u * (gauche + droite))
        }
      }
      # vp$u_mu[i,] <- grad_mu_u[i,] / model$nb_sampleVI
      vp$u_mu[i,] <- grad_mu_u[i,] / temp_grad - 1
      
    }
    
    Ut <- vp$u_mu
    logL <- likely(X_t,Ut,Ut,model$D)
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
  vp$u_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(0.1,model$K)))
  vp$v_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(0.1,model$K)))
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
      
      if (x < -13) {
        sig <- x 
      }else{
        if (x > 36) {
          sig <- 0
        }else{
          sig <- log(1/(1 + exp(-x)))
        }      
      }
      gauche <- X_t$P[i,j] * sig
      
      if (x > 13) {
        sig <- -x
      }else{
        if (x < -36) {
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
model <- init_model(100,2,1,1,10,3,100)
model$timePosition
X_t <- list()

# GenerateData
n <- 100
K = 3
alpha = rep(1/K, K)
intra <- 0.8
inter <- (1 - intra)/(K-1)
pi = matrix(inter, K, K)
diag(pi) = intra
e <- 50000
rSBM <- randomSBM(n,e,K,alpha,pi)

rSBM$Adj
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


