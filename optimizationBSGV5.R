
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
      denom <- 0
      for (samp in 1:model$nb_sampleVI) {
        
        # Esv <- mvrnorm(n = model$D,rep(0,model$K),diag(rep(1,model$K)))
        
        # U_sample <- do.call(rbind,(lapply(1:model$D,function(i) mvrnorm(n = 1,vp$u_mu[i,],diag(rep(1,model$K))))))
        U_sample <- Ut
        U_sample[i,] <- mvrnorm(n = 1,vp$u_mu[i,],diag(rep(1,model$K)))
        L_var <- likelyVar(X_t,U_sample,Ut,model$D)
        denom <- denom + L_var
        
        grad_mu_u[i,] <-  grad_mu_u[i,] + (U_sample[i,] * L_var)
      }
      grad_mu_u[i,] <- grad_mu_u[i,] / model$nb_sampleVI
      denom <- denom / model$nb_sampleVI - 1
      vp$u_mu[i,] <- grad_mu_u[i,] / rep(denom,2)
    }
    # grad_mu_u <- grad_mu_u / model$nb_sampleVI
    # denom <- denom / model$nb_sampleVI - 1
    # vp$u_mu <- grad_mu_u / rep(denom,2)
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
      
    }
    ll <- ll +  log(dmvnorm(U[i,],rep(0,model$K),diag(rep(1,model$K))))
    ll <- ll +  log(dmvnorm(V[i,],rep(0,model$K),diag(rep(1,model$K))))
  }
  return(ll)
}

likelyVar <- function(X_t,U,V,D){
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
    }
  }
  return(ll)
}

source("model.R")
source("utils.R")
# (D,K,sigma_t,sigma_0,nb_epochs,nb_MB,nb_sampleVI)
model <- init_model(100,2,1,1,10,3,10)
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
# debug(likelyVar)

blob <- optim(X_t,model,rSBM$cluster)
plot(blob$embedding,col = rSBM$cluster)
plot(blob$ll,type= "b")


