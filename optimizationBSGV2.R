
library(mixtools)

optim <- function(X_t,model,color){
  logLout <- c()
  learning_rate <- seq(0.1,0.01,length.out = model$nb_epochs)
  # print(learning_rate)
  
  vp <- draw_VP_withoutSigma(model) 
  # vp <- draw_VP(model) 
  U <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  V <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  
  for (epo in 1:model$nb_epochs) {
    print(paste("eopch :",epo))
    # epo <- 1
    grad_mu_u <- matrix(0,model$D,model$K)
    grad_mu_v <- matrix(0,model$D,model$K)
    
    # grad_sigma_u <- matrix(0,model$D,model$K)
    # grad_sigma_v <- matrix(0,model$D,model$K)
    
    
    for (s in 1:model$nb_sampleVI) {
      
      # print(paste("sample :",s))
      Esu <- mvrnorm(n = model$D,rep(0,model$K),diag(rep(1,model$K)))
      Esv <- mvrnorm(n = model$D,rep(0,model$K),diag(rep(1,model$K)))
      
      
      Us <- matrix(0,model$D,model$K)
      Vs <- matrix(0,model$D,model$K)
      
      Us <- vp$u_mu + Esu * vp$u_sigma
      Vs <- vp$v_mu + Esv * vp$v_sigma
      
      for (i in 1:model$D) {    
        u <- Us[i,]
        temp_grad <- 0
        for (j in 1:model$D) {    
          
          v <- V[j,]
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
        grad_mu_u[i,] <-  grad_mu_u[i,] +  (Us[i,] * temp_grad)
        grad_mu_u[i,] <- grad_mu_u[i,] / model$nb_sampleVI - (vp$u_mu[i,])
        # grad_mu_v <- grad_mu_v / model$nb_sampleVI - (vp$u_mu)
        
        vp$u_mu[i,] <- vp$u_mu[i,] + learning_rate[epo] * grad_mu_u[i,]
        # grad_mu_v[i,] <-  grad_mu_v[i,] +  (Vs[i,] * Esu * vp$v_sigma) * temp_grad
        
      }
      # 
      # grad_mu_u <- grad_mu_u / model$nb_sampleVI - (vp$u_mu)
      # # grad_mu_v <- grad_mu_v / model$nb_sampleVI - (vp$u_mu)
      # 
      # vp$u_mu <- vp$u_mu + learning_rate[epo] * grad_mu_u
      
      for (j in 1:model$D) {    
        v <- Vs[j,]
        temp_grad <- 0
        for (i in 1:model$D) {    
          
          u <- U[i,]
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
        grad_mu_v[j,] <-  grad_mu_v[j,] +  (Vs[j,] * temp_grad)
        grad_mu_v[j,] <- grad_mu_v[j,] / model$nb_sampleVI - (vp$v_mu[j,])
        # grad_mu_v <- grad_mu_v / model$nb_sampleVI - (vp$u_mu)
        
        vp$v_mu[j,] <- vp$v_mu[j,] + learning_rate[epo] * grad_mu_v[j,]
        # grad_mu_v[i,] <-  grad_mu_v[i,] +  (Vs[i,] * Esu * vp$v_sigma) * temp_grad
        
      }
    }
    # grad_mu_v <- grad_mu_v / model$nb_sampleVI - (vp$v_mu)
    # # grad_mu_v <- grad_mu_v / model$nb_sampleVI - (vp$u_mu)
    # 
    # vp$v_mu <- vp$v_mu + learning_rate[epo] * grad_mu_v
    # vp$v_mu <- vp$v_mu + learning_rate[epo] * grad_mu_v
    
    # grad_sigma_u <- grad_sigma_u / model$nb_sampleVI
    # grad_sigma_v <- grad_sigma_v / model$nb_sampleVI
    # vp$u_sigma <- vp$u_sigma + learning_rateS[epo] * grad_sigma_u
    # vp$v_sigma <- vp$v_sigma + learning_rateS[epo] * grad_sigma_v
    
    Ut <- vp$u_mu
    Vt <- vp$v_mu
    # model$U[[t]] <- Ut
    # model$V[[t]] <- Vt
    
    # debug(likely)
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

logprobaQI <- function(U,V,vp,D,I){
  ll <- 0
  i <- I
  u <- U[i,]
  v <- V[i,]
  
  mu_u <- vp$u_mu[i,]
  mu_v <- vp$v_mu[i,]
  
  sigma_u <- diag(vp$u_sigma[i,])
  sigma_v <- diag(vp$v_sigma[i,])
  
  gauche <- log(dmvnorm(u,mu_u,sigma_u))
  droite <- log(dmvnorm(v,mu_v,sigma_v))
  
  ll <- ll + gauche + droite
  if (is.infinite(ll)) {
    ll <- 0
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

source("model.R")
source("utils.R")
# (D,K,sigma_t,sigma_0,nb_epochs,nb_MB,nb_sampleVI)
model <- init_model(100,50,1,1,10,3,200)
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

# rSBM$Adj
# rSBM$cluster


X_t$P <- as.matrix(rSBM$Adj)
sum(X_t$N)
# X_t$P
X_t$N <- 10000 * round((rowSums(X_t$P) %*% t(colSums(X_t$P))) / sum(X_t$P))
# X_t$N

debug(optim)

blob <- optim(X_t,model,rSBM$cluster)
plot(blob$embedding,col = rSBM$cluster)
plot(blob$ll,type= "b")


