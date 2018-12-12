
library(mixtools)

optim <- function(X_t,model,color){
  logLout <- c()
  learning_rate <- seq(0.001,0.0001,length.out = model$nb_epochs)
  learning_rateS <- seq(0.00000005,0.00000005,length.out = model$nb_epochs)
  # print(learning_rate)
  
  vp <- draw_VP_withoutSigma(model) 
  # vp <- draw_VP(model) 
  
  for (epo in 1:model$nb_epochs) {
    # epo <- 1
    grad_mu_u <- matrix(0,model$D,model$K)
    grad_mu_v <- matrix(0,model$D,model$K)
    grad_sigma_u <- matrix(0,model$D,model$K)
    grad_sigma_v <- matrix(0,model$D,model$K)
    
    
    for (s in 1:model$nb_sampleVI) {
      
      Us <- matrix(0,model$D,model$K)
      Vs <- matrix(0,model$D,model$K)
      for (i in 1:model$D) {
        Us[i,] <- mvrnorm(n = 1,vp$u_mu[i,],diag(vp$u_sigma[i,]))
        Vs[i,] <- mvrnorm(n = 1,vp$v_mu[i,],diag(vp$v_sigma[i,]))
      }
      # PQ <- probaQ(Us,Vs,vp,model$D)
      LL <- likely(X_t,Us,Vs,model$D)
      # lPC <- LprobaCond(Us,Vs,model)
      # lPQ <- logprobaQ(Us,Vs,vp,model$D)
      
      for (i in 1:model$D) {    
        
        #Learning, on peut mettre ici une condition d'arret
        
        #BLACKBOX VI pour U
        
        #To verify
        # temp_gradU <- (LL - lPQ)*(-(vp$u_m[i,] - Us[i,])/vp$u_sigma[i,])
        temp_gradU <- LL*(-(vp$u_m[i,] - Us[i,])/vp$u_sigma[i,])

        grad_mu_u[i,] <-  grad_mu_u[i,] + temp_gradU

        # temp_gradV <- (LL - lPQ)*(-(vp$v_m[i,] - Vs[i,])/vp$v_sigma[i,])
        temp_gradV <- LL*(-(vp$v_m[i,] - Vs[i,])/vp$v_sigma[i,])

        grad_mu_v[i,] <-  grad_mu_v[i,] + temp_gradV
        # 
        # #To verify for sigma
        # 
        # temp_gradU <- (LL - lPQ)*((vp$u_m[i,] - Us[i,])^2/(2 * vp$u_sigma[i,]^2))
        # 
        # grad_sigma_u[i,] <-  grad_sigma_u[i,] + temp_gradU
        # 
        # temp_gradV <- (LL - lPQ)*((vp$v_m[i,] - Vs[i,])^2/(2 * vp$v_sigma[i,]^2))
        # 
        # grad_sigma_v[i,] <-  grad_sigma_v[i,] + temp_gradV

      }
    }
    grad_mu_u <- grad_mu_u / model$nb_sampleVI
    grad_mu_v <- grad_mu_v / model$nb_sampleVI
    
    vp$u_mu <- vp$u_mu + learning_rate[epo] * grad_mu_u
    vp$v_mu <- vp$v_mu + learning_rate[epo] * grad_mu_v
    
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
model <- init_model(100,5,1,1,50,3,1000)
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
e <- 1000
rSBM <- randomSBM(n,e,K,alpha,pi)

# rSBM$Adj
# rSBM$cluster


X_t$P <- as.matrix(rSBM$Adj)
# X_t$P
X_t$N <- 5 * round((rowSums(X_t$P) %*% t(colSums(X_t$P))) / sum(X_t$P))
# X_t$N

# debug(optim)

blob <- optim(X_t,model,rSBM$cluster)
plot(blob$embedding,col = rSBM$cluster)
plot(blob$ll,type= "b")


