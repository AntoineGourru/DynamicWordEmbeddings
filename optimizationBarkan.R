optimi <- function(X_t,model,color){
  logLout <- c()
  
  vp <- draw_VP(model) 
  
  Ut <- vp$u_mu
  Vt <- vp$v_mu
  
  plot(Ut,col = color)
  logL <- likely(X_t,Ut,Vt)
  print(logL)
  
  C <- X_t$P - X_t$N 
  
  
  A <- vp$u_sigma * vp$u_sigma + vp$u_mu * vp$u_mu
  
  for (epo in 1:model$nb_epochs) {
    print(paste("eopch :",epo))
    
    #XI
    B <- vp$v_sigma*vp$v_sigma + vp$v_mu*vp$v_mu
    
    vp$xi <- sqrt(A %*% t(B))
    
    
    u_mu_temp <-  vp$u_mu
    v_mu_temp <- vp$v_mu
    u_sigma_temp <- vp$u_sigma
    v_sigma_temp <- vp$v_sigma
    
    for (i in 1:model$D) {
      # R
      r_ui <- 0.5 * colSums(C[i,] * vp$v_mu)
      
      #P
      P <- matrix(0,model$K,model$K)
      for (j in 1:model$D) {
        E <- diag(vp$v_sigma[j,]) + vp$v_mu[j,] %*% t(vp$v_mu[j,])
        sig_xi <- 1/(1 + exp(-vp$xi[i,j]))
        deux_lambda_xi <- (1/sig_xi)  * (sig_xi  - 0.5)
        
        P <- P + (deux_lambda_xi *  E)
        
      }
      P <- P + diag(rep(1,model$K))
      
      Pm1 <- solve(P)
      u_mu_temp[i,] <- Pm1 %*% r_ui
      
      u_sigma_temp[i,] <- diag(Pm1)
      
    }
    
    vp$u_mu <- u_mu_temp 
    vp$u_sigma  <- u_sigma_temp
    
    #XI
    
    A <- vp$u_sigma * vp$u_sigma + vp$u_mu * vp$u_mu
    
    vp$xi <- sqrt(A %*% t(B))
    
    for (i in 1:model$D) {
      # R
      r_vi <- 0.5 * colSums(C[,i] * vp$u_mu)
      
      #P
      P <- matrix(0,model$K,model$K)
      for (j in 1:model$D) {
        E <- diag(vp$u_sigma[j,]) + vp$u_mu[j,] %*% t(vp$u_mu[j,])
        sig_xi <- 1/(1 + exp(-vp$xi[j,i]))
        deux_lambda_xi <- (1/sig_xi)  * (sig_xi  - 0.5)
        
        P <- P + (deux_lambda_xi *  E)
        
      }
      P <- P + diag(rep(1,model$K))
      
      Pm1 <- solve(P)
      v_mu_temp[i,] <- Pm1 %*% r_vi
      
      v_sigma_temp[i,] <- diag(Pm1)
      
    }
    
    
    vp$v_mu <- v_mu_temp
    vp$v_sigma <- v_sigma_temp
    
    Ut <- vp$u_mu
    Vt <- vp$v_mu
    
    logL <- likely(X_t,Ut,Vt)
    print(logL)
    logLout <- c(logLout,logL)
    plot(Ut,col = color)
  } # epo
  
  model$timePosition <- model$timePosition + 1 
  return(list(embedding = Ut,ll=logLout,sigma = vp$u_sigma))
}

draw_VP <- function(model){
  vp <- list()
  vp$u_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  vp$v_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1,model$K)))
  vp$u_sigma <- matrix(1,model$D,model$K)
  vp$v_sigma <- matrix(1,model$D,model$K)
  vp$xi <- matrix(0,model$D,model$D)
  return(vp)
}

likely <- function(X_t,U,V){
  require(mixtools)
  x <- U %*% t(V)
  
  sig <- -log(1 + exp(-x))
  
  gauche <- X_t$P * sig
  
  sig <- -log(1 + exp(x))
  
  droite <- X_t$N * sig
  
  ll <- sum(gauche + droite)
  
  ll <- ll +  sum(apply(U,1,function(x){log(dmvnorm(x,rep(0,model$K),diag(rep(1,model$K))))}))
  ll <- ll +  sum(apply(V,1,function(x){log(dmvnorm(x,rep(0,model$K),diag(rep(1,model$K))))}))
  return(ll)
}



source("model.R")
source("utils.R")
# (D,K,sigma_t,sigma_0,nb_epochs,nb_MB,nb_sampleVI)
n <- 300

model <- init_model(n,2,1,1,nb_epochs = 100,3,500)
model$timePosition
X_t <- list()

# GenerateData

K = 3
alpha = rep(1/K, K)
intra <- 0.8
inter <- (1 - intra)/(K-1)
pi = matrix(inter, K, K)
diag(pi) = intra
e <- n * 100
rSBM <- randomSBM(n,e,K,alpha,pi)

# rSBM$Adj
# rSBM$cluster


X_t$P <- as.matrix(rSBM$Adj)
X_t$N <- 4 * round((rowSums(X_t$P) %*% t(colSums(X_t$P))) / sum(X_t$P))
# debug(optimi)
blob <- optimi(X_t,model,rSBM$cluster)
plot(blob$ll)




