optimZero <- function(X_t,model,vp){
  logLout <- c()

  Ut <- vp$u_mu
  Vt <- vp$v_mu
  
  plot(Ut)
  logL <- likely(X_t,Ut,Vt)
  print(logL)
  
  C <- X_t$P - X_t$N 

  A <- vp$u_sigma * vp$u_sigma + vp$u_mu * vp$u_mu
  
  for (epo in 1:model$nb_epochs) {
  # for (epo in 1:) {
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
    plot(Ut)
  } # epo
  
  model$U[[model$timePosition]] <- Ut
  model$V[[model$timePosition]] <- Vt
  model$LL[[model$timePosition]] <- logLout

  model$timePosition <- model$timePosition + 1 
  
  return(list(model = model,vp = vp))
}

optimT <- function(X_t,model,vp){
  logLout <- c()
  
  Ut <- vp$u_mu
  Vt <- vp$v_mu
  
  # plot(Ut,col = color)
  logL <- likely(X_t,Ut,Vt)
  print(logL)
  
  C <- X_t$P - X_t$N 
  
  Vm1 <- vp$v_mu
  Um1 <- vp$u_mu
  sigma_U_Pm1 <- vp$u_sigma
  sigma_V_Pm1 <- vp$v_sigma
  
  for(i in 1:model$D){
    sigma_U_Pm1[i,] <- (1/(vp$u_sigma[i,] + model$sigma_t) + 1/model$tau)
    Um1[i,] <- sigma_U_Pm1[i,] * (1/(vp$u_sigma[i,] + model$sigma_t)) * vp$u_mu[i,]
  }
  
  for(i in 1:model$D){
    sigma_V_Pm1[i,] <- (1/(vp$v_sigma[i,] + model$sigma_t) + 1/model$tau)
    Vm1[i,] <- sigma_V_Pm1[i,] * (1/(vp$v_sigma[i,] + model$sigma_t)) * vp$v_mu[i,]
  }
  
  
  
  A <- vp$u_sigma * vp$u_sigma + vp$u_mu * vp$u_mu
  for (epo in 1:model$nb_epochs) {
    print(paste("eopch :",epo))
    
    #XI
    B <- vp$v_sigma*vp$v_sigma + vp$v_mu*vp$v_mu
    for (i in 1:model$D) {
      
      bob <- A[i,] * B
      
      vp$xi[i,] <- sqrt(rowSums(bob))
      
    }
    
    u_mu_temp <-  vp$u_mu
    v_mu_temp <- vp$v_mu
    u_sigma_temp <- vp$u_sigma
    v_sigma_temp <- vp$v_sigma
    
    for (i in 1:model$D) {
      # R
      r_ui <- 0.5 * colSums(C[i,] * vp$v_mu) + (sigma_U_Pm1[i,] * Um1[i,])
      
      #P
      P <- matrix(0,model$K,model$K)
      for (j in 1:model$D) {
        E <- diag(vp$v_sigma[j,]) + vp$v_mu[j,] %*% t(vp$v_mu[j,])
        sig_xi <- 1/(1 + exp(-vp$xi[i,j]))
        deux_lambda_xi <- (1/sig_xi)  * (sig_xi  - 0.5)
        
        P <- P + (deux_lambda_xi *  E)
        
      }
      P <- P + diag(rep(1,model$K)) + sigma_U_Pm1[i,]
      
      Pm1 <- solve(P)
      u_mu_temp[i,] <- Pm1 %*% r_ui
      
      u_sigma_temp[i,] <- diag(Pm1)
      
    }
    
    vp$u_mu <- u_mu_temp 
    vp$u_sigma  <- u_sigma_temp
    
    #XI
    
    A <- vp$u_sigma * vp$u_sigma + vp$u_mu * vp$u_mu
    for (i in 1:model$D) {
      
      bob <- A[i,] * B
      
      vp$xi[i,] <- sqrt(rowSums(bob))
      # print(vp$xi[i,j])
      
    }
    
    for (i in 1:model$D) {
      # R
      r_vi <- 0.5 * colSums(C[,i] * vp$u_mu) + (sigma_V_Pm1[i,] * Vm1[i,])
      
      #P
      P <- matrix(0,model$K,model$K)
      for (j in 1:model$D) {
        E <- diag(vp$u_sigma[j,]) + vp$u_mu[j,] %*% t(vp$u_mu[j,])
        sig_xi <- 1/(1 + exp(-vp$xi[j,i]))
        deux_lambda_xi <- (1/sig_xi)  * (sig_xi  - 0.5)
        
        P <- P + (deux_lambda_xi *  E)
        
      }
      P <- P + diag(rep(1,model$K)) + sigma_V_Pm1[i,]
      
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
  } # epo
  
  model$U[[model$timePosition]] <- Ut
  model$V[[model$timePosition]] <- Vt
  model$LL[[model$timePosition]] <- logLout
  
  model$timePosition <- model$timePosition + 1 
  
  return(list(model = model,vp = vp))
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


init_model <- function(D,K,sigma_t,tau,nb_epochs){
  require(MASS)
  
  model <- list()
  model$sigma_t <- sigma_t
  # model$sigma_0 <- sigma_0
  model$tau <- tau
  
  model$K <- K
  model$D <- D
  
  model$nb_epochs <- nb_epochs
  
  model$U <- list()
  model$V <- list()
  
  model$LL <- list()
  
  model$timePosition <- 1
  return(model)
}

# Main

library(readr)
data <- read_delim("export_articles_EGC_2004_2018.csv", 
                                            "\t", escape_double = FALSE, na = "empty", 
                                            trim_ws = TRUE)

data <- data[,c(3,4,5)]
data[,4] <- apply(data,1,function(x){paste(x[2],x[3])})
data <- data[,c(1,4)]

# Create X
X_t <- list()
library(text2vec)
# iterator <- itoken(data$V4[which(data$year == 2004)], tokenizer=space_tokenizer, progressbar=FALSE)
iterator <- itoken(data$V4, tokenizer=space_tokenizer, progressbar=FALSE)
vocabulary <- create_vocabulary(iterator)
print(sum(vocabulary$term_count))
pruned_vocabulary <- prune_vocabulary(vocabulary,term_count_min = 20)
vocab <- pruned_vocabulary$term
N <- length(vocab)

# vectorizer <- vocab_vectorizer(pruned_vocabulary)
vectorizer <- vocab_vectorizer(pruned_vocabulary)
l <- 5
X_t$P <- as.matrix(create_tcm(iterator, vectorizer, skip_grams_window=l,skip_grams_window_context = "symmetric", weights=rep(1, l)))
X_t$P <- X_t$P + t(X_t$P)
sum(X_t$P)
colSums(X_t$P)

freq <- colSums(X_t$P)
freq <- freq/sum(freq)
bob <- round(apply(X_t$P,1,function(x){x - freq*sum(x)}))
bob[bob<0] <- 0
sum(X_t$P)
sum(bob)
which(rowSums(bob) == 0)
which(colSums(bob) == 0)
X_t$P <- bob
X_t$N <- 5 * round((rowSums(X_t$P) %*% t(colSums(X_t$P))) / sum(X_t$P))
sum(X_t$P)
which(rowSums(X_t$P)==0)

# T 1
model <- init_model(length(vocab),10,1,1,nb_epochs = 20)
vp <- draw_VP(model) 
# plot(vp$u_mu)
# debug(optimZero)
out <- optimZero(X_t,model,vp)
model <- out$model
vp <- out$vp

emb <- model$U[[1]] / apply(model$U[[1]],1,function(x){sqrt(sum(x * x))})
# apply(emb,1,function(x){sqrt(sum(x * x))})
cosine <- emb %*% t(emb)
res <- c()
for(i in 1:N){
  cso <- cosine[i,]
  cso[i] <- -Inf
  res <- c(res,which.max(cso))
}
resultat <- cbind(1:N,vocab,res,vocab[res])
View(resultat)
# # T 2
# iterator <- itoken(data$V4[which(data$year == 2005)], tokenizer=space_tokenizer, progressbar=FALSE)
# l <- 5
# X_t$P <- as.matrix(create_tcm(iterator, vectorizer, skip_grams_window=l, weights=rep(1, l)))
# X_t$N <- 5 * round((rowSums(X_t$P) %*% t(colSums(X_t$P))) / sum(X_t$P))
# 
# # debug(optimT)
# out <- optimT(X_t,model,cl,vp)
# model <- out$model
# vp <- out$vp
# 
# plot(model$LL[[2]])
# 
# 
# 
