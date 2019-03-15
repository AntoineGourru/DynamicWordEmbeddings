optim <- function(X_t,model,color){
  library(mixtools)
  logLout <- c()
  
  vp <- draw_VP_withoutSigma(model) 
  
  Ut <- vp$u_mu
  colnames(Ut) <- c("V1","V2")
  p <- plot_ly(as.data.frame(Ut), x = ~V1, y = ~V2, type = 'scatter', mode = 'markers',
               text = ~color,hoverinfo = 'text', marker=list(opacity=0.5) )
  p
  logL <- likely(X_t,Ut,Ut,model$D)
  print(logL)
  # Vt <- vp$v_mu
  print(X_t$P[1,2])
  for (epo in 1:model$nb_epochs) {
    print(paste("eopch :",epo))
    # epo <- 1
    grad_mu_u <- matrix(0,model$D,model$K)
    
    for (i in 1:model$D) {
      # i <- 1
      denom <- 0
      for (samp in 1:model$nb_sampleVI) {
        
        u <- mvrnorm(n = 1,vp$u_mu[i,],diag(rep(3,model$K)))
        
        ll <- 0
        for (j in setdiff(1:model$D,i)){
          v <- vp$u_mu[j,]
          x <- t(u) %*% v
          
          sig <- 1/(1 + exp(-x))
          
          gauche <- X_t$P[i,j] * sig
          
          sig <- 1/(1 + exp(x))
          
          droite <- X_t$N[i,j] * sig
          # droite <- 0
          ll <- ll + gauche + droite
          # ll <- ll + gauche 
        }
        
        denom <- denom + ll
        grad_mu_u[i,] <-  grad_mu_u[i,] + (u * ll)
      }
      
      
      grad_mu_u[i,] <- grad_mu_u[i,] / model$nb_sampleVI
      denom <- denom / model$nb_sampleVI - 1
      vp$u_mu[i,] <- grad_mu_u[i,] / rep(denom,2)
    }
    
    Ut <- vp$u_mu
    logL <- likely(X_t,Ut,Ut,model$D)
    print(logL)
    logLout <- c(logLout,logL)
    p <- plot_ly(as.data.frame(Ut), x = ~V1, y = ~V2, type = 'scatter', mode = 'markers',
                 text = ~color,hoverinfo = 'text', marker=list(opacity=0.5) )
    p
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
  print(X_t$P[1,2])
  for (i in 1:D) {
    u <- U[i,]
    for (j in setdiff(1:D,i)){
      v <- V[j,]
      x <- t(u) %*% v
      
      sig <- 1/(1 + exp(-x))
      
      gauche <- X_t$P[i,j] * sig
      
      sig <- 1/(1 + exp(x))
      
      droite <- X_t$N[i,j] * sig
      # droite <- 0
      ll <- ll + gauche + droite
      # ll <- ll + gauche 
    }
    ll <- ll +  log(dmvnorm(U[i,],rep(0,model$K),diag(rep(1,model$K))))
    ll <- ll +  log(dmvnorm(V[i,],rep(0,model$K),diag(rep(1,model$K))))
  }
  return(ll)
}

source("model.R")
source("utils.R")
# (D,K,sigma_t,sigma_0,nb_epochs,nb_MB,nb_sampleVI)
n <- 136

model <- init_model(n,2,1,1,100,3,100)
model$timePosition
X_t <- list()

StructuralHashtags <- read.csv("input/StructuralHashtags.csv", header=FALSE, stringsAsFactors=FALSE)
nodesHashtag <- read.table("input/nodesHashtag.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)

X_t$P <- as.matrix(StructuralHashtags)

# X_t$P
X_t$N <- 5 * round((rowSums(X_t$P) %*% t(colSums(X_t$P))) / sum(X_t$P))
sum(X_t$N)
library("plotly")
blob <- optim(X_t,model,nodesHashtag)
plot(blob$ll)
p <- plot_ly(as.data.frame(blob$embedding), x = ~V1, y = ~V2, type = 'scatter', mode = 'markers',
             text = ~nodesHashtag,hoverinfo = 'text', marker=list(opacity=0.5) )
p

