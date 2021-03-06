sigmoid <- function(x){
  return(1/(1 + exp(-x)))
}

optimZero <- function(X_t,model,vp){
  logLout <- c()
  
  Ut <- vp$u_mu
  Vt <- vp$v_mu
  
  # logL <- likely(X_t$P,X_t$N,Ut,Vt)
  # print(logL)
  
  freq <- colSums(X_t$P)
  freq <- freq /sum(freq)
  freq <-  1 - sqrt(0.0001/freq)
  freq[freq<0] <- 0
  print(freq)
  N <- length(vocab)
  
  #Init P
  oldPu <- list()
  oldPv <- list()
  oldRu <- list()
  oldRv <- list()
  beta <- 0
  for (epo in 1:model$nb_epochs) {
    print(paste("eopch :",epo))
    Pos <- X_t$P
    for (j in 1:N){
      for(i in 1:N){
        Pos[i,j] <- Pos[i,j] - sum(rbinom(Pos[i,j],1,freq[j]))
      }
    }
    
    freqNeg <- (rowSums(Pos)/sum(rowSums(Pos)))^(3/4)
    freqNeg <- freqNeg /sum(freqNeg)
    
    NN <- rowSums(Pos)
    Neg <- Pos
    for(i in 1:nrow(X_t$P)){
      Neg[i,] <- rmultinom(1,NN[i],freqNeg)
    }
    
    C <- Pos - Neg
    
    
    for (i in 1:model$D) {
      
      a_i = vp$u_sigma[i,] + vp$u_mu[i,] * vp$u_mu[i,]
      b_i = vp$v_sigma[i,] + vp$v_mu[i,] * vp$v_mu[i,]
      
      P_u <- matrix(0,model$K,model$K)
      R_u <- rep(0,model$K)
      
      P_v <- matrix(0,model$K,model$K)
      R_v <- rep(0,model$K)
      

      for(j in 1:model$D){
        
        
        a_j <- vp$u_sigma[j,] + vp$u_mu[j,] * vp$u_mu[j,]
        b_j <- vp$v_sigma[j,] + vp$v_mu[j,] * vp$v_mu[j,]
        
        xij <- sqrt(a_i %*% b_j)
        deux_lambda_xi <- (1/xij)  * (sigmoid(xij)  - 0.5)
        
        E <- diag(vp$v_sigma[j,]) + vp$v_mu[j,] %*% t(vp$v_mu[j,])
        
        P_u <- P_u + abs(C[i,j]) * (matrix(deux_lambda_xi,model$K,model$K) *  E)
        R_u <- R_u + 0.5 * C[i,j] * vp$v_mu[j,]
        
        
        xji <- sqrt(b_i %*% a_j)
        deux_lambda_xi <- (1/xji)  * (sigmoid(xji)  - 0.5)
        
        E <- diag(vp$u_sigma[j,]) + vp$u_mu[j,] %*% t(vp$u_mu[j,])
        
        P_v <- P_v + abs(C[j,i]) * (matrix(deux_lambda_xi,model$K,model$K) *  E)
        
        R_v <- R_v + 0.5 * C[j,i] * vp$u_mu[j,]
      }
      
      P_u <- P_u + diag(rep(model$tau,model$K))
      P_v <- P_v + diag(rep(model$tau,model$K))
      
      
      
      if(epo>5){
        beta <- (epo-5)^(-0.7)
        P_u <-  beta * P_u + (1-beta) * oldPu[[i]]
        P_v <-  beta * P_v + (1-beta) * oldPv[[i]]
        R_u <-  beta * R_u + (1-beta) * oldRu[[i]]
        R_v <-  beta * R_v + (1-beta) * oldRv[[i]]
      }
      # Update U
      Pm1 <- solve(P_u)
      vp$u_mu[i,] <- Pm1 %*% R_u
      vp$u_sigma[i,] <- diag(Pm1)
      
      # Update V
      Pm1 <- solve(P_v)
      vp$v_mu[i,] <- Pm1 %*% R_v
      vp$v_sigma[i,] <- diag(Pm1)
      
      
      oldPu[[i]] <- P_u
      oldPv[[i]] <- P_v
      oldRu[[i]] <- R_u
      oldRv[[i]] <- R_v
      
      
    }

    print(beta)
    Ut <- vp$u_mu
    Vt <- vp$v_mu
    
    logL <- likely(X_t$P,X_t$N,Ut,Vt)
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


draw_VP <- function(model){
  vp <- list()
  vp$u_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1/model$tau,model$K)))
  vp$v_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1/model$tau,model$K)))
  vp$u_sigma <- matrix(1/model$tau,model$D,model$K)
  vp$v_sigma <- matrix(1/model$tau,model$D,model$K)
  vp$xi <- matrix(0,model$D,model$D)
  return(vp)
}

likely <- function(Pos,Neg,U,V){
  require(mixtools)
  x <- U %*% t(V)
  
  sig <- -log(1 + exp(-x))
  
  gauche <- Pos * sig
  
  sig <- -log(1 + exp(x))
  
  droite <- Neg * sig
  
  ll <- sum(gauche + droite)
  
  # ll <- ll +  sum(apply(U,1,function(x){log(dmvnorm(x,rep(0,model$K),diag(rep(1/model$tau,model$K))))}))
  # ll <- ll +  sum(apply(V,1,function(x){log(dmvnorm(x,rep(0,model$K),diag(rep(1/model$tau,model$K))))}))
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
data <- read_delim("../../data/09_12.csv",
                   "\t", escape_double = FALSE, na = "empty",
                   trim_ws = TRUE)
# ind <- c(which(data$year == 2017),which(data$year == 2018))
ind <- 1:nrow(data)
data <- data[,c(3,4,5)]
data[,4] <- apply(data,1,function(x){paste(x[2],x[3])})
data <- data[,c(1,4)]

data <- data[ind,]
# Create X
X_t <- list()
library(tm)
tmCorpus <- Corpus(VectorSource(data$V4))
tmCorpus <- tm_map(tmCorpus, stripWhitespace)
tmCorpus <- tm_map(tmCorpus, content_transformer(tolower))
tmCorpus <- tm_map(tmCorpus, removeNumbers)
tmCorpus <- tm_map(tmCorpus, removePunctuation, preserve_intra_word_dashes = TRUE)
tmCorpus <- tm_map(tmCorpus, removeWords, stopwords("french"))
tmCorpus <- tm_map(tmCorpus, removeWords, stopwords("english"))
txt <- data.frame(text = get("content", tmCorpus),stringsAsFactors = FALSE)$text


###########Creating cooccurnece matrix 
library(text2vec)
iterator <- itoken(txt, tokenizer=space_tokenizer, progressbar=FALSE)
vocabulary <- create_vocabulary(iterator)
# print(sum(vocabulary$term_count))
pruned_vocabulary <- prune_vocabulary(vocabulary,  term_count_min = 20,doc_proportion_max = 0.8)
# pruned_vocabulary <- prune_vocabulary(vocabulary,  term_count_min = 10)
vectorizer <- vocab_vectorizer(pruned_vocabulary)
l <- 4
X_t$P <- as.matrix(create_tcm(iterator, vectorizer, skip_grams_window=l,skip_grams_window_context = "symmetric", weights=rep(1, l)))
X_t$P <- X_t$P + t(X_t$P)
vocab <- pruned_vocabulary$term
View(cbind(vocab,colSums(X_t$P)))

freqNeg <- (rowSums(X_t$P)/sum(rowSums(X_t$P)))^(3/4)
freqNeg <- freqNeg /sum(freqNeg)

NN <- rowSums(X_t$P)
X_t$N <- X_t$P
for(i in 1:nrow(X_t$P)){
  X_t$N[i,] <- rmultinom(1,NN[i],freqNeg)
}

sum(X_t$P)
sum(X_t$N)

# debug(optimZero)
model <- init_model(length(vocab),50,1,tau = 1,nb_epochs = 10)
vp <- draw_VP(model) 

out <- optimZero(X_t,model,vp)
model <- out$model
plot(model$LL[[1]])
vp <- out$vp

N <- length(vocab)
emb <- model$U[[1]] / apply(model$U[[1]],1,function(x){sqrt(sum(x * x))})
# emb <- model$U[[1]]
# embV <- model$V[[1]] / apply(model$V[[1]],1,function(x){sqrt(sum(x * x))})
# embV <- model$V[[1]]
# apply(emb,1,function(x){sqrt(sum(x * x))})
cosine <- emb %*% t(emb)
# cosine <- emb %*% t(embV)
res <- matrix("",N,4)
for(i in 1:N){
  cso <- cosine[i,]
  cso[i] <- -Inf
  res[i,1] <- vocab[i]
  res[i,2:4] <- vocab[sort(cso,decreasing = T,index.return=T)$ix[1:3]]
}
View(res)

# library(plotly)
# Sys.setenv("plotly_username"="AntoineGourru")
# Sys.setenv("plotly_api_key"="Vhnv9xPEeRYTicy0cenf")
# ##Ploty
# embedding <- as.data.frame(model$U[[1]] )
# 
# p <- plot_ly(embedding, x = ~embedding$V1, y = ~embedding$V2, type = 'scatter', mode = 'markers',
#              text = ~vocab,hoverinfo = 'text')
# p
# 
# api_create(p, filename = "Barkan")
# T 2
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



