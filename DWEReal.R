optimZero <- function(X_t,model,vp){
  logLout <- c()

  Ut <- vp$u_mu
  Vt <- vp$v_mu
  
  plot(Ut)
  logL <- likely(X_t,Ut,Vt)
  print(logL)
  
  C <- X_t$P - X_t$N 
  
  for (epo in 1:model$nb_epochs) {
  # for (epo in 1:) {
    print(paste("eopch :",epo))
    
    for (ep2 in 1:10){
    #XI
    A <- vp$u_sigma * vp$u_sigma + vp$u_mu * vp$u_mu
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
        deux_lambda_xi <- (1/vp$xi[i,j])  * (sig_xi  - 0.5)
        
        P <- P + (deux_lambda_xi *  E)
        
      }
      P <- P + diag(rep(model$tau,model$K))

      Pm1 <- solve(P)
      u_mu_temp[i,] <- Pm1 %*% r_ui
      
      u_sigma_temp[i,] <- diag(Pm1)
      
    }
  
    vp$u_mu <- u_mu_temp
    vp$u_sigma  <- u_sigma_temp
    }

    #XI
    
    for (ep2 in 1:10){
    A <- vp$u_sigma * vp$u_sigma + vp$u_mu * vp$u_mu
    B <- vp$v_sigma*vp$v_sigma + vp$v_mu*vp$v_mu
    vp$xi <- sqrt(A %*% t(B))
    
    for (i in 1:model$D) {
      # R
      r_vi <- 0.5 * colSums(C[,i] * vp$u_mu)
      
      #P
      P <- matrix(0,model$K,model$K)
      for (j in 1:model$D) {
        E <- diag(vp$u_sigma[j,]) + vp$u_mu[j,] %*% t(vp$u_mu[j,])
        sig_xi <- 1/(1 + exp(-vp$xi[j,i]))
        deux_lambda_xi <- (1/vp$xi[j,i])  * (sig_xi  - 0.5)
        
        P <- P + (deux_lambda_xi *  E)
        
      }
      P <- P + diag(rep(model$tau,model$K))
      
      Pm1 <- solve(P)
      v_mu_temp[i,] <- Pm1 %*% r_vi
      
      v_sigma_temp[i,] <- diag(Pm1)
      
    }
    

    vp$v_mu <- v_mu_temp
    vp$v_sigma <- v_sigma_temp
    
    }
    
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

## Attention il faut refaire le copier coler prpre tahu
# optimT <- function(X_t,model,vp){
#   logLout <- c()
#   
#   Ut <- vp$u_mu
#   Vt <- vp$v_mu
#   
#   # plot(Ut,col = color)
#   logL <- likely(X_t,Ut,Vt)
#   print(logL)
#   
#   C <- X_t$P - X_t$N 
#   
#   Vm1 <- vp$v_mu
#   Um1 <- vp$u_mu
#   sigma_U_Pm1 <- vp$u_sigma
#   sigma_V_Pm1 <- vp$v_sigma
#   
#   for(i in 1:model$D){
#     sigma_U_Pm1[i,] <- (1/(vp$u_sigma[i,] + model$sigma_t) + 1/model$tau)
#     Um1[i,] <- sigma_U_Pm1[i,] * (1/(vp$u_sigma[i,] + model$sigma_t)) * vp$u_mu[i,]
#   }
#   
#   for(i in 1:model$D){
#     sigma_V_Pm1[i,] <- (1/(vp$v_sigma[i,] + model$sigma_t) + 1/model$tau)
#     Vm1[i,] <- sigma_V_Pm1[i,] * (1/(vp$v_sigma[i,] + model$sigma_t)) * vp$v_mu[i,]
#   }
#   
#   
#   
#   A <- vp$u_sigma * vp$u_sigma + vp$u_mu * vp$u_mu
#   for (epo in 1:model$nb_epochs) {
#     print(paste("eopch :",epo))
#     
#     #XI
#     B <- vp$v_sigma*vp$v_sigma + vp$v_mu*vp$v_mu
#     for (i in 1:model$D) {
#       
#       bob <- A[i,] * B
#       
#       vp$xi[i,] <- sqrt(rowSums(bob))
#       
#     }
#     
#     u_mu_temp <-  vp$u_mu
#     v_mu_temp <- vp$v_mu
#     u_sigma_temp <- vp$u_sigma
#     v_sigma_temp <- vp$v_sigma
#     
#     for (i in 1:model$D) {
#       # R
#       r_ui <- 0.5 * colSums(C[i,] * vp$v_mu) + (sigma_U_Pm1[i,] * Um1[i,])
#       
#       #P
#       P <- matrix(0,model$K,model$K)
#       for (j in 1:model$D) {
#         E <- diag(vp$v_sigma[j,]) + vp$v_mu[j,] %*% t(vp$v_mu[j,])
#         sig_xi <- 1/(1 + exp(-vp$xi[i,j]))
#         deux_lambda_xi <- (1/sig_xi)  * (sig_xi  - 0.5)
#         
#         P <- P + (deux_lambda_xi *  E)
#         
#       }
#       P <- P + diag(rep(1,model$K)) + sigma_U_Pm1[i,]
#       
#       Pm1 <- solve(P)
#       u_mu_temp[i,] <- Pm1 %*% r_ui
#       
#       u_sigma_temp[i,] <- diag(Pm1)
#       
#     }
#     
#     vp$u_mu <- u_mu_temp 
#     vp$u_sigma  <- u_sigma_temp
#     
#     #XI
#     
#     A <- vp$u_sigma * vp$u_sigma + vp$u_mu * vp$u_mu
#     for (i in 1:model$D) {
#       
#       bob <- A[i,] * B
#       
#       vp$xi[i,] <- sqrt(rowSums(bob))
#       # print(vp$xi[i,j])
#       
#     }
#     
#     for (i in 1:model$D) {
#       # R
#       r_vi <- 0.5 * colSums(C[,i] * vp$u_mu) + (sigma_V_Pm1[i,] * Vm1[i,])
#       
#       #P
#       P <- matrix(0,model$K,model$K)
#       for (j in 1:model$D) {
#         E <- diag(vp$u_sigma[j,]) + vp$u_mu[j,] %*% t(vp$u_mu[j,])
#         sig_xi <- 1/(1 + exp(-vp$xi[j,i]))
#         deux_lambda_xi <- (1/sig_xi)  * (sig_xi  - 0.5)
#         
#         P <- P + (deux_lambda_xi *  E)
#         
#       }
#       P <- P + diag(rep(1,model$K)) + sigma_V_Pm1[i,]
#       
#       Pm1 <- solve(P)
#       v_mu_temp[i,] <- Pm1 %*% r_vi
#       
#       v_sigma_temp[i,] <- diag(Pm1)
#       
#     }
#     
#     
#     vp$v_mu <- v_mu_temp
#     vp$v_sigma <- v_sigma_temp
#     
#     Ut <- vp$u_mu
#     Vt <- vp$v_mu
#     
#     logL <- likely(X_t,Ut,Vt)
#     print(logL)
#     logLout <- c(logLout,logL)
#   } # epo
#   
#   model$U[[model$timePosition]] <- Ut
#   model$V[[model$timePosition]] <- Vt
#   model$LL[[model$timePosition]] <- logLout
#   
#   model$timePosition <- model$timePosition + 1 
#   
#   return(list(model = model,vp = vp))
# }


draw_VP <- function(model){
  vp <- list()
  vp$u_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1/model$tau,model$K)))
  vp$v_mu <- mvrnorm(n = model$D, rep(0,model$K),diag(rep(1/model$tau,model$K)))
  vp$u_sigma <- matrix(1/model$tau,model$D,model$K)
  vp$v_sigma <- matrix(1/model$tau,model$D,model$K)
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
  
  ll <- ll +  sum(apply(U,1,function(x){log(dmvnorm(x,rep(0,model$K),diag(rep(1/model$tau,model$K))))}))
  ll <- ll +  sum(apply(V,1,function(x){log(dmvnorm(x,rep(0,model$K),diag(rep(1/model$tau,model$K))))}))
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
data <- read_delim("input/export_articles_EGC_2004_2018.csv", 
                                            "\t", escape_double = FALSE, na = "empty", 
                                            trim_ws = TRUE)
ind <- which(data$year == 2017)
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
pruned_vocabulary <- prune_vocabulary(vocabulary,  term_count_min = 3)
vectorizer <- vocab_vectorizer(pruned_vocabulary)
l <- 10
X_t$P <- as.matrix(create_tcm(iterator, vectorizer, skip_grams_window=l,skip_grams_window_context = "symmetric", weights=rep(1, l)))
X_t$P <- X_t$P + t(X_t$P)
vocab <- pruned_vocabulary$term
View(cbind(vocab,colSums(X_t$P)))
# # Custom
# Bow <- lapply(txt,Boost_tokenizer)
# word <- unlist(Bow)
# tab <- table(word)
# 
# freq <- as.vector(table(word))
# vocab <- names(tab)
# 
# naze <- which(freq<3)
# freq <- freq[-naze]
# vocab <- vocab[-naze]
# 
# sort(freq)
# dT <- matrix(0,nrow = length(Bow), ncol = length(vocab))
# 
# for (i in 1:length(Bow)) {
#   for (mot in Bow[[i]]) {
#     j <- which(vocab == tolower(mot))
#     dT[i,j] <- dT[i,j] + 1
#   }
# }
# 
# 
# X_t$P <- t(dT) %*% dT
# diag(X_t$P) <- 0
# sum(X_t$P)

save <- X_t$P

X_t$P <- save 
# Sub sampling
freq <- colSums(X_t$P)
freq <- freq /sum(freq)
# freq <-  1 - sqrt(0.00001/freq)
freq <-  1 - sqrt(0.0007/freq)

#avec fenetre
# freq <-  1 - sqrt(0.00005/freq)

X_t$P <- t(round(apply(X_t$P,1,function(x){x - freq*x})))
X_t$P
sum(X_t$P) 
X_t$N <- 5 * round((rowSums(X_t$P) %*% t(colSums(X_t$P))) / sum(X_t$P))
which(rowSums(X_t$P) == 0)
View(cbind(vocab,colSums(X_t$P)))
# X_t$P <- bob

don <- which(vocab =="données")

X_t$P <- X_t$P[-don,-don]
X_t$N <- X_t$N[-don,-don]
vocab <- vocab[-don]
# T 1

# (D,K,sigma_t,tau,nb_epochs)
model <- init_model(length(vocab),100,1,tau = 0.5,nb_epochs = 12)
vp <- draw_VP(model) 
# debug(optimZero)
out <- optimZero(X_t,model,vp)
model <- out$model
plot(model$LL[[1]])
vp <- out$vp

N <- length(vocab)
emb <- model$U[[1]] / apply(model$U[[1]],1,function(x){sqrt(sum(x * x))})
embV <- model$V[[1]] / apply(model$V[[1]],1,function(x){sqrt(sum(x * x))})
# apply(emb,1,function(x){sqrt(sum(x * x))})
cosine <- emb %*% t(embV)
res <- matrix("",N,4)
for(i in 1:N){
  cso <- cosine[i,]
  cso[i] <- -Inf
  res[i,1] <- vocab[i]
  res[i,2:4] <- vocab[sort(cso,decreasing = T,index.return=T)$ix[1:3]]
}
View(res)



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
