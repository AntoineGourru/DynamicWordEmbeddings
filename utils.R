randomSBM <- function(n,e,K,alpha,pi){
  Z = t(rmultinom(n, 1, alpha))
  cl = apply(Z, 1, which.max)
  Adj = matrix(0, n, n)
  # directed network
  cou <- 0
  while(cou < e) {
    i <- sample(1:n,1)
    j <- sample(1:n,1) 
    if(i != j) { # no self-loop
      k = cl[i] # cluster of i
      l = cl[j] # cluster of j
      add <- rbinom(1, 1, pi[k, l])
      Adj[i, j] = Adj[i, j] + add
      Adj[j, i] = Adj[j, i] + add
    }
    cou <- sum(Adj)
  }
  return(list(Adj = Adj,clusters = cl))
}