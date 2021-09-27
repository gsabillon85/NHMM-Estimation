

rDiscrete<-function(prob){
  u<-runif(1)
  P<-cumsum(prob)
  val<-sum(P<u)+1
  return(val)
}

DataGeneratorNHMM <- function(X, Beta, K=2, T=100, theta, Nt){
  S<-NULL
  Y<-NULL
  S[1]<-rDiscrete(rep(1/K,K)) #O valor para o primeiro estado oculto
  Y[1]<-rbinom(1,Nt[1],theta[S[1]])# O valor para o primeiro valor observavel
  #Geramos a Sequencia S de cumprimento T e um valor observavel Y para cada S.
  for (t in 2:T){
    prob<-NULL
    for (i in 1:K) prob[i]<-exp(X[t,]%*%matrix(Beta[i,,S[t-1]],ncol=1))
    prob<-prob/sum(prob)
    S[t]<-rDiscrete(prob)
    Y[t]<-rbinom(1,Nt[t],theta[S[t]])
  }
  data<-cbind(Y,S)
  return(data)
}

Beta.Names <- function(K=2, D){
  BetasNames <- NULL
  for (i in 1:K){
    for (j in 2:K){
      for (d in 1:D){
        BetasNames <- rbind(BetasNames, paste("Beta[",toString(j),",", toString(d),",",toString(i),"]", sep = "") )
      }
    }
  }
  BetasNames<-t(BetasNames)[1,]
  return(BetasNames)
}

PresetThetas <- function(K=2){
  offset <- 0.90/(K-1)
  vec<-NULL
  vec[1]<-0.05 
  for (i in 2:K) {
    vec[i] <- vec[1] + offset*(i-1)
  }
  return(vec)
}
