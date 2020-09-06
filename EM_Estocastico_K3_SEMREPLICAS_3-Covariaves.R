library('label.switching')
library("r2excel")
library("xlsx")
library(e1071)

options(digits=8)
options(scipen=999)
num<-sample(1:10000000,1) # Seed number para conseguer os mesmos valores simulados
set.seed(num)

tempo_inicial<-proc.time()#Calcularemos quanto demoro todo o proceso

N=500 #Tamanho da amostra Binomial
T=1200 #Cumprimento da cadeia simulada
K=3   #Numero de estados ocultos
D=3   #Quantidade de Covariaveis
tol<-0.0000001 #Nivel de tolerancia que estabelecemos como criterio de parada do EM Est
tolval=NULL
tolval[1]=1

path = "/home/gustavo/Documentos/Programas PosMestrado/NHMM PosMestrado/K3/EM Estocastico/K3/Sem_Replicas/Output"
Output_Number = 4
path_completo = paste(path, toString(Output_Number),sep = "")
ifelse(!dir.exists(path_completo), dir.create(path_completo), FALSE)
VerProx<-NULL
VerAct<-NULL
arqname<-c(path_completo,"/Output_Estocastico_",toString(T),".xlsx")
nome_arquivo<-paste(arqname,collapse = "") 
nomefolha<-paste("Saida T=",toString(T))
header<-paste("Resultados para K=",toString(K), "e T=",toString(T), collapse = "")

####### Simulação #######
TransCount <- matrix(data = c(0,0,0,0,0,0,0,0,0), nrow = K, ncol = K)
P0=rep(1/K,K) #Inicializamos vetor de probabilidades inciais para o HMM
theta=c(0.20,0.5,0.80) # vetor com a probabilidade de sucesso das 2 distribuiçoes Binomiais
Nt=rep(N,T) # número de ensaios de Bernoulli associado a dada uma das T variáveis Binomiais. Cada coloquei tudo igual mas eles podem diferentes.
theta_hat = NULL #Variavel para estimar os thetas em cada iteração do EM Estocastico
BetaArray = array(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), dim = c(K,D,K))
Betas=array(dim=c(K,D,K)) # valores de beta utilizados na geração dos valores (consideranda intercepto e duas covariáveis)

#Fazemos o valor iniciais dos Betas especificos igual a 0, para ter a função
#de Ligação mlogit
Betas[1,1,1]=0
Betas[1,2,1]=0
Betas[1,3,1]=0
Betas[2,1,1]=-1.3
Betas[2,2,1]=-1.9
Betas[2,3,1]=-1.1
Betas[3,1,1]=-2.3
Betas[3,2,1]=-2.9
Betas[3,3,1]=-2.2

Betas[1,1,2]=0
Betas[1,2,2]=0
Betas[1,3,2]=0
Betas[2,1,2]=2.4
Betas[2,2,2]=2.3
Betas[2,3,2]=2.1
Betas[3,1,2]=1.1
Betas[3,2,2]=1.2
Betas[3,3,2]=2.2

Betas[1,1,3]=0
Betas[1,2,3]=0
Betas[1,3,3]=0
Betas[2,1,3]=-1.4
Betas[2,2,3]=-1.2
Betas[2,3,3]=-1.3
Betas[3,1,3]=1.8
Betas[3,2,3]=1.9
Betas[3,3,3]=1.7

#####Função para gerar valores uniformes discretos
rDiscreta<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  return(val)}
#####

#######   Escrevemos as funções que serão o objetivo da optimização   ######
# Com o temos um array de Betas, utilizaremos tres funções para achar os valores otimos
# Uma para a matriz Betas[,,1] uma para a matriz Betas[,,2] e uma para 
# a matriz Betas[,,3]

FSM1 <-function(params){#função a maximizar para achar os Betas_1
  resp <- sum(1 - log(1 + exp(params[1]*Xtemp11[,1] + params[2]*Xtemp11[,2] + params[3]*Xtemp11[,3])+ exp(params[4]*Xtemp11[,1]+params[5]*Xtemp11[,2] + params[6]*Xtemp11[,3]))) + sum(params[1]*Xtemp12[,1] + params[2]*Xtemp12[,2] + params[3]*Xtemp12[,3] - log(1 + exp(params[1]*Xtemp12[,1] + params[2]*Xtemp12[,2] + params[3]*Xtemp12[,3]) + exp(params[4]*Xtemp12[,1]+params[5]*Xtemp12[,2] + params[6]*Xtemp12[,3] ))) + sum(params[4]*Xtemp13[,1] + params[5]*Xtemp13[,2] + params[6]*Xtemp13[,3] - log(1 + exp(params[1]*Xtemp13[,1] + params[2]*Xtemp13[,2] + params[3]*Xtemp13[,3]) + exp(params[4]*Xtemp13[,1]+params[5]*Xtemp13[,2] + params[6]*Xtemp13[,3])))
}

FSM2 <-function(params){#função a maximizar para achar os Betas_2
  resp <- sum(1 - log(1 + exp(params[1]*Xtemp21[,1] + params[2]*Xtemp21[,2] + params[3]*Xtemp21[,3])+ exp(params[4]*Xtemp21[,1]+params[5]*Xtemp21[,2] + params[6]*Xtemp21[,3]))) + sum(params[1]*Xtemp22[,1] + params[2]*Xtemp22[,2] + params[3]*Xtemp22[,3] - log(1 + exp(params[1]*Xtemp22[,1] + params[2]*Xtemp22[,2] + params[3]*Xtemp22[,3]) + exp(params[4]*Xtemp22[,1]+params[5]*Xtemp22[,2] + params[6]*Xtemp22[,3] ))) + sum(params[4]*Xtemp23[,1] + params[5]*Xtemp23[,2] + params[6]*Xtemp23[,3] - log(1 + exp(params[1]*Xtemp23[,1] + params[2]*Xtemp23[,2] + params[3]*Xtemp23[,3]) + exp(params[4]*Xtemp23[,1]+params[5]*Xtemp23[,2] + params[6]*Xtemp23[,3])))  
}

FSM3 <-function(params){#função a maximizar para achar os Betas_3
  resp <- sum(1 - log(1 + exp(params[1]*Xtemp31[,1] + params[2]*Xtemp31[,2] + params[3]*Xtemp31[,3])+ exp(params[4]*Xtemp31[,1]+params[5]*Xtemp31[,2] + params[6]*Xtemp31[,3]))) + sum(params[1]*Xtemp32[,1] + params[2]*Xtemp32[,2] + params[3]*Xtemp32[,3] - log(1 + exp(params[1]*Xtemp32[,1] + params[2]*Xtemp32[,2] + params[3]*Xtemp32[,3]) + exp(params[4]*Xtemp32[,1]+params[5]*Xtemp32[,2] + params[6]*Xtemp32[,3] ))) + sum(params[4]*Xtemp33[,1] + params[5]*Xtemp33[,2] + params[6]*Xtemp33[,3] - log(1 + exp(params[1]*Xtemp33[,1] + params[2]*Xtemp33[,2] + params[3]*Xtemp33[,3]) + exp(params[4]*Xtemp33[,1]+params[5]*Xtemp33[,2] + params[6]*Xtemp33[,3])))
}

################################
###        SIMULAÇÃO         ###
################################

S<-NULL #Inicializamos a sequência de estados não observaveis
Y<-NULL #Inicializamos a sequência de valores observaveis
X1<-NULL #Inicializamos as covariaves
X2<-NULL #Inicializamos as covariaves
X<-NULL #Inicializamos o vetor de covariaveis
S_treino<-NULL #Inicializamos a sequência S de treinamento
X1<-rnorm(T,mean = 0, sd = 1) # covariável 1
X2<-rnorm(T,mean = 0, sd = 1) # covariável 1
X<-cbind(1,X1,X2) # Creamos uma matriz de covariaveis com X0 sendo 1 e X1 sendo 


S[1]<-rDiscreta(P0) #O valor para o primeiro estado oculto
Y[1]<-rbinom(1,Nt[1],theta[S[1]])# O valor para o primeiro valor observavel
#Geramos a Sequencia S de cumprimento T e um valor observavel Y para cada S.
for (t in 2:T){
  prob<-NULL
  for (i in 1:K) prob[i]<-exp(X[t,]%*%matrix(Betas[i,,S[t-1]],ncol=1))
  prob<-prob/sum(prob)
  S[t]<-rDiscreta(prob)
  Y[t]<-rbinom(1,Nt[t],theta[S[t]])
}

#########################################
#####   Procedimento de Estimação   #####
#########################################

#######Primeiro geramos uma sequência não observavel de treinamento######
P_Treino=rep(1/K,K) #Vetor de probabilidade utilizadas para gerar a sequência de treino
S_treino<-NULL # Inicializamos a sequência oculta de treinamento
init1 = c(-5.0,2.85,-3.8,-1.9,-1.5,3)#Valores iniciais para os Betas_1
init2 = c(3.2,-1.2,0.3,-2.5,-2,-2)#Valores iniciais para os Betas_2
init3 = c(0.1,-3.1,4,-1,-3, -1)#Valores iniciais para os Betas_3

#Geramos uma sequência de treinamento
for (i in 1:T) {
  S_treino[i] = rDiscreta(P_Treino)
}

## Escrevemos a função para recalcular a matriz de transição em cada iteração
## do algoritmo EM Estocástico.
Mat_trans <-function(covar){#função para calcular MAtriz de Transição em cada instante de tempo
  B = matrix(nrow=K, ncol=K)
  for (j in 1:K) {
    for (i in 1:K){ 
      B[i,j] = exp(BetaArray[i,1,j] + BetaArray[i,2,j]*covar[1] + BetaArray[i,3,j]*covar[2])/(exp(BetaArray[1,1,j] + BetaArray[1,2,j]*covar[1] + BetaArray[1,3,j]*covar[2])+exp(BetaArray[2,1,j] + BetaArray[2,2,j]*covar[1] + BetaArray[2,3,j]*covar[2])+exp(BetaArray[3,1,j] + BetaArray[3,2,j]*covar[1] + BetaArray[3,3,j]*covar[2]))
    }  
  }
  return(B)
}

val=1
tolval[1]=1
#Agora executamos o Algoritmo EM Estocástico
while ( (abs(tolval[val]))>tol ){
  #VeroSimActual=VeroSimProxima
  #Aqui devemos calcular a diferença entre a L.V. em na iteração atual e na anterior  
  acVS1=0
  acVS2=0
  acVS3=0
  VeroSimActual=0
  
  for (k in 1:K) {id = S_treino == k
  theta_hat[k] = sum(id*Y)/sum(id*Nt)
  }
  
  for (i in 1:T) {#Calculo do primeiro segmento da LL
     acVS1 = acVS1 + Y[i]*log(theta_hat[S_treino[i]])
  }
  for (i in 1:T) {#Calculo do segundo segmento da LL
    acVS2 = acVS2 + (Nt[i]-Y[i])*(log(1-theta_hat[S_treino[i]]))
  }
  temp=NULL
  for (i in 2:T) {#Calculo do terceiro segmento da LL
    for (g in 1:K) {
      temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[t-1]],ncol=1))
  }
    acVS3 = acVS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
  }
  VeroSimActual <- sum(log(choose(Nt,Y), base = exp(1))) + log(P0[S_treino[1]]) + acVS1 + acVS2 + acVS3 #calculo da LogVerosim
  
  #Calculamos a sequência S_treino utilizando os Betas
  #Atualizados na iteração passada e os valores observados Y
  S_treino[1]=which.max(dbinom(Y[1], Nt[1], theta_hat))
  for (i in 2:T) {
    A_hat_t = Mat_trans(X[i,])
    prob<-(A_hat_t[S_treino[i], ]*dbinom(Y[i], Nt[i], theta_hat))/sum(A_hat_t[S_treino[i], ]*dbinom(Y[i], Nt[i], theta_hat))
    S_treino[i]=rDiscreta(prob)
    #S_treino[i]=which.max(A_hat_t[S_treino[i-1], ]*dbinom(Y[i], Nt[i], theta_hat))
  }
  
  #####################################
  #Este segmento de codigo testa se aconteceram todas as transições possiveis
  #No caso que elas não tinham acontecido, as que
  #não aconteceram são forçadas a acontecer
  TransCount <- matrix(data = c(0,0,0,0,0,0,0,0,0), nrow = K, ncol = K)
  for (i in 2:T) {
    for (j in 1:K) {
      for (k in 1:K) {
        if (S_treino[i]==j && S_treino[i-1]==k)
          TransCount[k,j]=TransCount[k,j]+1
      }
    }
  }
  
  for (j in 1:K) {
    for (k in 1:K) {
      if (TransCount[k,j]==0){
        positions = sample(2:T, 4)
        for (d in 1:4) {
          S_treino[positions[d]]=j
          S_treino[positions[d]-1]=k
        }
      }
    }
  }
  
  #### Aqui inicia a filtragem dos dados para cada iteração
  Xtemp11<-NULL
  Xtemp12<-NULL
  Xtemp13<-NULL
  Xtemp21<-NULL
  Xtemp22<-NULL
  Xtemp23<-NULL
  Xtemp31<-NULL
  Xtemp32<-NULL
  Xtemp33<-NULL
  
  for (t in 2:T) {
    #filtros indo para o Estado # 1
    if(S_treino[t]%in%1 && S_treino[t-1]%in%1)
      Xtemp11<-rbind(Xtemp11, X[t,])
    
    if(S_treino[t]%in%1 && S_treino[t-1]%in%2)
      Xtemp21<-rbind(Xtemp21, X[t,])
    
    if(S_treino[t]%in%1 && S_treino[t-1]%in%3)
      Xtemp31<-rbind(Xtemp31, X[t,])
    
    #Filtros indo para o Estado # 2
    if(S_treino[t]%in%2 && S_treino[t-1]%in%1)
      Xtemp12<-rbind(Xtemp12, X[t,])
    
    if(S_treino[t]%in%2 && S_treino[t-1]%in%2)
      Xtemp22<-rbind(Xtemp22, X[t,])
    
    if(S_treino[t]%in%2 && S_treino[t-1]%in%3)
      Xtemp32<-rbind(Xtemp32, X[t,])
    
    #Filtros indo para o Estado # 3
    if(S_treino[t]%in%3 && S_treino[t-1]%in%1)
      Xtemp13<-rbind(Xtemp13, X[t,])
    
    if(S_treino[t]%in%3 && S_treino[t-1]%in%2)
      Xtemp23<-rbind(Xtemp23, X[t,])
    
    if(S_treino[t]%in%3 && S_treino[t-1]%in%3)
      Xtemp33<-rbind(Xtemp33, X[t,])
  }
  Xtemp11
  Xtemp12
  Xtemp13
  Xtemp21
  Xtemp22
  Xtemp23
  Xtemp31
  Xtemp32
  Xtemp33
  
  ##O ajuste para estimar os parâmetros de transição é
  ##feito aqui usando a função optim e os valores das
  #covariaveis filtradas
  fit1 <- optim(par = init1, fn = FSM1, control = list(fnscale=-1), method = "Nelder-Mead", hessian = FALSE)
  fit2 <- optim(par = init2, fn = FSM2, control = list(fnscale=-1), method = "Nelder-Mead", hessian = FALSE)
  fit3 <- optim(par = init3, fn = FSM3, control = list(fnscale=-1), method = "Nelder-Mead", hessian = FALSE)
  
  # Aqui atribuimos os valores estimados dos parâmetros de 
  # transição a um array que sera utilizado para recalcular 
  # a sequência S_treino na seguinte iteração do EM Est. 
  # Em outras palavras, aqui acontece a ATUALIZAÇÃO dos parâmetros de transição.
  
  BetaArray[1,1,1]=0
  BetaArray[1,2,1]=0
  BetaArray[1,2,1]=0
  BetaArray[2,1,1]=fit1$par[1]
  BetaArray[2,2,1]=fit1$par[2]
  BetaArray[2,3,1]=fit1$par[3]
  BetaArray[3,1,1]=fit1$par[4]
  BetaArray[3,2,1]=fit1$par[5]
  BetaArray[3,3,1]=fit1$par[6]
  
  BetaArray[1,1,2]=0
  BetaArray[1,2,2]=0
  BetaArray[1,3,2]=0
  BetaArray[2,1,2]=fit2$par[1]
  BetaArray[2,2,2]=fit2$par[2]
  BetaArray[2,3,2]=fit2$par[3]
  BetaArray[3,1,2]=fit2$par[4]
  BetaArray[3,2,2]=fit2$par[5]
  BetaArray[3,3,2]=fit2$par[6]
  
  BetaArray[1,1,3]=0
  BetaArray[1,2,3]=0
  BetaArray[1,3,3]=0
  BetaArray[2,1,3]=fit3$par[1]
  BetaArray[2,2,3]=fit3$par[2]
  BetaArray[2,3,3]=fit3$par[3]
  BetaArray[3,1,3]=fit3$par[4]
  BetaArray[3,2,3]=fit3$par[5]
  BetaArray[3,3,3]=fit3$par[6]
  
  ac2VS1=0
  ac2VS2=0
  ac2VS3=0
  VeroSimProxima=0
  
  for (i in 1:T) {#Calculo do primeiro segmento da LL
    ac2VS1 = ac2VS1 + Y[i]*log(theta_hat[S_treino[i]])
  }
  for (i in 1:T) {#Calculo do segundo segmento da LL
    ac2VS2 = ac2VS2 + (Nt[i]-Y[i])*(log(1-theta_hat[S_treino[i]]))
  }
  temp=NULL
  for (i in 2:T) {#Calculo do terceiro segmento da LL
    for (g in 1:K) {
      temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[t-1]],ncol=1))
    }
    ac2VS3 = ac2VS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
  }
  VeroSimProxima <- sum(log(choose(Nt,Y), base = exp(1))) + log(P0[S_treino[1]]) + ac2VS1 + ac2VS2 + ac2VS3 #calculo da LogVerosim
  
  val=val+1
  VerAct[val]<-VeroSimActual
  VerProx[val]<-VeroSimProxima
  tolval[val]<-VeroSimProxima - VeroSimActual
  print(tolval[val])
  
}#######Fim da primeira rodada do EM Estocastico#######

#Criar algumas matrizes para fazer calculos e manipular a saida MCMC
#nestas matrizes, as estimativas serão reordenadas usando o metodo ECR
mat_thetar<-matrix(nrow = 1, ncol = K)
reorder_S<-matrix(nrow = 1, ncol = T)
mat_S<-matrix(nrow = 1, ncol = T)
mat_S[1,]<-S_treino
zpvt_S = S #Como pivot para o metodo ECR usamos o S original
perms_S = ecr(zpivot = zpvt_S, z = mat_S, K = 3)# aplicamos o metodo ECR que retornara as permutações das dos estados ocultos que devem ser utilizadas para reordenar a saida do algoritmo bayesiano

# Reordenamos a saido do algoritmo EMEst usando as 
# permutações fornecidas pelo ECR para K=3 
# só rerotulamos a Sequência S_treino, e reordenamos os Thetas
# Os Betas serão estimados usando a sequência S_Treino rerotulada
# e os Thetas, na segunda etapa do EMEst

for (i in 1:1) {
  for (j in 1:T) {
    if(S_treino[j]!=S[j] && ((perms_S$permutations[i,1]==2 && perms_S$permutations[i,2]==3 && perms_S$permutations[i,3]==1) | (perms_S$permutations[i,1]==3 && perms_S$permutations[i,2]==1 && perms_S$permutations[i,3]==2))){
      S_treino[j]=perms_S$permutations[i,perms_S$permutations[i,S_treino[j]]]
    }
    
    else {
      S_treino[j]=perms_S$permutations[i,perms_S$permutations[i,S[j]]]
    }
  }
  theta_hat<-theta_hat[perms_S$permutations[i,]]
}

# repetimos o EM Estocastico, porque para K=3
# Entraremos com a sequência S_treino estimada na rodada 
# anterior E ja rotulada corretamente usando o ECR Para resolver
# o problema dos Parametros Fantasmas e a troca de rotulos

VeroSimProxima=1
VeroSimActual=0
val=1
tolval[1]=1
while ( abs(tolval[val])>tol ) {
  #Aqui devemos calcular a diferença entre a L.V. em na iteração atual e na anterior  
  acVS1=0
  acVS2=0
  acVS3=0
  VeroSimActual=0
  
  for (i in 1:T) {#Calculo do primeiro segmento da LL
    acVS1 = acVS1 + Y[i]*log(theta_hat[S_treino[i]])
  }
  for (i in 1:T) {#Calculo do segundo segmento da LL
    acVS2 = acVS2 + (Nt[i]-Y[i])*(log(1-theta_hat[S_treino[i]]))
  }
  temp=NULL
  for (i in 2:T) {#Calculo do terceiro segmento da LL
    for (g in 1:K) {
      temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[t-1]],ncol=1))
    }
    acVS3 = acVS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
  }
  VeroSimActual <- sum(log(choose(Nt,Y), base = exp(1))) + log(P0[S_treino[1]]) + acVS1 + acVS2 + acVS3 #calculo da LogVerosim
  
  #filtragem dos dados
  Xtemp11<-NULL
  Xtemp12<-NULL
  Xtemp13<-NULL
  Xtemp21<-NULL
  Xtemp22<-NULL
  Xtemp23<-NULL
  Xtemp31<-NULL
  Xtemp32<-NULL
  Xtemp33<-NULL
  
  for (t in 2:T) {
    #filtros indo para o Estado # 1
    if(S_treino[t]%in%1 && S_treino[t-1]%in%1)
      Xtemp11<-rbind(Xtemp11, X[t,])
    
    if(S_treino[t]%in%1 && S_treino[t-1]%in%2)
      Xtemp21<-rbind(Xtemp21, X[t,])
    
    if(S_treino[t]%in%1 && S_treino[t-1]%in%3)
      Xtemp31<-rbind(Xtemp31, X[t,])
    
    #Filtros indo para o Estado # 2
    if(S_treino[t]%in%2 && S_treino[t-1]%in%1)
      Xtemp12<-rbind(Xtemp12, X[t,])
    
    if(S_treino[t]%in%2 && S_treino[t-1]%in%2)
      Xtemp22<-rbind(Xtemp22, X[t,])
    
    if(S_treino[t]%in%2 && S_treino[t-1]%in%3)
      Xtemp32<-rbind(Xtemp32, X[t,])
    
    #Filtros indo para o Estado # 3
    if(S_treino[t]%in%3 && S_treino[t-1]%in%1)
      Xtemp13<-rbind(Xtemp13, X[t,])
    
    if(S_treino[t]%in%3 && S_treino[t-1]%in%2)
      Xtemp23<-rbind(Xtemp23, X[t,])
    
    if(S_treino[t]%in%3 && S_treino[t-1]%in%3)
      Xtemp33<-rbind(Xtemp33, X[t,])
  }
  
  fit1 <- optim(par = init1, fn = FSM1, control = list(fnscale=-1), method = "Nelder-Mead", hessian = FALSE)
  fit2 <- optim(par = init2, fn = FSM2, control = list(fnscale=-1), method = "Nelder-Mead", hessian = FALSE)
  fit3 <- optim(par = init3, fn = FSM3, control = list(fnscale=-1), method = "Nelder-Mead", hessian = FALSE)
 
  BetaArray[1,1,1]=0
  BetaArray[1,2,1]=0
  BetaArray[1,2,1]=0
  BetaArray[2,1,1]=fit1$par[1]
  BetaArray[2,2,1]=fit1$par[2]
  BetaArray[2,3,1]=fit1$par[3]
  BetaArray[3,1,1]=fit1$par[4]
  BetaArray[3,2,1]=fit1$par[5]
  BetaArray[3,3,1]=fit1$par[6]
  
  BetaArray[1,1,2]=0
  BetaArray[1,2,2]=0
  BetaArray[1,3,2]=0
  BetaArray[2,1,2]=fit2$par[1]
  BetaArray[2,2,2]=fit2$par[2]
  BetaArray[2,3,2]=fit2$par[3]
  BetaArray[3,1,2]=fit2$par[4]
  BetaArray[3,2,2]=fit2$par[5]
  BetaArray[3,3,2]=fit2$par[6]
  
  BetaArray[1,1,3]=0
  BetaArray[1,2,3]=0
  BetaArray[1,3,3]=0
  BetaArray[2,1,3]=fit3$par[1]
  BetaArray[2,2,3]=fit3$par[2]
  BetaArray[2,3,3]=fit3$par[3]
  BetaArray[3,1,3]=fit3$par[4]
  BetaArray[3,2,3]=fit3$par[5]
  BetaArray[3,3,3]=fit3$par[6] 
  
  ac2VS1=0
  ac2VS2=0
  ac2VS3=0
  VeroSimProxima=0
  
  for (i in 1:T) {#Calculo do primeiro segmento da LL
    ac2VS1 = ac2VS1 + Y[i]*log(theta_hat[S_treino[i]])
  }
  for (i in 1:T) {#Calculo do segundo segmento da LL
    ac2VS2 = ac2VS2 + (Nt[i]-Y[i])*(log(1-theta_hat[S_treino[i]]))
  }
  temp=NULL
  for (i in 2:T) {#Calculo do terceiro segmento da LL
    for (g in 1:K) {
      temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[t-1]],ncol=1))
    }
    ac2VS3 = ac2VS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
  }
  VeroSimProxima <- sum(log(choose(Nt,Y), base = exp(1))) + P0[S_treino[1]] + ac2VS1 + ac2VS2 + ac2VS3 #calculo da LogVerosim
  VerAct[val]<-VeroSimActual
  VerProx[val]<-VeroSimProxima
  tolval[val]<-VeroSimProxima-VeroSimActual
  print(tolval[val])
  val=val+1
}###fim da segunda rodada do EM Estocastico###

Beta_Post_Array=array(round(runif((K*D*K),0,1),1),c(K,D,K))#Utilizamos este arreglo para reorganizar os Betas numa estrutura conveniente na hora de calcular a LL e o BIC
#Inicializamos algumas variaveis para
#almacenar calculos finais importantes
Porcentagem_Acertos=NULL
Vies_Betas=NULL
Vies_Thetas=NULL
Thetas_Finais=NULL
Betas_Finais=NULL

#Comparamos S com S_treino para saber a porcentagem de acertos
acertos=0
for (i in 1:T) {
  if (S[i]==S_treino[i]){
    acertos = acertos + 1
  }
}

Porcentagem_Acertos<-acertos/T

#Almaceno os valores finais das estimativas
Beta_Post_Array[2,1,1]= fit1$par[1]
Beta_Post_Array[2,2,1]= fit1$par[2]
Beta_Post_Array[2,3,1]= fit1$par[3]
Beta_Post_Array[3,1,1]= fit1$par[4]
Beta_Post_Array[3,2,1]= fit1$par[5]
Beta_Post_Array[3,3,1]= fit1$par[6]

Beta_Post_Array[2,1,2]= fit2$par[1]
Beta_Post_Array[2,2,2]= fit2$par[2]
Beta_Post_Array[2,3,2]= fit2$par[3]
Beta_Post_Array[3,1,2]= fit2$par[4]
Beta_Post_Array[3,2,2]= fit2$par[5]
Beta_Post_Array[3,3,2]= fit2$par[6]

Beta_Post_Array[2,1,3]= fit3$par[1]
Beta_Post_Array[2,2,3]= fit3$par[2]
Beta_Post_Array[2,3,3]= fit3$par[3]
Beta_Post_Array[3,1,3]= fit3$par[4]
Beta_Post_Array[3,2,3]= fit3$par[5]
Beta_Post_Array[3,3,3]= fit3$par[6]

Thetas_Finais<-theta_hat



#Calculo do AICc e BIC

acumLL1=0 #usaremos para calcular 1era parte da LL
acumLL2=0 #usaremos para calcular segunda parte da LL
acumLL3=0 #usaremos para calcular 3era parte da LL
LL=0 #usaremos para calcular a logverosim
AICc=0
BIC=0

for (i in 1:T) {#Calculo do primeiro segmento da LL
  acumLL1 = acumLL1 + Y[i]*log(Thetas_Finais[S_treino[i]])
}
for (i in 1:T) {#Calculo do segundo segmento da LL
  acumLL2 = acumLL2 + (Nt[i]-Y[i])*(log(1-Thetas_Finais[S_treino[i]], base = exp(1)))
}
temp=NULL
for (i in 2:T) {#Calculo do terceiro segmento da LL
  for (g in 1:K) {
    temp[g]<-exp(X[i,]%*%matrix(Beta_Post_Array[g,,S_treino[t-1]],ncol=1))
  }
  acumLL3 = acumLL3 + (X[i,]%*%matrix(Beta_Post_Array[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
}
LL <- sum(log(choose(Nt,Y), base = exp(1))) + sum(log(P0[])) + acumLL1 + acumLL2 + acumLL3 #calculo da LogVerosim

BIC<--2*LL + 8*log(T) #Calculo do BIC para cada Replica. O 8 aqui representa a quantidade de parametros do modelo
AICc<--2*(LL) + 2*8 + (2*8^2 + 2*8)/(T - 8 - 1) #Calculo do AIC corregido
#####################################################################
#####################################################################


#Calculo do Vies
Vies_Betas[1]<-Betas[2,1,1]-Beta_Post_Array[2,1,1]
Vies_Betas[2]<-Betas[2,2,1]-Beta_Post_Array[2,2,1]
Vies_Betas[3]<-Betas[2,3,1]-Beta_Post_Array[2,3,1]
Vies_Betas[4]<-Betas[3,1,1]-Beta_Post_Array[3,1,1]
Vies_Betas[5]<-Betas[3,2,1]-Beta_Post_Array[3,2,1]
Vies_Betas[6]<-Betas[3,3,1]-Beta_Post_Array[3,3,1]

Vies_Betas[7]<-Betas[2,1,2]-Beta_Post_Array[2,1,2]
Vies_Betas[8]<-Betas[2,2,2]-Beta_Post_Array[2,2,2]
Vies_Betas[9]<-Betas[2,3,2]-Beta_Post_Array[2,3,2]
Vies_Betas[10]<-Betas[3,1,2]-Beta_Post_Array[3,1,2]
Vies_Betas[11]<-Betas[3,2,2]-Beta_Post_Array[3,2,2]
Vies_Betas[12]<-Betas[3,3,2]-Beta_Post_Array[3,3,2]

Vies_Betas[13]<-Betas[2,1,3]-Beta_Post_Array[2,1,3]
Vies_Betas[14]<-Betas[2,2,3]-Beta_Post_Array[2,2,3]
Vies_Betas[15]<-Betas[2,3,3]-Beta_Post_Array[2,3,3]
Vies_Betas[16]<-Betas[3,1,3]-Beta_Post_Array[3,1,3]
Vies_Betas[17]<-Betas[3,2,3]-Beta_Post_Array[3,2,3]
Vies_Betas[18]<-Betas[3,3,3]-Beta_Post_Array[3,3,3]

for (i in 1:K) {
  Vies_Thetas[i]=theta[i]-theta_hat[i]
}

Parametro<-c("Beta121", "Beta122", "Beta123", "Beta131", "Beta132", "Beta133", "Beta221", "Beta222", "Beta223", "Beta231", "Beta232", "Beta233", "Beta321", "Beta322", "Beta323", "Beta331", "Beta332", "Beta333","theta1","theta2","theta3")
Real<-c(Betas[2,1,1],Betas[2,2,1],Betas[2,3,1],Betas[3,1,1],Betas[3,2,1],Betas[3,3,1],Betas[2,1,2],Betas[2,2,2],Betas[2,3,2],Betas[3,1,2],Betas[3,2,2],Betas[3,3,2],Betas[2,1,3],Betas[2,2,3],Betas[2,3,3],Betas[3,1,3],Betas[3,2,3],Betas[3,3,3],theta[1],theta[2],theta[3])
Estimado<-c(Beta_Post_Array[2,1,1],Beta_Post_Array[2,2,1],Beta_Post_Array[2,3,1],Beta_Post_Array[3,1,1],Beta_Post_Array[3,2,1],Beta_Post_Array[3,3,1],Beta_Post_Array[2,1,2],Beta_Post_Array[2,2,2],Beta_Post_Array[2,3,2],Beta_Post_Array[3,1,2],Beta_Post_Array[3,2,2],Beta_Post_Array[3,3,2],Beta_Post_Array[2,1,3],Beta_Post_Array[2,2,3],Beta_Post_Array[2,3,3],Beta_Post_Array[3,1,3],Beta_Post_Array[3,2,3],Beta_Post_Array[3,3,3],Thetas_Finais[1],Thetas_Finais[2],Thetas_Finais[3])
Vies<-c(Vies_Betas[1],Vies_Betas[2],Vies_Betas[3],Vies_Betas[4],Vies_Betas[5],Vies_Betas[6],Vies_Betas[7],Vies_Betas[8],Vies_Betas[9],Vies_Betas[10],Vies_Betas[11],Vies_Betas[12],Vies_Betas[13],Vies_Betas[14],Vies_Betas[15],Vies_Betas[16],Vies_Betas[17],Vies_Betas[18],Vies_Thetas[1],Vies_Thetas[2],Vies_Thetas[3])


df1<-data.frame(Parametro,Real,Estimado,Vies)

tempo_final<-proc.time()
tempoTOTAL<-(tempo_final - tempo_inicial)/60

Indicadores<-c("% de Acertos em S", "BIC", "AICc", "Tempo de Procesamento", "Seed Value")
Valor<-c(Porcentagem_Acertos, BIC, AICc, tempoTOTAL[1], num)
df2<-data.frame(Indicadores,Valor)

wb<-createWorkbook(type = "xlsx")#Criamos um livro de trabalho de excel
sheet<-createSheet(wb, sheetName = nomefolha) #criamos uma folha no livro de trabalh

xlsx.addHeader(wb, sheet, value=header,level=1, color="black", underline=1) # escrevemos um titulo na folha 
xlsx.addLineBreak(sheet, 1) #Um linha em branco embaixo do titulo
xlsx.addTable(wb, sheet, df1, startCol=2)#Criamos a primeira tabela usando o primeiro dataframe
xlsx.addLineBreak(sheet, 2)#insertamos 2 linhas em branco
xlsx.addTable(wb, sheet, df2, startCol=2)#criamos a segunda Tabela usando o segundo dataframe
saveWorkbook(wb, nome_arquivo)#guardamos o arquivo de excel.

