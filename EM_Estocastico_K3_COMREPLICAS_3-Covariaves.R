library('label.switching')
library("ggplot2")
library("r2excel")
library("xlsx")
library(e1071)


tempo_inicial<-proc.time()#Calcularemos quanto demoro todo o proceso
SERVER=FALSE
COMPARE_DATA=FALSE
N=300  #Tamanho da amostra Binomial
T=300 #Cumprimento da cadeia simulada
K=3    #Numero de estados ocultos
D=3    #Quantidade de Covariaveis
R<-30  #Numero de Replicas que serão executadas
tol<-0.001 #Nivel de tolerancia que estabelecemos como criterio de parada do EM Est

data_num = 1
output_num = 17
direc = "/home/gustavo/Documentos/Programas PosMestrado/NHMM PosMestrado/K3/EM Estocastico/Com_Replicas"
folder_name = paste("/",toString(output_num),"_Output_EM_Estocastico_","K",toString(K),"_","T",toString(T), sep = "")
path_completo = paste(direc, folder_name,sep = "")
ifelse(!dir.exists(path_completo), dir.create(path_completo), FALSE)
data_path = "/home/gustavo/Documentos/Programas PosMestrado/NHMM PosMestrado/TestData/K3"

arqname<-c(path_completo,"/Output_Estocastico_",toString(T),".xlsx")
nome_arquivo<-paste(arqname,collapse = "") 
nomefolha<-paste("Saida T=",toString(T))
header<-paste("Resultados para K=",toString(K), "e T=",toString(T), collapse = "")

set.seed(2011)
seeds<-sample(1:100000,R)
options(digits=8)
options(scipen=999)

#Inicializamos algumas estruturas para
#almacenar valores importantes ao longo das
#"R" replicas que serão feitas
AICc_Acum<-NULL
BIC_Acum<-NULL
Betas_Rep<-matrix(nrow = R, ncol = 18)
thetas_Rep<-matrix(nrow = R, ncol = 3)
acertos<-rep.int(0,R)

S_Table<-matrix(nrow = T, ncol = R)
S_Post_Table<-matrix(nrow = T, ncol = R)
Vies_Table<-matrix(nrow = 21, ncol = R)
CritInfo_Table<-matrix(nrow = 2, ncol = R)

#Aqui vãi iniciar o ciclo das replicas
for (p in 1:R) {
  set.seed(seeds[p]) 
  print(paste("Numero de Replica:",toString(p), collapse = ""))
  ################################
  ###        SIMULAÇÃO         ###
  ################################
  TransCount <- matrix(data = c(0,0,0,0,0,0,0,0,0), nrow = K, ncol = K)
  P0=rep(1/K,K) #Inicializamos vetor de probabilidades inciais para o HMM
  theta=c(0.40,0.5,0.60) # vetor com a probabilidade de sucesso das 2 distribuiçoes Binomiais
  Nt=rep(N,T) # número de ensaios de Bernoulli associado a cada uma das T variáveis Binomiais. Cada coloquei tudo igual mas eles podem diferentes.
  theta_hat = NULL #Variavel para estimar os thetas em cada iteração do EM Estocastico
  BetaArray = array(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), dim = c(K,D,K))
  Betas=array(dim=c(K,D,K)) # valores de beta utilizados na geração dos valores (consideranda intercepto e duas covariáveis)
  S_treino<-NULL #Inicializamos a sequência S de treinamento
  tolval=NULL #Inicializamos a variavel que almacenara a diferença nas verosim no paso i e i+1
  tolval[1]=1 #Para garantir que a primeira iteração do EMEst rode, atribuimos 1 a esse valor
  
  
  #Fazemos o valor iniciais dos Betas especificos igual a 0, para ter a função
  #de Ligação mlogit
  Betas[1,1,1]=0
  Betas[1,2,1]=0
  Betas[1,3,1]=0
  Betas[2,1,1]=-3.4
  Betas[2,2,1]=2.5
  Betas[2,3,1]=-4.2
  Betas[3,1,1]=-2.7
  Betas[3,2,1]=1.6
  Betas[3,3,1]=-3.5
  
  Betas[1,1,2]=0
  Betas[1,2,2]=0
  Betas[1,3,2]=0
  Betas[2,1,2]=-4.3
  Betas[2,2,2]=2.1
  Betas[2,3,2]=-3.4
  Betas[3,1,2]=3.2
  Betas[3,2,2]=-4.4
  Betas[3,3,2]=2.6
  
  Betas[1,1,3]=0
  Betas[1,2,3]=0
  Betas[1,3,3]=0
  Betas[2,1,3]=3.2
  Betas[2,2,3]=-2.4
  Betas[2,3,3]=4.3
  Betas[3,1,3]=-5.7
  Betas[3,2,3]=3.6
  Betas[3,3,3]=-2.8
  
  
  #####Função para gerar valores uniformes discretos
  rDiscreta<-function(v){
    u<-runif(1)
    P<-cumsum(v)
    val<-sum(P<u)+1
    return(val)}
  #####
  if (!COMPARE_DATA){
    S<-NULL #Inicializamos a sequência de estados não observaveis
    Y<-NULL #Inicializamos a sequência de valores observaveis
    X1<-NULL #Inicializamos as covariaves
    X2<-NULL #Inicializamos as covariaves
    X<-NULL #Inicializamos o vetor de covariaveis
    X1<-rnorm(T,mean = 0, sd = 1) #Geramos a covariável 1
    X2<-rnorm(T,mean = 0, sd = 1) #Geramos a covariável 2
    X<-cbind(1,X1,X2) # Creamos uma matriz de covariaveis com X0 sendo 1, e X1 e X2 .
    
    ### Finalmente, Aqui simulamos um NHMM com os 
    ### parâmetros de simulação anteriormente definidos
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
    }else{
      df1<-readRDS(paste(data_path, "/", toString(data_num),"_XTable_K",toString(K),"_T",toString(T),".rds",sep = ""))
      df2<-readRDS(paste(data_path, "/", toString(data_num),"_STable_K",toString(K),"_T",toString(T),".rds",sep = ""))
      df3<-readRDS(paste(data_path, "/", toString(data_num),"_YTable_K",toString(K),"_T",toString(T),".rds",sep = ""))
      X<-df1[,,p]
      S<-df2[,p]
      Y<-df3[,p]
  }
  
  print(S)
  
  #########################################
  #####   Procedimento de Estimação   #####
  #########################################
  
  # primeiro, escrevemos as funções que são objetivos da optimização   
  # Com o temos um array de Betas, utilizaremos tres funções 
  # para achar os valores otimos uma para a matriz Betas[,,1] 
  # uma para a matriz Betas[,,2] e uma para a matriz Betas[,,3]
  
  FSM1 <-function(params){#função a maximizar para achar os Betas_1
    resp <- sum(1 - log(1 + exp(params[1]*Xtemp11[,1] + params[2]*Xtemp11[,2] + params[3]*Xtemp11[,3])+ exp(params[4]*Xtemp11[,1]+params[5]*Xtemp11[,2] + params[6]*Xtemp11[,3]))) + sum(params[1]*Xtemp12[,1] + params[2]*Xtemp12[,2] + params[3]*Xtemp12[,3] - log(1 + exp(params[1]*Xtemp12[,1] + params[2]*Xtemp12[,2] + params[3]*Xtemp12[,3]) + exp(params[4]*Xtemp12[,1]+params[5]*Xtemp12[,2] + params[6]*Xtemp12[,3] ))) + sum(params[4]*Xtemp13[,1] + params[5]*Xtemp13[,2] + params[6]*Xtemp13[,3] - log(1 + exp(params[1]*Xtemp13[,1] + params[2]*Xtemp13[,2] + params[3]*Xtemp13[,3]) + exp(params[4]*Xtemp13[,1]+params[5]*Xtemp13[,2] + params[6]*Xtemp13[,3])))
  }
  
  FSM2 <-function(params){#função a maximizar para achar os Betas_2
    resp <- sum(1 - log(1 + exp(params[1]*Xtemp21[,1] + params[2]*Xtemp21[,2] + params[3]*Xtemp21[,3])+ exp(params[4]*Xtemp21[,1]+params[5]*Xtemp21[,2] + params[6]*Xtemp21[,3]))) + sum(params[1]*Xtemp22[,1] + params[2]*Xtemp22[,2] + params[3]*Xtemp22[,3] - log(1 + exp(params[1]*Xtemp22[,1] + params[2]*Xtemp22[,2] + params[3]*Xtemp22[,3]) + exp(params[4]*Xtemp22[,1]+params[5]*Xtemp22[,2] + params[6]*Xtemp22[,3] ))) + sum(params[4]*Xtemp23[,1] + params[5]*Xtemp23[,2] + params[6]*Xtemp23[,3] - log(1 + exp(params[1]*Xtemp23[,1] + params[2]*Xtemp23[,2] + params[3]*Xtemp23[,3]) + exp(params[4]*Xtemp23[,1]+params[5]*Xtemp23[,2] + params[6]*Xtemp23[,3])))  
  }
  
  FSM3 <-function(params){#função a maximizar para achar os Betas_3
    resp <- sum(1 - log(1 + exp(params[1]*Xtemp31[,1] + params[2]*Xtemp31[,2] + params[3]*Xtemp31[,3])+ exp(params[4]*Xtemp31[,1]+params[5]*Xtemp31[,2] + params[6]*Xtemp31[,3]))) + sum(params[1]*Xtemp32[,1] + params[2]*Xtemp32[,2] + params[3]*Xtemp32[,3] - log(1 + exp(params[1]*Xtemp32[,1] + params[2]*Xtemp32[,2] + params[3]*Xtemp32[,3]) + exp(params[4]*Xtemp32[,1]+params[5]*Xtemp32[,2] + params[6]*Xtemp32[,3] ))) + sum(params[4]*Xtemp33[,1] + params[5]*Xtemp33[,2] + params[6]*Xtemp33[,3] - log(1 + exp(params[1]*Xtemp33[,1] + params[2]*Xtemp33[,2] + params[3]*Xtemp33[,3]) + exp(params[4]*Xtemp33[,1]+params[5]*Xtemp33[,2] + params[6]*Xtemp33[,3])))
  }
  
  # geramos uma sequência não observavel de treinamento
  P_Treino=rep(1/K,K) #Vetor de probabilidade utilizadas para gerar a sequência de treino
  S_treino<-NULL # Inicializamos a sequência oculta de treinamento
    
  #Geramos uma sequência de treinamento
  for (i in 1:T) {
    S_treino[i] = rDiscreta(P_Treino)
  }
  
  # Escrevemos a função para recalcular a matriz de transição 
  # em cada iteração do algoritmo EM Estocástico.
  Mat_trans <-function(covar){#função para calcular MAtriz de Transição em cada instante de tempo
    B = matrix(nrow=K, ncol=K)
    for (j in 1:K) {
      for (i in 1:K){ 
        B[i,j] = exp(BetaArray[i,1,j] + BetaArray[i,2,j]*covar[1] + BetaArray[i,3,j]*covar[2])/(exp(BetaArray[1,1,j] + BetaArray[1,2,j]*covar[1] + BetaArray[1,3,j]*covar[2])+exp(BetaArray[2,1,j] + BetaArray[2,2,j]*covar[1] + BetaArray[2,3,j]*covar[2])+exp(BetaArray[3,1,j] + BetaArray[3,2,j]*covar[1] + BetaArray[3,3,j]*covar[2]))
      }  
    }
    return(B)
  }
  
  init1 = c(-5.0,2.85,-3.8,-1.9,-1.5,3)#Valores iniciais para os Betas_1
  init2 = c(3.2,-1.2,0.3,-2.5,-2,-2)#Valores iniciais para os Betas_2
  init3 = c(0.1,-3.1,4,-1,-3, -1)#Valores iniciais para os Betas_3
  
  val=1 # Esse é um contador para controlar o calculo do tolval em cada iteração
  
  # O EMEst sera rodado em duas etapas. Uma para estimar 
  # S e theta, e outra  para estimar os Betas
  
  #Agora executamos a primeira etapa do EM Estocástico
  while ( (tolval[val])>tol ){
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
        temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
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
    BetaArray[1,3,1]=0
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
        temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
      }
      ac2VS3 = ac2VS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
    }
    VeroSimProxima <- sum(log(choose(Nt,Y), base = exp(1))) + log(P0[S_treino[1]]) + ac2VS1 + ac2VS2 + ac2VS3 #calculo da LogVerosim
    
    val=val+1
    tolval[val]<-VeroSimProxima - VeroSimActual
    print(tolval[val])
  }
  
  # Criar algumas matrizes para fazer calculos e manipular 
  # a saida MCMC nestas matrizes, as estimativas serão 
  # reordenadas usando o metodo ECR
  
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
  
  while ( tolval[val]>tol ) {
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
    Betas_Rep[p,1] = BetaArray[2,1,1]=fit1$par[1]
    Betas_Rep[p,2] = BetaArray[2,2,1]=fit1$par[2]
    Betas_Rep[p,3] = BetaArray[2,3,1]=fit1$par[3]
    Betas_Rep[p,4] = BetaArray[3,1,1]=fit1$par[4]
    Betas_Rep[p,5] = BetaArray[3,2,1]=fit1$par[5]
    Betas_Rep[p,6] = BetaArray[3,3,1]=fit1$par[6]
    
    BetaArray[1,1,2]=0
    BetaArray[1,2,2]=0
    BetaArray[1,3,2]=0
    Betas_Rep[p,7] = BetaArray[2,1,2]=fit2$par[1]
    Betas_Rep[p,8] = BetaArray[2,2,2]=fit2$par[2]
    Betas_Rep[p,9] = BetaArray[2,3,2]=fit2$par[3]
    Betas_Rep[p,10] = BetaArray[3,1,2]=fit2$par[4]
    Betas_Rep[p,11] = BetaArray[3,2,2]=fit2$par[5]
    Betas_Rep[p,12] = BetaArray[3,3,2]=fit2$par[6]
    
    BetaArray[1,1,3]=0
    BetaArray[1,2,3]=0
    BetaArray[1,3,3]=0
    Betas_Rep[p,13] = BetaArray[2,1,3]=fit3$par[1]
    Betas_Rep[p,14] = BetaArray[2,2,3]=fit3$par[2]
    Betas_Rep[p,15] = BetaArray[2,3,3]=fit3$par[3]
    Betas_Rep[p,16] = BetaArray[3,1,3]=fit3$par[4]
    Betas_Rep[p,17] = BetaArray[3,2,3]=fit3$par[5]
    Betas_Rep[p,18] = BetaArray[3,3,3]=fit3$par[6] 
    
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
   
    tolval[val]<-VeroSimProxima-VeroSimActual
    print(tolval[val])
    val=val+1
  }###fim da segunda rodada do EM Estocastico###
  
  #############
  # Agora almacenaremos todas as quantidades de interesse 
  # em cada Replica, produtos do processo de estimação
  
  #Calculamos os Acertos em S em cada replica
  for (i in 1:T) {
    if (S[i]==S_treino[i]){
      acertos[p] = acertos[p] + 1
    }
  }
  
  # Almacenamos A sequencia S e S_treino estimada em cada
  # replica do algoritmo, para consultar depois se for necessario
  S_Post_Table[,p]<-S_treino
  S_Table[,p]<-S
  
  #Almacenamos os valores das estimativas em cada Replica
  #do algoritmo EMEst
  Betas_Rep[p,1] = BetaArray[2,1,1]
  Betas_Rep[p,2] = BetaArray[2,2,1]
  Betas_Rep[p,3] = BetaArray[2,3,1]
  Betas_Rep[p,4] = BetaArray[3,1,1]
  Betas_Rep[p,5] = BetaArray[3,2,1]
  Betas_Rep[p,6] = BetaArray[3,3,1]
  
  Betas_Rep[p,7] = BetaArray[2,1,2]
  Betas_Rep[p,8] = BetaArray[2,2,2]
  Betas_Rep[p,9] = BetaArray[2,3,2]
  Betas_Rep[p,10] = BetaArray[3,1,2]
  Betas_Rep[p,11] = BetaArray[3,2,2]
  Betas_Rep[p,12] = BetaArray[3,3,2]
  
  Betas_Rep[p,13] = BetaArray[2,1,3]
  Betas_Rep[p,14] = BetaArray[2,2,3]
  Betas_Rep[p,15] = BetaArray[2,3,3]
  Betas_Rep[p,16] = BetaArray[3,1,3]
  Betas_Rep[p,17] = BetaArray[3,2,3]
  Betas_Rep[p,18] = BetaArray[3,3,3]
  
  thetas_Rep[p,1] = theta_hat[1]
  thetas_Rep[p,2] = theta_hat[2]
  thetas_Rep[p,3] = theta_hat[3]
  
  # Calculamos os Vieses para cada um dos parametros 
  # em cada uma das replicas, e almacenamos eles
  
  Vies_Table[1,p]<-Betas[2,1,1]-BetaArray[2,1,1]
  Vies_Table[2,p]<-Betas[2,2,1]-BetaArray[2,2,1]
  Vies_Table[3,p]<-Betas[2,3,1]-BetaArray[2,3,1]
  Vies_Table[4,p]<-Betas[3,1,1]-BetaArray[3,1,1]
  Vies_Table[5,p]<-Betas[3,2,1]-BetaArray[3,2,1]
  Vies_Table[6,p]<-Betas[3,3,1]-BetaArray[3,3,1]
  
  Vies_Table[7,p]<-Betas[2,1,2]-BetaArray[2,1,2]
  Vies_Table[8,p]<-Betas[2,2,2]-BetaArray[2,2,2]
  Vies_Table[9,p]<-Betas[2,3,2]-BetaArray[2,3,2]
  Vies_Table[10,p]<-Betas[3,1,2]-BetaArray[3,1,2]
  Vies_Table[11,p]<-Betas[3,2,2]-BetaArray[3,2,2]
  Vies_Table[12,p]<-Betas[3,3,2]-BetaArray[3,3,2]
  
  Vies_Table[13,p]<-Betas[2,1,3]-BetaArray[2,1,3]
  Vies_Table[14,p]<-Betas[2,2,3]-BetaArray[2,2,3]
  Vies_Table[15,p]<-Betas[2,3,3]-BetaArray[2,3,3]
  Vies_Table[16,p]<-Betas[3,1,3]-BetaArray[3,1,3]
  Vies_Table[17,p]<-Betas[3,2,3]-BetaArray[3,2,3]
  Vies_Table[18,p]<-Betas[3,3,3]-BetaArray[3,3,3]
  
  Vies_Table[19,p]<-theta[1]-theta_hat[1]
  Vies_Table[20,p]<-theta[2]-theta_hat[2]
  Vies_Table[21,p]<-theta[3]-theta_hat[3]
  
  #Calculo do AICc e BIC
  acumLL1=0 #usaremos para calcular 1era parte da LL
  acumLL2=0 #usaremos para calcular segunda parte da LL
  acumLL3=0 #usaremos para calcular 3era parte da LL
  LL=0 #usaremos para calcular a logverosim
  
  for (i in 1:T) {#Calculo do primeiro segmento da LL
    acumLL1 = acumLL1 + Y[i]*log(theta_hat[S_treino[i]])
  }
  for (i in 1:T) {#Calculo do segundo segmento da LL
    acumLL2 = acumLL2 + (Nt[i]-Y[i])*(log(1-theta_hat[S_treino[i]], base = exp(1)))
  }
  temp=NULL
  for (i in 2:T) {#Calculo do terceiro segmento da LL
    for (g in 1:K) {
      temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[t-1]],ncol=1))
    }
    acumLL3 = acumLL3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
  }
  LL <- sum(log(choose(Nt,Y), base = exp(1))) + sum(log(P0[])) + acumLL1 + acumLL2 + acumLL3 #calculo da LogVerosim
  
  #Almacenamos os valores calculados dos criterios de info
  #Nas CritInfo_Tables, as almacenamos para poder ter acesso
  #A cada valor em cada replica, e não so à media.
  CritInfo_Table[1,p]<-BIC_Acum[p]<--2*LL + 8*log(T) #Calculo do BIC para cada Replica. O 8 aqui representa a quantidade de parametros do modelo
  CritInfo_Table[2,p]<-AICc_Acum[p]<--2*(LL) + 2*8 + (2*8^2 + 2*8)/(T - 8 - 1) #Calculo do AIC corregido
 
}

#Inicializamos algumas variaveis para
#almacenar calculos finais importantes
Porcentagem_Acertos=NULL
EQM_Thetas=NULL
EQM_Betas=NULL
Vies_Betas=NULL
Vies_Thetas=NULL
Thetas_Finais=NULL
Betas_Finais=NULL
SD_Betas=NULL
SD_Thetas=NULL
AICc_Final=NULL
BIC_Final=NULL

#Guardar as tabelas de informação de cada replica
#Para cada replcia guardamos S, S_posterior, EQM, Vies e SD 
#Para analise futuro, caso for necessario
saveRDS(S_Table,paste(path_completo, "/S_Table.rds",sep = ""))
saveRDS(S_Post_Table,paste(path_completo, "/S_Post_Table.rds",sep = ""))
saveRDS(Vies_Table,paste(path_completo, "/Vies_Table.rds",sep = ""))
saveRDS(CritInfo_Table ,paste(path_completo, "/CritInfo_Table.rds",sep = ""))

#Calculo da porcentagem de acertos na sequência S a posteriori em cada Replica
Porcentagem_Acertos<-acertos/T

#Calculo das Estimativas(medias a psoteriori) finais usando os valores das R replicas
for(i in 1:18){
  Betas_Finais[i]=mean(Betas_Rep[,i])
}
for (i in 1:K) {
  Thetas_Finais[i]=mean(thetas_Rep[,i])
}

#Calculo dos Desvios padrões usando os valores das R replicas
for (i in 1:18) {
  SD_Betas[i]=sd(Betas_Rep[,i])
}
for (i in 1:K) {
  SD_Thetas[i]=sd(thetas_Rep[,i])
}

#Calculo do Vies
Vies_Betas[1]<-Betas[2,1,1]-Betas_Finais[1]
Vies_Betas[2]<-Betas[2,2,1]-Betas_Finais[2]
Vies_Betas[3]<-Betas[2,3,1]-Betas_Finais[3]
Vies_Betas[4]<-Betas[3,1,1]-Betas_Finais[4]
Vies_Betas[5]<-Betas[3,2,1]-Betas_Finais[5]
Vies_Betas[6]<-Betas[3,3,1]-Betas_Finais[6]

Vies_Betas[7]<-Betas[2,1,2]-Betas_Finais[7]
Vies_Betas[8]<-Betas[2,2,2]-Betas_Finais[8]
Vies_Betas[9]<-Betas[2,3,2]-Betas_Finais[9]
Vies_Betas[10]<-Betas[3,1,2]-Betas_Finais[10]
Vies_Betas[11]<-Betas[3,2,2]-Betas_Finais[11]
Vies_Betas[12]<-Betas[3,3,2]-Betas_Finais[12]

Vies_Betas[13]<-Betas[2,1,3]-Betas_Finais[13]
Vies_Betas[14]<-Betas[2,2,3]-Betas_Finais[14]
Vies_Betas[15]<-Betas[2,3,3]-Betas_Finais[15]
Vies_Betas[16]<-Betas[3,1,3]-Betas_Finais[16]
Vies_Betas[17]<-Betas[3,2,3]-Betas_Finais[17]
Vies_Betas[18]<-Betas[3,3,3]-Betas_Finais[18]

for (i in 1:K) {
  Vies_Thetas[i]=theta[i]-Thetas_Finais[i]
}

#Calculo9 dos EQMs
for (i in 1:18) {
  EQM_Betas[i]=Vies_Betas[i]^2 + SD_Betas[i]^2
}
for (i in 1:K) {
  EQM_Thetas[i]=Vies_Thetas[i]^2 + SD_Thetas[i]^2
}

#Calculo da media dos AICc para o modelo
AICc_Final=mean(AICc_Acum)

#Calculo da media dos AICc para o modelo
BIC_Final=mean(BIC_Acum)

#calculo dos quanties para os parametros
quantiles<-matrix(nrow = 21, ncol = 2)
quantiles[1,]<-quantile(Betas_Rep[,1],probs = c(0.025, 0.975))
quantiles[2,]<-quantile(Betas_Rep[,2],probs = c(0.025, 0.975))
quantiles[3,]<-quantile(Betas_Rep[,3],probs = c(0.025, 0.975))
quantiles[4,]<-quantile(Betas_Rep[,4],probs = c(0.025, 0.975))
quantiles[5,]<-quantile(Betas_Rep[,5],probs = c(0.025, 0.975))
quantiles[6,]<-quantile(Betas_Rep[,6],probs = c(0.025, 0.975))
quantiles[7,]<-quantile(Betas_Rep[,7],probs = c(0.025, 0.975))
quantiles[8,]<-quantile(Betas_Rep[,8],probs = c(0.025, 0.975))
quantiles[9,]<-quantile(Betas_Rep[,9],probs = c(0.025, 0.975))
quantiles[10,]<-quantile(Betas_Rep[,10],probs = c(0.025, 0.975))
quantiles[11,]<-quantile(Betas_Rep[,11],probs = c(0.025, 0.975))
quantiles[12,]<-quantile(Betas_Rep[,12],probs = c(0.025, 0.975))
quantiles[13,]<-quantile(Betas_Rep[,13],probs = c(0.025, 0.975))
quantiles[14,]<-quantile(Betas_Rep[,14],probs = c(0.025, 0.975))
quantiles[15,]<-quantile(Betas_Rep[,15],probs = c(0.025, 0.975))
quantiles[16,]<-quantile(Betas_Rep[,16],probs = c(0.025, 0.975))
quantiles[17,]<-quantile(Betas_Rep[,17],probs = c(0.025, 0.975))
quantiles[18,]<-quantile(Betas_Rep[,18],probs = c(0.025, 0.975))
quantiles[19,]<-quantile(thetas_Rep[,1],probs = c(0.025, 0.975))
quantiles[20,]<-quantile(thetas_Rep[,2],probs = c(0.025, 0.975))
quantiles[21,]<-quantile(thetas_Rep[,3],probs = c(0.025, 0.975))

#Agora construimos os dataframes que serão exportados ao
#reporte de excel
Assimetria<-c(skewness(Betas_Rep[,1]), skewness(Betas_Rep[,2]), skewness(Betas_Rep[,3]), skewness(Betas_Rep[,4]), skewness(Betas_Rep[,5]), skewness(Betas_Rep[,6]), skewness(Betas_Rep[,7]), skewness(Betas_Rep[,8]), skewness(Betas_Rep[,9]), skewness(Betas_Rep[,10]), skewness(Betas_Rep[,11]), skewness(Betas_Rep[,12]), skewness(Betas_Rep[,13]), skewness(Betas_Rep[,14]), skewness(Betas_Rep[,15]), skewness(Betas_Rep[,16]), skewness(Betas_Rep[,17]), skewness(Betas_Rep[,18]), skewness(thetas_Rep[,1]), skewness(thetas_Rep[,2]), skewness(thetas_Rep[,3]))
Parametro<-c("Beta121", "Beta122", "Beta123", "Beta131", "Beta132", "Beta133", "Beta221", "Beta222", "Beta223", "Beta231", "Beta232", "Beta233", "Beta321", "Beta322", "Beta323", "Beta331", "Beta332", "Beta333","theta1","theta2","theta3")
Real<-c(Betas[2,1,1],Betas[2,2,1],Betas[2,3,1],Betas[3,1,1],Betas[3,2,1],Betas[3,3,1],Betas[2,1,2],Betas[2,2,2],Betas[2,3,2],Betas[3,1,2],Betas[3,2,2],Betas[3,3,2],Betas[2,1,3],Betas[2,2,3],Betas[2,3,3],Betas[3,1,3],Betas[3,2,3],Betas[3,3,3],theta[1],theta[2],theta[3])
Estimado<-c(Betas_Finais[1],Betas_Finais[2],Betas_Finais[3],Betas_Finais[4],Betas_Finais[5],Betas_Finais[6],Betas_Finais[7],Betas_Finais[8],Betas_Finais[9],Betas_Finais[10],Betas_Finais[11],Betas_Finais[12],Betas_Finais[13],Betas_Finais[14],Betas_Finais[15],Betas_Finais[16],Betas_Finais[17],Betas_Finais[18],Thetas_Finais[1],Thetas_Finais[2],Thetas_Finais[3])
SD<-c(SD_Betas[1],SD_Betas[2],SD_Betas[3],SD_Betas[4],SD_Betas[5],SD_Betas[6],SD_Betas[7],SD_Betas[8],SD_Betas[9],SD_Betas[10],SD_Betas[11],SD_Betas[12],SD_Betas[13],SD_Betas[14],SD_Betas[15],SD_Betas[16],SD_Betas[17],SD_Betas[18],SD_Thetas[1],SD_Thetas[2],SD_Thetas[3])
Vies<-c(Vies_Betas[1],Vies_Betas[2],Vies_Betas[3],Vies_Betas[4],Vies_Betas[5],Vies_Betas[6],Vies_Betas[7],Vies_Betas[8],Vies_Betas[9],Vies_Betas[10],Vies_Betas[11],Vies_Betas[12],Vies_Betas[13],Vies_Betas[14],Vies_Betas[15],Vies_Betas[16],Vies_Betas[17],Vies_Betas[18],Vies_Thetas[1],Vies_Thetas[2],Vies_Thetas[3])
EQM<-c(EQM_Betas[1],EQM_Betas[2],EQM_Betas[3],EQM_Betas[4],EQM_Betas[5],EQM_Betas[6],EQM_Betas[7],EQM_Betas[8],EQM_Betas[9],EQM_Betas[10],EQM_Betas[11],EQM_Betas[12],EQM_Betas[13],EQM_Betas[14],EQM_Betas[15],EQM_Betas[16],EQM_Betas[17],EQM_Betas[18],EQM_Thetas[1],EQM_Thetas[2],EQM_Thetas[3])
IC_95<-c(toString(quantiles[1,]),toString(quantiles[2,]),toString(quantiles[3,]), toString(quantiles[4,]), toString(quantiles[5,]), toString(quantiles[6,]),toString(quantiles[7,]),toString(quantiles[8,]),toString(quantiles[9,]),toString(quantiles[10,]),toString(quantiles[11,]), toString(quantiles[12,]), toString(quantiles[13,]), toString(quantiles[14,]),toString(quantiles[15,]),toString(quantiles[16,]),toString(quantiles[17,]),toString(quantiles[18,]),toString(quantiles[19,]),toString(quantiles[20,]),toString(quantiles[21,]))

df1<-data.frame(Parametro,Real,Estimado,SD,Vies,EQM,Assimetria,IC_95)

tempo_final<-proc.time()
tempoTOTAL<-(tempo_final - tempo_inicial)/60

Indicadores<-c("% de Acertos em S", "BIC", "AICc", "Tempo de Procesamento")
Valor<-c(mean(Porcentagem_Acertos), BIC_Final, AICc_Final, tempoTOTAL[1])
df2<-data.frame(Indicadores,Valor)

if(!SERVER){
  library("r2excel")
  wb<-createWorkbook(type = "xlsx")   #Criamos um livro de trabalho de excel
  sheet<-createSheet(wb, sheetName = nomefolha) #criamos uma folha no livro de trabalho que tinhamos criado
  xlsx.addHeader(wb, sheet, value=header,level=1, color="black", underline=1) # escrevemos um titulo na folha 
  xlsx.addLineBreak(sheet, 1) #Um linha em ranco embaixo do titulo
  xlsx.addTable(wb, sheet, df1, startCol=2)#Criamos a primeira tabela usando o primeiro dataframe
  xlsx.addLineBreak(sheet, 2)#insertamos 2 linhas em branco
  xlsx.addTable(wb, sheet, df2, startCol=2)#criamos a segunda Tabela usando o segundo dataframe
  saveWorkbook(wb, nome_arquivo)#guardamos o arquivo de excel.
}else{
  saveRDS(df1,paste(path_completo, "/df1.rds",sep = ""))
  saveRDS(df2,paste(path_completo, "/df2.rds",sep = ""))
}