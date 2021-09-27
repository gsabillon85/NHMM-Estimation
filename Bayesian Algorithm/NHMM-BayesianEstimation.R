require(xlsx)
library('rjags')
library('coda')
library('label.switching')
library(e1071)
library(gridExtra)
library(gridBase)
library(grid)
library("ggmcmc")
library("bayesplot")
library("ggplot2")
library(ggplotify)

source("Auxiliar/AuxiliarFunctions.R")
InitialTime <- proc.time()

# Reading the data
#-----
DataPath = "Bayesian Estimation/"
Name = "RainData.xlsx"
FileName = paste(DataPath, Name, sep = "")
df<-read.xlsx(file = FileName, sheetIndex = 1)
#----

# Some preleminary treatment of the data
#------
df$Month<-as.numeric(as.character(df$Month))
df$Year<-as.numeric(as.character(df$Year))
df$Temperature<-as.numeric(as.character(df$Temperature))
df$RelativeHumidity<-as.numeric(as.character(df$RelativeHumidity))
df$RainyDays<-as.numeric(as.character(df$RainyDays))
#----

# Set the number of hidden states, covariates and Binomial sample size
# NOTE: For this application with rainfall data, we tested multpĺie situations
# The situation that produced the best fitted model was K=3 and D=3.
# You must be able to pŕovide a value for K when the Estimation process 
# is beginning
#-------
P0 <- rep(1/K,K) # Initial probabilities for the first state in the hidden chain
N=30  #Binomial Sample Size
K=3   #Number of hidden states
D=3   #number of covariates in the model including the intercept
param_number <- K*D*(K-1) + K # Calculates the number of parameters of the model. This value will be used frequently in various calculations.
Beta_num <- K*D*(K-1)
#----

# Now we will assign the data to our variables 
#-------
Y<-NULL #Initialize the sequence of observable values here
X1<-NULL #Covariate # 1, Temperature
X2<-NULL #Covariate # 2, Relative Humidity
X<-NULL #Matrix of input covariates

tempVar<-df$Temperature 
tempVar2<-df$RelativeHumidity
Y<-df$RainyDays 
X1<-(tempVar-mean(tempVar))/sd(tempVar) # All variables will be normalized for ease of calculations
X2<-(tempVar2-mean(tempVar2))/sd(tempVar2) # All variables will be normalized for ease of calculation
X<-cbind(1,X1,X2)
T=length(Y) 
Nt=rep(N,T) # número de ensaios de Bernoulli associado a dada uma das T variáveis Binomiais. Cada coloquei tudo igual mas eles podem diferentes.
#----

# Initialize MCMC settings for RJAGS
#------
M=1000 # Number of samples to be drawn by the sampler
thin=4  # Thinning interval in which samples will be taken and stored
burn_in=0.10*M # Amount of values to be discarded at the beginning of the MCMC chain.
Sample_Size = M/thin # Final sample size for the MCMC chain which will be used for inference.
#----


#########################################
#####      Estimation Algorithm     #####
#########################################

# Initialize some structures to store estimation values and relabeling values in each Replica of the simulation
#--------
S_post<-NULL # Used to store the simulated posterior hidden sequence 
S_Posterior_Probs<-matrix(nrow=K,ncol=T) # Used to store the posterior probabilities of each of the hidden states
samples <- NULL# variable used to store raw mcmc output
samples2 <- NULL# variable used to store raw mcmc output for the Beta coefficients
S_Relabel <- matrix(nrow=Sample_Size,ncol=T) # this matrix will be used to store the relabeled S sequences
samples_matrix<-NULL # This matrix accepts parsed MCMC outputs.
samples_matrix2<-NULL # This matrix accepts parsed MCMC outputs.
S_matrix <- NULL #This matrix will store predictions for the S sequence in each iteration
paramlist<-NULL # Used to store a list of strings with parameter names that will generated and passed to the MCMC sampler
theta_matrix <- NULL # This matrix will store estimations for the success probabilities for each hidden state.
mcmc.pars <- array(data = NA, dim = c(Sample_Size, K, T)) # Create a temporary array to use in the permute.mcmc function. This reorders the theta matrix. 
Beta_FinalNames = array("",c(K,D,K))
Names_Beta_Final=NULL
Final.Betas=NULL
SD.Betas=NULL
CI.Betas=NULL


Betas_Post <- array(0,c(K,D,K)) # array used to store the posterior means for the Beta coefficients
Betas_SD <- array(0,c(K,D,K))
Betas_CI.Length<-array(0,c(K,D,K))

Thetas_post<-NULL # Used to store the posterior means for the success probabilities for each hidden state
Thetas_SD <- NULL
Thetas_CI.Length<-NULL

InfoCriteria_Table<-NULL

#----

# Phase 1 of MCMC sampling using RJAGS - Sample S sequence
#--------
modeloNHMM <- jags.model("Auxiliar/NHMM.bug", data = list('x' = X, 'y' = Y, 'T' = T, 'N'= N, 'K'= K, 'D' = D), n.chains = 1, n.adapt = burn_in) # A função para gerar o modelo que o RJAGS usara no amostrador
samples<- coda.samples(modeloNHMM, c('S'), thin = thin, n.iter = M)
#----

# Manipulate the output from the MCMC sampler
#-------
# Here we create some structures that will be used to manipulate 
# the ouput of the mcmc sampler. The parameter families need to be separated 
# and relabeled individually, e.g. he S 
# sequence must be relabeled. In this first phase of sampling, the S sequences are 
# sampled, and the relabeled and the posterior S sequence is calculated. It is then 
# used as input in a second phase of sampling to estimate the transition parameters and the Thetas

samples_matrix = as.matrix(samples, iters = FALSE) # the raw output from the mcmc sampler is converted to a matrix
S_matrix = samples_matrix[,1:T] # extract the values of the S sequence for all iterations of the sampler

zpvt_S = S_matrix[1,] # Set the first posterior S sequence as pivot for the relabelling algorithm
perms_S = ecr(zpivot = zpvt_S, z = S_matrix, K = K) # Run the ECR algorithm to find the permutations for relabeling

# This loop relabels the S sequence according to the permutations
# created by the ECR algorithm
for (i in 1:Sample_Size) {
  for (j in 1:T) {
    if(S_matrix[i,j]!=S_matrix[1,j]){
      S_Relabel[i,j] = which(perms_S$permutations[i,] %in% S_matrix[i,j])
    } else {
      S_Relabel[i,j] = S_matrix[i,j]
    } 
  }
}
#-----

#Estimate the posterior S sequence
#------
# the following calculates the posterior probabilities
# For each hidden state at each position in the hidden markov chain
for (i in 1:T) {
  for (j in 1:K) {
    S_Posterior_Probs[j,i]=sum(S_Relabel[,i]==j)/Sample_Size 
  }
}

for (i in 1:T) {
  S_post[i]<-rDiscrete(S_Posterior_Probs[,i])
}
#----


# Phase 2 of MCMC Sampling - Use Posterior S Sequence as Data to sample Thetas and the Beta transition coefficients
#-------
# Here we will use the estimates the predicted posterior S sequence 
# as data along with Y and X to sample the transition coefficients, 
# without having to deal with the label switching problem in this phase. 

modeloNHMM2 <- jags.model('Auxiliar/NHMM.bug', data = list('x' = X, 'y' = Y, 'S'=S_post, 'T' = T, 'N'= N, 'K' = K, 'D'= D), n.chains = 1, n.adapt = burn_in) # Run the model

#Create lists of names for the Betas and Thetas
NamesBetas <- Beta.Names(K, D) #Create a list of names of the Beta parameters to be sampled
for (i in 1:K) {
  paramlist[i] <- paste("theta[", toString(i),"]", sep = "") # We create a list of names of Thetas as characters to pass to the sampler as parameter names
}
samples2 <- coda.samples(modeloNHMM2, c(NamesBetas,paramlist), thin = thin, n.iter = M) # Sample values for the transition coefficients

#Calculate parameter estimates for Beta coefficients, store
#-------
samples_matrix2<-as.matrix(samples2, iters = FALSE)
mcmc.pars[,,1] = samples_matrix2[,(Beta_num+1):(Beta_num+K)] # extract the values of the theta vector for all iterations fo the sample
theta_matrix <- permute.mcmc(mcmc.pars,perms_S$permutations) # store the reordered thetas in a matrix in order to calculate posterior means in each Replica.


# Create the list of names for the Beta coefficients in the correct order
for (i in 1:K){
  for (d in 1:D){
    for (j in 2:K){
      Beta_FinalNames[j,d,i] <- paste("Beta",toString(i),toString(j),toString(d), sep = "")
    }
  }
}
for (i in 1:K) {
  paramlist[i] <- paste("Theta", toString(i), sep = "") # We create a list of thetas as character to pass as parameter names to the sampler
}



for (i in 1:K){
  for (j in 2:K){
    for (d in 1:D){
      Betas_Post[j,d,i] <- mean(samples_matrix2[,paste("Beta[",toString(j),",", toString(d),",",toString(i),"]", sep = "")])
      Betas_SD[j,d,i] <- sd(samples_matrix2[,paste("Beta[",toString(j),",", toString(d),",",toString(i),"]", sep = "")])
      Betas_CI.Length[j,d,i] <- abs(quantile(samples_matrix2[,paste("Beta[",toString(j),",", toString(d),",",toString(i),"]", sep = "")], probs = c(0.025))-quantile(samples_matrix2[,paste("Beta[",toString(j),",", toString(d),",",toString(i),"]", sep = "")], probs = c(0.975)))
    }
  }
}

#store those quantities in variables to create a dataframe
for (i in 1:K){
  for (j in 2:K){
    for (d in 1:D){
      Names_Beta_Final <- rbind(Names_Beta_Final, Beta_FinalNames[j,d,i])
      Final.Betas <- rbind(Final.Betas, Betas_Post[j,d,i])
      SD.Betas <- rbind(SD.Betas, Betas_SD[j,d,i])
      CI.Betas <- rbind(CI.Betas, Betas_CI.Length[j,d,i])
    }
  }
}
for (i in 1:K) {
  Thetas_post[i] <- mean(theta_matrix$output[,i,1])
  Thetas_SD[i] <- sd(theta_matrix$output[,i,1])
  Thetas_CI.Length[i] <- abs(quantile(theta_matrix$output[,i,1], probs = c(0.025))-quantile(theta_matrix$output[,i,1], probs = c(0.975)))
}

#----


# Calculate information criteria for model selection
#------

acumLL1=0 # Will be used to calculate the first part of the LL
acumLL2=0 # Will be used to calculate the second part of the LL
acumLL3=0 # Will be used to calculate the third part of the LL
LL=0 # will be used to calculate the LL
AICc=0 
BIC=0

for (i in 1:T) {
  acumLL1 = acumLL1 + Y[i]*log(Thetas_post[S_post[i]])
}
for (i in 1:T) {
  acumLL2 = acumLL2 + (Nt[i]-Y[i])*(log(1-Thetas_post[S_post[i]], base = exp(1)))
}
temp=NULL
for (i in 2:T) {
  for (g in 1:K) {
    temp[g]<-exp(X[i,]%*%matrix(Betas_Post[g,,S_post[i-1]],ncol=1))
  }
  acumLL3 = acumLL3 + (X[i,]%*%matrix(Betas_Post[g,,S_post[i-1]]) - log(sum(temp), base = exp(1)))
}
LL <- sum(log(choose(Nt,Y), base = exp(1))) + log(P0[S_post[1]]) + acumLL1 + acumLL2 + acumLL3 #Calculates the LL 

InfoCriteria_Table[1] <- -2*LL + param_number*log(T) # Calculation of the BIC for each replica 
InfoCriteria_Table[2] <- -2*(LL) + 2*param_number + (2*param_number^2 + 2*param_number)/(T - param_number - 1) # Calculation of the corrected AIC


# Create a dataframes with all the calculated values to export to pdf
#-----
Parameter<-c(Names_Beta_Final, paramlist)
Estimate<-c(Final.Betas, Thetas_post)
SD<-c(SD.Betas,Thetas_SD)
MeanCI_Length<-c(CI.Betas,Thetas_CI.Length)

df1<-data.frame(Parameter,Estimate,SD,MeanCI_Length)

FinalTime<-proc.time()
TotalTime<-(FinalTime - InitialTime)/60
Indicator <- c("BIC", "AICc", "Processing Time (mins)")
Value <- c(InfoCriteria_Table[1],InfoCriteria_Table[2], TotalTime[1])
df2 <- data.frame(Indicator, Value)

post<-as.array(samples2)
special<-list()
for (i in 1:Beta_num) {
  b<-mcmc_hist(post, pars = c(NamesBetas[i])) + xlab(Names_Beta_Final[i])
  c<-mcmc_acf(post, pars = c(NamesBetas[i])) + ggtitle(Names_Beta_Final[i])
  d<-mcmc_trace(post, pars = c(NamesBetas[i])) + xlab("Iteration") + ylab("Value")
  
  b1<-as.grob(b)
  c1<-as.grob(c)
  d1<-as.grob(d)
  special[[i]] <- list(b,c,d)
}

special 

View(df1)
View(df2)
