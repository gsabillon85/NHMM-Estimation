#####
#This file holds code used to perform a simulation study with controlled parameters for
#the non-homogeneous hidden Markov model. All parameter
####


installed.packages("gridbase")
source("Auxiliar/AuxiliarFunctions.R")
library('rjags')
library('coda')
library('label.switching')
library(e1071)
library(gridExtra)
library(gridBase)
library(grid)
library(ggplot2)
InitialTime <- proc.time()




# Initialize values to generate data
#------
### INITIALIZATION ###

#Initialization of all the values needed to simulate data
#for a non-homogeneous hidden Markov model

T=450   #Set the length of NHMM to be generated
K=3     #Set the amount of hidden states that the NHMM will contain   
D=3     #Set the amount of covariates to be used to calculate transition probabilities

# Generate Thetas at equally spaced intervals in the (0,1) parameter space
PRESET_THETAS = TRUE # Set to true if you want Thetas to be automatically preset with values.
  if(PRESET_THETAS){
    theta <- PresetThetas(K) # Initialize the vector of sucess probabilities with generated presets. 
  } else {
    theta = c(0.10,0.35,0.60,0.85)  #Manually Initialize a vector of success probabilities, each related to a hidden state.
  }

P0 <- rep(1/K,K) # Initial probabilities for the first state in the hidden chain
N=500   # Set a sample size for each Binomial observation at every instant of time 
Nt=rep(N,T) # Create a vector with T sample sizes of size N
param_number <- K*D*(K-1) + K # Calculates the number of parameters of the model. This value will be used frequently in various calculations.
R <- 3  # Amount of replicas of the simulation to be run
#----

# Initialize MCMC settings for RJAGS
#------
M=1000 # Number of samples to be drawn by the sampler
thin=4  # Thinning interval in which samples will be taken and stored
burn_in=0.10*M # Amount of values to be discarded at the beginning of the MCMC chain.
Sample_Size = M/thin # Final sample size for the MCMC chain which will be used for inference.
#----

#Create directory and filenames for outputs
#------
Directory = paste("Bayesian Simulation/Output/Simulation-",Sys.Date(),"_K",toString(K),"_D",toString(D),"_T",toString(T), sep = "")
ifelse(!dir.exists(Directory), dir.create(Directory), FALSE)
MCMC_Dir =paste(Directory, "MCMC-Objects")
ifelse(!dir.exists(MCMC_Dir), dir.create(MCMC_Dir), FALSE)
name <- paste("BayesianOutput-",Sys.Date(),"_K",toString(K),"_D",toString(D),"_T",toString(T),".pdf", sep = "") # Create the filename
filename <- paste(Directory,"/",name,sep = "")

#----
  
# Generate Beta - Transition coefficients
#------ 
### GENERATE BETAS ### 

# This section of code generates an array of transition parameters
# This array will have K*K positions. Each position will be a
# vector of D coefficients, which will be used to calculate transition probabilities

# Two options are available: Generate the Beta transition coefficients randomly
# Or set their values manually. Both cases contemplate the mlogit link function
# in which transition coefficients for the first transition for all states are
# fixed and set equal to 0.

# NOTE:setting the values for the Betas manually allows you to control
# the dynamics of the hidden Markov Chain. A graphing application
# is recommended to visualize how different values of the transition
# coefficients affect the transition probabilities that are generated. 
RandomBetas = TRUE

if (RandomBetas){ #Option 1: Generates the Betas randomly, RandomBetas must be set to TRUE.
  Betas=array(round(runif((K*D*K),-5,5),1),c(K,D,K))
  for (i in 1:K) {
    for (j in 1:K) {
      for (d in 1:D) {
        if (j == 1) {
          Betas[j,d,i] = 0
        }
      }
    }
  }  
} else { # Option 2: Assign values to the Betas manually. 
  # NOTE: You must add or remove array positions to match your 
  # number of hidden states and number of covariates. 
  # The structure provided here contemplates an NHMM with
  # 2 covariates and an intercept, along with 3 hidden states. 
  
  Betas[1,1,1]=0
  Betas[1,2,1]=0
  Betas[1,3,1]=0
  Betas[2,1,1]=-1.3
  Betas[2,2,1]=-1.4
  Betas[2,3,1]=-1.1
  Betas[3,1,1]=-1.5
  Betas[3,2,1]=-1.6
  Betas[3,3,1]=-1.2
  
  Betas[1,1,2]=0
  Betas[1,2,2]=0
  Betas[1,3,2]=0
  Betas[2,1,2]=1.4
  Betas[2,2,2]=1.1
  Betas[2,3,2]=1.2
  Betas[3,1,2]=-1.1
  Betas[3,2,2]=-1.2
  Betas[3,3,2]=-1.3
  
  Betas[1,1,3]=0
  Betas[1,2,3]=0
  Betas[1,3,3]=0
  Betas[2,1,3]=-1.4
  Betas[2,2,3]=-1.2
  Betas[2,3,3]=-1.3
  Betas[3,1,3]=1.3
  Betas[3,2,3]=1.5
  Betas[3,3,3]=1.2
  
}
#----

# Initialize some important structures for storage of generated information
#-------
# It would be convenient to store the values for Standard Deviaton, Bias
# Mean Square Error so that we could later access them, in case we would
# want to analyze an individual replica. The following structures will be
# used to store: Y, X, S, S_post, Standard Deviation, Bias, MSE and size of 
# credibility interval for all parameters, as well as Information Criteria 
# for each of the replicas.

Y_table <- matrix(nrow = T, ncol = R)
X_Table <- array(0,c(T,D,R))
S_table <- matrix(nrow = T, ncol = R)
S_Post_Probs_table <- array(0, c(T,K,R))
S_Post_Table <- matrix(nrow = T, ncol = R)
success <- rep(0,R) 

Betas_Post <- array(0,c(K,D,K,R)) # array used to store the posterior means for the Beta coefficients
Betas_SD_Table <- array(0,c(K,D,K,R))
Betas_Bias_Table<-array(0,c(K,D,K,R))
Betas_MSE_Table<-array(0,c(K,D,K,R))
Betas_CI.Length_Table<-array(0,c(K,D,K,R))
Betas_Coverage_Count <- array(0,c(K,D,K)) 

Thetas_post<-matrix(nrow = K, ncol = R) # Used to store the posterior means for the success probabilities for each hidden state
Thetas_SD_Table <- matrix(nrow = K, ncol = R)
Thetas_Bias_Table<-matrix(nrow = K, ncol = R)
Thetas_MSE_Table<-matrix(nrow = K, ncol = R)
Thetas_CI.Length_Table<-matrix(nrow = K, ncol = R)
Thetas_Coverage_Count <- c(rep(0,K)) 

InfoCriteria_Table<-matrix(nrow = 2, ncol = R)
#----


 for (p in 1:R) {
  #Generate covariates
  #------
  # For the model there will be D covariates
  X <- matrix(0, nrow = T, ncol = D)
  for (d in 1:D) {
    if (d == 1){
      X[,d] = rep(1,T)
    }else{
      X[,d]<- X1<-rnorm(T,mean = 0, sd = 1)
    }
  }
  #----
  
  #Generate Data
  #--------
  # Input the values for the reviously initialized variables
  # into the DataGenerator function to generate data in the
  # structure of a NHMM
  data <- DataGeneratorNHMM(X, Betas, K, T, theta, Nt)
  Y <- data[,1] # The observations generated 
  S <- data[,2] # The hidden state sequence
  #----
  
  #########################################
  #####      Estimation Algorithm     #####
  #########################################
  
  # Initialize some structures to store estimation values and relabeling values in each Replica of the simulation
  #--------
  S_post<-NULL # Used to store the simulated posterior hidden sequence 
  S_Posterior_Probs<-matrix(nrow=K,ncol=T) # Used to store the posterior probabilities of each of the hidden states
  Beta_Post_Array=array(round(runif((K*D*K),0,1),1),c(K,D,K)) # This array is used to store estimations for the transition coefficients in a structure which is convenient to calculate LL and information critertia
  samples <- NULL# variable used to store raw mcmc output
  samples2 <- NULL# variable used to store raw mcmc output for the Beta coefficients
  paramlist<-NULL # Used to store a list of strings with parameter names that will generated and passed to the MCMC sampler
  samples_matrix<-NULL # This matrix accepts parsed MCMC outputs.
  samples_matrix2<-NULL # This matrix accepts parsed MCMC outputs.
  S_matrix <- NULL #This matrix will store predictions for the S sequence in each iteration
  theta_matrix <- NULL # This matrix will store estimations for the success probabilities for each hidden state.
  S_Relabel <- matrix(nrow=Sample_Size,ncol=T) # this matrix will be used to store the relabeled S sequences
  #----
  
  # Phase 1 of MCMC sampling using RJAGS - Sample S sequence and success probabilities
  #--------
  print(paste("Replica number:",toString(p), collapse = "")) # A line of code to visualize the number of the current replica
  modeloNHMM <- jags.model("Auxiliar/NHMM.bug", data = list('x' = X, 'y' = Y, 'T' = T, 'N'= N, 'K'= K, 'D' = D), n.chains = 1, n.adapt = burn_in) # A função para gerar o modelo que o RJAGS usara no amostrador
  for (i in 1:K) {
    paramlist[i] <- paste("theta[", toString(i),"]", sep = "") # Criamos uma lista dos thetas como character para passar pro amostrador
  }
  samples <- coda.samples(modeloNHMM, c(paramlist,'S'), thin = thin, n.iter = M) # the raw mcmc output is stored into an object of type .mcmc
  save(samples, file=paste(MCMC_Dir,"/A_Replica",toString(p),".rda", sep = "" ))
  
  #----
  
  # Manipulate the output from the MCMC sampler
  #-------
  # Here we create some structures that will be used to manipulate 
  # the ouput of the mcmc sampler. The parameter families need to be separated 
  # and relabeled individually, e.g. the K thetas must be relabeled, the S 
  # sequence must be relabeled. In this first phase of sampling, the S sequences are 
  # sampled, and the relabeled and the posterior S sequence is calculated. It is then 
  # used as input in a second phase of sampling to estimate the transition parameters
  
  samples_matrix = as.matrix(samples, iters = FALSE) # the raw output from the mcmc sampler is converted to a matrix
  S_matrix = samples_matrix[,1:T] # extract the values of the S sequence for all iterations of the sampler
  mcmc.pars <- array(data = NA, dim = c(Sample_Size, K, T)) # Create a temporary array to use in the permute.mcmc function. This reorders the theta matrix. 
  mcmc.pars[,,1] = samples_matrix[,(T+1):(T+K)] # extract the values of the theta vector for all iterations fo the sample
 
  zpvt_S = S # Set the original S sequence as pivot for the relabelling algorithm
  perms_S = ecr(zpivot = zpvt_S, z = S_matrix, K = K) # Run the ECR algorithm to find the permutations for relabeling
  theta_matrix <- permute.mcmc(mcmc.pars,perms_S$permutations) # store the reordered thetas in a matrix in order to calculate posterior means in each Replica.
  
  for (i in 1:Sample_Size) {
    for (j in 1:T) {
      if(S_matrix[i,j]!=S[j]){
        S_Relabel[i,j] = which(perms_S$permutations[i,] %in% S_matrix[i,j])
        j=j+1
      } else {
        S_Relabel[i,j] = S_matrix[i,j]
      } 
    }
  }
  table(S_Relabel[5,], S)
  #----
  
  # Estimation of Posterior S Sequence and Posterior theta vector
  #-------
  
  # the following for calculates the posterior probabilities
  # For each hidden state at each position in the hidden markov chain
  for (i in 1:T) {
    for (j in 1:K) {
      S_Posterior_Probs[j,i]=sum(S_Relabel[,i]==j)/Sample_Size 
    }
  }
  
  # Using the posterior probabilities, we will randomly generate a posterior 
  # S sequence, which will then be used as data to sample the transition parameters
  for (i in 1:T) {
    S_post[i]<-rDiscrete(S_Posterior_Probs[,i])
    if (S[i]==S_post[i]){
      success[p] = success[p] + 1 # Here we count the successes in estimating the value of the hidden state at each position in the hidden chain. 
    }
  }
  
  # And we calculate the posterior means for the theta vector
  for (i in 1:K) {
    Thetas_post[i,p]<-mean(theta_matrix$output[,i,1])
  }
  #----
  
  # Phase 2 of MCMC Sampling - Use Posterior S Sequence and 
  # Posterior Theta as Data to sample the Beta transition coefficients
  #-------
  # Here we will use the estimates for the success parameters(theta) and
  # the predicted posterior S sequence as data along with Y and X to sample
  # the transition coefficients, without having to deal with the label switching 
  # problem in this phase. 
  
  modeloNHMM2 <- jags.model('Auxiliar/NHMM.bug', data = list('x' = X, 'y' = Y, 'S'=S_post, 'theta'= Thetas_post[,p], 'T' = T, 'N'= N, 'K' = K, 'D'= D), n.chains = 1, n.adapt = burn_in) # Run the model
  NamesBetas <- Beta.Names(K, D) #Create a list of names of the Beta parameters to be sampled
  samples2 <- coda.samples(modeloNHMM2, NamesBetas, thin = thin, n.iter = M) # Sample values for the transition coefficients
  save(samples2, file=paste(MCMC_Dir,"/B_Replica",toString(p),".rda", sep = "" ))
  #----
  
  #Calculate parameter estimates for Beta coefficients, store
  #-------
  samples_matrix2<-as.matrix(samples2, iters = FALSE)
  
  for (i in 1:K){
    for (j in 2:K){
      for (d in 1:D){
        Betas_Post[j,d,i,p] <- mean(samples_matrix2[,paste("Beta[",toString(j),",", toString(d),",",toString(i),"]", sep = "")])
        Betas_SD_Table[j,d,i,p] <- sd(samples_matrix2[,paste("Beta[",toString(j),",", toString(d),",",toString(i),"]", sep = "")])
        Betas_Bias_Table[j,d,i,p] <- Betas[j,d,i]-Betas_Post[j,d,i,p]
        Betas_MSE_Table[j,d,i,p] <- (Betas_SD_Table[j,d,i,p])^2 + (Betas_Bias_Table[j,d,i,p])^2
        Betas_CI.Length_Table[j,d,i,p] <- abs(quantile(samples_matrix2[,paste("Beta[",toString(j),",", toString(d),",",toString(i),"]", sep = "")], probs = c(0.025))-quantile(samples_matrix2[,paste("Beta[",toString(j),",", toString(d),",",toString(i),"]", sep = "")], probs = c(0.975)))
        if(data.table::between(Betas[j,d,i], quantile(samples_matrix2[,paste("Beta[",toString(j),",", toString(d),",",toString(i),"]", sep = "")],probs = 0.025),quantile(samples_matrix2[,paste("Beta[",toString(j),",", toString(d),",",toString(i),"]", sep = "")],probs = 0.975))){Betas_Coverage_Count[j,d,i]=Betas_Coverage_Count[j,d,i]+1}
      }
    }
  }
  for (i in 1:K) {
    Thetas_SD_Table[i,p] <- sd(theta_matrix$output[,i,1])
    Thetas_Bias_Table[i,p] <- theta[i] - Thetas_post[i,p] 
    Thetas_MSE_Table[i,p] <- (Thetas_SD_Table[i,p])^2 + (Thetas_Bias_Table[i,p])^2
    Thetas_CI.Length_Table[i,p] <- abs(quantile(theta_matrix$output[,i,1], probs = c(0.025))-quantile(theta_matrix$output[,i,1], probs = c(0.975)))
    if(data.table::between(theta[i], quantile(theta_matrix$output[,i,1],probs = 0.025),quantile(theta_matrix$output[,i,1],probs = 0.975))){Thetas_Coverage_Count[i]=Thetas_Coverage_Count[i]+1}
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
    acumLL1 = acumLL1 + Y[i]*log(Thetas_post[S_post[i],p])
  }
  for (i in 1:T) {
    acumLL2 = acumLL2 + (Nt[i]-Y[i])*(log(1-Thetas_post[S_post[i],p], base = exp(1)))
  }
  temp=NULL
  for (i in 2:T) {
    for (g in 1:K) {
      temp[g]<-exp(X[i,]%*%matrix(Betas_Post[g,,S_post[i-1],p],ncol=1))
    }
    acumLL3 = acumLL3 + (X[i,]%*%matrix(Betas_Post[g,,S_post[i-1],p]) - log(sum(temp), base = exp(1)))
  }
  LL <- sum(log(choose(Nt,Y), base = exp(1))) + log(P0[S_post[1]]) + acumLL1 + acumLL2 + acumLL3 #Calculates the LL 
  
  InfoCriteria_Table[1,p] <- -2*LL + param_number*log(T) # Calculation of the BIC for each replica 
  InfoCriteria_Table[2,p] <- -2*(LL) + 2*param_number + (2*param_number^2 + 2*param_number)/(T - param_number - 1) # Calculation of the corrected AIC
  
  Y_table[,p] <- Y
  S_table[,p] <- S
  X_Table[,,p] <- X
  S_Post_Probs_table[,,p] <- t(S_Posterior_Probs) # Store the posterior probabilities for S at each Replica in a table
  S_Post_Table[,p] <- S_post
  #----
  
} ########## End of Simulation Replicas#########
 
  options(digits=4)
  options(scipen=999)
  
  #Define some arrays and vectors to store final calculations for estimates
  #------
  Betas_Final = array(0,c(K,D,K))
  SD_Betas = array(0,c(K,D,K))
  Bias_Betas = array(0,c(K,D,K))
  MSE_Betas = array(0,c(K,D,K))
  Beta_FinalNames = array("",c(K,D,K))
  Betas_Cover_Probs = array(0,c(K,D,K))
  Mean_CI_Betas = array(0,c(K,D,K))
  Names_Beta_Final=NULL
  Final.Betas=NULL
  SD.Betas=NULL
  Bias.Betas=NULL
  MSE.Betas=NULL
  Cover.Probs.Betas=NULL
  Mean.CI.Betas=NULL
  Real.Betas=NULL
  
  Thetas_Final=NULL
  SD_Thetas=NULL
  Bias_Thetas=NULL
  MSE_Thetas=NULL
  Cover.Probs.Thetas=NULL
  Mean.CI.Thetas=NULL
  
  Success_Rate=NULL
  AICc_Final=NULL
  BIC_Final=NULL
  #----
  
  
  # Calculate final values for each of the parameters related to the model
  # and for the performance indicators
  Success_Rate <- success/T # Calculate the success rate in predicting S in every replica
  
  Betas_Cover_Probs <- Betas_Coverage_Count/R
  Cover.Probs.Thetas <- Thetas_Coverage_Count/R
  
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
  
  #Calculate parameter values and performance indicators for the Beta coefficients
  for (i in 1:K){
    for (j in 2:K){
      for (d in 1:D){
        Betas_Final[j,d,i] <- mean(Betas_Post[j,d,i,])
        SD_Betas[j,d,i] <- sd(Betas_Post[j,d,i,])
        Bias_Betas[j,d,i] <- Betas[j,d,i] - Betas_Final[j,d,i]
        MSE_Betas[j,d,i] <- (Bias_Betas[j,d,i])^2 + (SD_Betas[j,d,i])^2
        Mean_CI_Betas[j,d,i] <- mean(Betas_CI.Length_Table[j,d,i,])
      }
    }
  }
  
  #Calculate parameter values and performance indicators for the Theta success probabilities
  for (i in 1:K) {
    Thetas_Final[i] <- mean(Thetas_post[i,])
    SD_Thetas[i] <- sd(Thetas_post[i,])
    Bias_Thetas[i] <- Thetas_Final[i] - theta[i]
    MSE_Thetas[i] <- (SD_Thetas[i])^2 + (Bias_Thetas[i])^2
    Mean.CI.Thetas[i] <- mean(Thetas_CI.Length_Table[i,])
  }
  
  #store those quantities in variables to create a dataframe
  for (i in 1:K){
    for (j in 2:K){
      for (d in 1:D){
        Names_Beta_Final <- rbind(Names_Beta_Final, Beta_FinalNames[j,d,i])
        Final.Betas <- rbind(Final.Betas, Betas_Final[j,d,i])
        SD.Betas <- rbind(SD.Betas, SD_Betas[j,d,i])
        Bias.Betas <- rbind(Bias.Betas, Bias_Betas[j,d,i])
        MSE.Betas <- rbind(MSE.Betas, MSE_Betas[j,d,i])
        Mean.CI.Betas <- rbind(Mean.CI.Betas ,Mean_CI_Betas[j,d,i])
        Cover.Probs.Betas <- rbind(Cover.Probs.Betas, Betas_Cover_Probs[j,d,i])
        Real.Betas <- rbind(Real.Betas, Betas[j,d,i])
      }
    }
  }
  
  # Create a dataframes with all the calculated values to export to pdf
  #-----
  Parameter<-c(Names_Beta_Final, paramlist)
  Real<-c(Real.Betas,theta)
  Estimate<-c(Final.Betas, Thetas_Final)
  SD<-c(SD.Betas,SD_Thetas)
  Bias<-c(Bias.Betas,Bias_Thetas)
  MSE<-c(MSE.Betas, MSE_Thetas)
  Coverage<-c(Cover.Probs.Betas,Cover.Probs.Thetas)
  MeanCI_Length<-c(Mean.CI.Betas, Mean.CI.Thetas)

  DF1<-data.frame(Parameter,Real, Estimate, SD, Bias,MSE,Coverage,MeanCI_Length)
  
  MeanSuccessRate <- mean(Success_Rate)
  BIC_Final = mean(InfoCriteria_Table[1,])
  AICc_Final = mean(InfoCriteria_Table[2,])
  FinalTime<-proc.time()
  TotalTime<-(FinalTime - InitialTime)/60
  
  Indicator <- c("Success Rate", "Mean BIC", "Mean AICc", "Processing Time (mins)")
  Value <- c(MeanSuccessRate, BIC_Final, AICc_Final, TotalTime[1])
  DF2 <- data.frame(Indicator, Value)
  #----
  
  #Exporting the two dataframse to pdf, using the gridextra package
  #----
  saveRDS(Y_table ,paste(Directory, "/Y_Table.rds",sep = ""))
  saveRDS(X_Table,paste(Directory, "/X_Table.rds",sep = ""))
  saveRDS(S_Table,paste(Directory, "/S_Table.rds",sep = ""))
  saveRDS(S_Post_Table,paste(Directory, "/S_Post_Table.rds",sep = ""))
  
  saveRDS(Betas_MSE_Table ,paste(Directory, "/Beta_MSE_Table.rds",sep = ""))
  saveRDS(Betas_SD_Table ,paste(Directory, "/Beta_SD_Table.rds",sep = ""))
  saveRDS(Betas_Bias_Table ,paste(Directory, "/Beta_Bias_Table.rds",sep = ""))
  
  saveRDS(Thetas_MSE_Table, paste(Directory, "/Theta_MSE_Table.rds",sep = ""))
  saveRDS(Thetas_SD_Table, paste(Directory, "/Theta_SD_Table.rds",sep = ""))
  saveRDS(Thetas_Bias_Table, paste(Directory, "/Theta_Bias_Table.rds",sep = ""))
  
  saveRDS(InfoCriteria_Table ,paste(Directory, "/InfoCrit_Table.rds",sep = ""))
  
  
  #Convert the two dataframes to Grobtables
  Table1 <- tableGrob(DF1)
  Table2 <- tableGrob(DF2)
  lay <- rbind(c(1),c(2)) #Establish a layout matrix for the grobtables
  
  g<-grid.arrange(Table1,Table2) #Arrange the tables
  pdf(filename, height = 4.5+param_number/2.20, width = 8.5) #Set the format you will save to along with file dimension
  ggsave(filename, g) #save the file
 
  
  
