rm(list=ls())
library(deSolve)
library(rootSolve)
library(reshape2)
library(ggplot2)
library(ggstance)#not needed?
library(cowplot)
library(tidyverse)

###Simulation model based on Dwyer et al host-baculovirus models
###Intended to model evolution of disease resistance alleles in a real dynamical
### host-disease system.

##Disease outbreak model parameters (value can be changed in simulation function below)
Tmax = 100 # Time horizon for disease model - could be increased for burnout
TimeStep = 1 # Integration time step
Time <- seq(0, Tmax, by = TimeStep) #Vector of times (in days)
v1<-0.75 #Transmission parameter for homozygote wild type (0 copies of resistance allele)
v2<-0.5 #Transmission parameter for heterozygotes (1 copy)
v3<-0.25 # Transmission parameter for homozygotes with resistance allele
r1<-20 #Reproduction rate for homozygote wild type (0 copies of resistance allele)
r2<-15 #Reproduction rate for heterozygotes (1 copy)
r3<-10 #Reproduction rate for homozygote wild type (0 copies of resistance allele)
mu<-0.4 #Decay rate of infectious cadavers
tau<-12 #Delay from infection to death of infected individuals, in days
V<-0.86^2 #Variance in transmission rate between individuals
X0<-0.5 #Starting frequency of disease resistance allele
N0<-1 # Starting host population density, expressed in terms of outbreak threshold density (1=threshold for spread)(possibly should be adjusted to equilibrium value)
S10<-N0*(1-X0)^2 # Starting abundance of  homozygote wild type 
S20<-N0*2*X0*(1-X0)# Starting abundance of  heterozygote 
S30<-N0*X0^2# Starting abundance of  homozygote for resistance allele
P0<-.1 # Starting pathogen density (possibly should be adjusted to equilibrium value)

Zt<-P0
parms <- c(v1,v2,v3,mu,tau,V,S10,S20,S30) #parameter vector for disease outbreak ODE
parms
Y0 <- c(S10,S20,S30,Zt) # starting state value vector
Y0
gamma<-0.01 # parameter for long-term baculovirus persistence in the environment
phi<-3 #parameter for persistence of baculovirus from cadavers to the next year, scaled by threshold population size


###Disease outbreak ODE, used in simulation function
outbreak_model <- function(t, y, parms) {
  v1=parms[1]; v2=parms[2] ; v3= parms[3]; mu=parms[4]; tau=parms[5]; 
  V=parms[6]; S1 <- y[1]; S2 <- y[2]; S3 <- y[3] ; P <- y[4];
  S10 <- parms[7]; S20 <- parms[8]; S30 <- parms[9] ;
  ###Lag function - there are no cadavers produced until the incubation time has elapsed, and the cadaver production rate is a function of susceptible pop size 12 days previously
  if(t<tau){
    S1lag<-0
    S2lag<-0
    S3lag<-0
    Plag<-0
  }
  else{
    S1lag<-lagvalue(t-tau, 1)
    S2lag<-lagvalue(t-tau, 2)
    S3lag<-lagvalue(t-tau, 3)
    Plag<-lagvalue(t-tau, 4)
  }
  dS1dt <- -v1*S1*P*(((S1+S2+S3)/(S10+S20+S30))^V)    
  dS2dt <- -v2*S2*P*(((S1+S2+S3)/(S10+S20+S30))^V)
  dS3dt <- -v3*S3*P*(((S1+S2+S3)/(S10+S20+S30))^V)         
  dPdt <- v1*S1lag*Plag*(((S1lag+S2lag+S3lag)/(S10+S20+S30))^V)+v2*S2lag*Plag*(((S1lag+S2lag+S3lag)/(S10+S20+S30))^V)+v3*S3lag*Plag*(((S1lag+S2lag+S3lag)/(S10+S20+S30))^V)-mu*P    
  return(list(c( dS1dt, dS2dt, dS3dt,dPdt)));
}


ode_out <- dede(Y0, Time,outbreak_model,parms,method='lsoda')
ode_out
par(las=1,bty="l")
matplot(ode_out[,1],ode_out[,-1], type='l', xlab="time", ylab="density") 

###Function to simulation population dynamics, with single runs to see population dynamics
ecoevo_sims<-function(model=outbreak_model,vs=c(0.75,0.5,0.25),rs=c(15,10,5),mu=0.4,tau=12, V=0.86^2,X0=0.01, P0=0.1, N0=0.5,Tmax=100,Years=10,phi=3,nsims=10, popfactor=1000){
    S<-c(N0*(1-X0)^2,2*N0*X0*(1-X0),N0*X0^2)##Vector of susceptible genotype frequencies
    
    parms <- c(vs[1],vs[2],vs[3],mu,tau,V,S[1],S[2],S[3])
    parms
    
    Y0 <- c(S[1],S[2],S[3],P0)
    
    datamat<-matrix(data=NA,nrow=Years,ncol=7) ##data matrix to take results of one population simulation
    datamat[1,1:4]<-Y0
    datamat[,5]<-X0
    datamat[,6]<-seq(from=1,to=Years,by=1)
    
    TimeStep = 1 # integration time step
    Time <- seq(0, Tmax, by = TimeStep)
    
    for(i in 1:(Years-1)){
      
      Y<-datamat[i,1:4] # Extract starting values for current year
      ode_out <- dede(Y, Time,model,parms) #Run disease outbreak
      
      Nt<-ode_out[Tmax,2:4] #Extract abundances of surviving hosts by genotype after outbreak
      Zt<-ode_out[Tmax,5] #Extract abundance of pathogens at end of outbreak 
      
      Nt1prime<-rs*Nt #Gamete production by genotype, before mating
      Nt1primetotal<-sum(Nt1prime) #Total population size
      
      xtmean<-0.5*(Nt1prime[2]/Nt1primetotal)+Nt1prime[3]/Nt1primetotal  #Mean frequency of disease resistance allele
      chromosomes<-round(2*Nt1primetotal)*popfactor ##Count number of alleles, for size of binomial draw. Popfactor multiplies density because it is scaled by threshold density, not by actual number of individuals
      xt1<-rbinom(1,chromosomes,xtmean)/chromosomes ## The stochastic frequency of the disease resistance allele, drawn from a binomial
      Nt1<-c((1-xt1)^2*Nt1primetotal, 2*(1-xt1)*xt1*Nt1primetotal,xt1^2*Nt1primetotal) ##Population vector of genotypes, for the following year
      Zt1<-phi*Nt1primetotal+gamma*Zt #Pathogen density for the following year
      
      #Inputting values for the following year
      datamat[i+1,1:3]<-Nt1 
      datamat[i+1,4]<-Zt1
      datamat[i+1,5]<-xt1
    }
    out<-datamat
    return(out)
  }



s1<-ecoevo_sims(Years=300,vs=c(0.7,0.66, 0.62),rs=c(20,18,16), phi=3,X0=1/1000,P0=1, N0=2.5)
s1
outbreak1<-data.frame(N=rowSums(s1[,1:3]),Z=s1[,4],t=s1[,6],x=s1[,5])
str(outbreak1)
ggplot(outbreak1,aes(x=t,y=N))+geom_line()

outbreak2<-melt(outbreak1,id.vars='t')
str(outbreak2)
ggplot(outbreak2,aes(x=t,y=value,color=variable))+geom_line()+xlim(50,75)+scale_y_log10()

freq<-data.frame(t=s1[,6],X=s1[,5])
str(freq)
ggplot(freq,aes(x=t,y=X))+geom_line()

