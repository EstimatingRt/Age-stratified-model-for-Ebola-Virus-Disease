
###################################################
# TWO-HOST AGE-STRATIFIED MODEL VS ONE-HOST MODEL #
###################################################

# Cleanup - CAUTION clears R environment ####
rm(list=ls())

# Setup - load packages and define plotting functions ####
library("EpiEstim")
library("sparsevar")
library("ggplot2")
library("dplyr")
library("tidyverse")
library("readxl")
library("zoo")


#The one-host model R_t inference ####

one_host_model <- function(n_days,cases,w){
  
  a<-c() #shape parameter
  p<-c() #infection potential
  b<-c()#scale parameter
  alpha <- 1 #shape hyperparamter 
  beta <- 5 #scale hyperparameter
  for (i in 1:(n_days-1)){
    p[1] <- NaN
    p[i+1] <- sum(cases[i:1]*w[1:i]) #total infectivity potential on day t starting from day t=2
  }
  p<-p
  for (j in 2:(n_days-tau)){
    a[j-1] <- alpha+sum(cases[j:(j+tau)])
    b[j-1] <- 1/(sum(p[j:(j+tau)])+1/beta) 
  }
  return(data.frame(shape=a,scale=b,R_mean=a*b))
}


#The two-host model inference ####

two_host_model <- function(I1,I2,n_days,cases,phi,m){
  
  w_adults <- discr_si(m,mu,sd) # normalised discretised adult serial interval distribution
  w_adults <-w_adults/sum(w_adults)
  
  w_children <- discr_si(m,mu*phi,sd) #phi sets a difference between age-specific serial intervals
  w_children<- w_children/sum(w_children) # normalised discretised child serial interval distribution
  
  
  a <- c() #shape parameter
  b <- c() #scale parameter
  alpha <- 1 #shape hyperparameter
  beta <- 5 #scale hyperparamter 
  eta <- c()
  for (i in 1:(n_days-1)){
    eta[1] <- NaN
    eta[i+1] <- ((C_eff[1,1]+C_eff[1,2])*sum(I1[i:1]*w_children[1:i])+(C_eff[2,2]+C_eff[2,1])*sum(I2[i:1]*w_adults[1:i]))/lambda ##total infectivity potential on day t starting from day t=2
                                                                                
  }
  eta <- eta
  for (j in 2:(n_days-tau)){
    a[j-1] <- alpha+sum(cases[j:(j+tau)])
    b[j-1] <- 1/(sum(eta[j:(j+tau)])+1/beta)
  }
  
  return(data.frame(shape=a,scale=b,R_mean=a*b))
  
}

#parameterising the population serial interval distribution when age-specific serial intervals are introduced

pop_si_dist_corrected <- function(phi,m){
  
  w_adults <- discr_si(m,mu,sd) 
  w_adults <-w_adults/sum(w_adults) 
  
  w_children <- discr_si(m,mu*phi,sd) 
  w_children<- w_children/sum(w_children) 
  
  w_corrected <- ((C_eff[1,1]+C_eff[1,2])*(v11/(v11+v12))*w_children+(C_eff[2,1]+C_eff[2,2])*(v12/(v11+v12))*w_adults)/((v11/(v11+v12))*(C_eff[1,1]+C_eff[1,2])+(C_eff[2,1]+C_eff[2,2])*(v12/(v12+v11))) #normalised population serial interval distribution 
  
  return(w_corrected) #returns the corrected population serial interval distribution when serial intervals are age-specific
  
}

#Testing the models for simulated data ####

#Loading population and contact data

load("./poptotal.rdata")
load("./contact_all.rdata")
index <- which(poptotal$iso3c=="COD") #index corresponding to DRC
P <- unlist(poptotal[index,],use.names = FALSE)[4:20] #population age structure in 5-year age intervals for DRC
M_all <- contact_all$COD #extracting the standard contact matrix for DRC for all age groups at all locations

#Aggregating age groups to build a 2x2 contact matrix 

c11<- c()
c12<- c()
c22 <- c()
c21<-c()
for (i in 1:4){
  c11[i] <- P[i]*sum(M_all[i,1:4])/sum(P[1:4])
  c12[i] <- P[i]*sum(M_all[i,5:16])/sum(P[1:4])
}
for(i in 5:16){
  c21[i-4] <- P[i]*sum(M_all[i,1:4])/sum(P[5:16])
  c22[i-4] <- P[i]*sum(M_all[i,5:16])/sum(P[5:16])
}
C11 <- sum(c11)
C12 <- sum(c12)
C21 <- sum(c21)
C22 <- sum(c22)

C <- matrix(c(C11,C12,C21,C22),nrow=2,ncol=2) #standard contact matrix weighted to contain only the relevant age groups (< 20 and >=20)

p_da <- 0.7 # adult case fatality ratio for community cases 
p_dc <- 0.58 # child case fatality ratio for community cases 
p_f <- 0.571 #probability of unsafe burial given death in the community
C[1,2] <- C[1,2]*(1+5.67*p_f*p_dc) #scaling up the child-to-adult contacts
C[2,2] <- C[2,2]*(1+5.67*p_f*p_da) #scaling up the adult-to-adult contacts
C_eff <- C #effective contact matrix 
k<-1 #k sets a difference between child and adult pathogen transmissibility (see Supplementary File E)
C_eff[1,] <- k*C_eff[1,]
lambda<-spectralRadius(C_eff) # the dominant eigenvalue of the effective contact matrix 

#entries of the eigenvector corresponding to lambda (used further in the code)

v11 <- eigen(t(C_eff))$vectors[,1][1]
v12 <- eigen(t(C_eff))$vectors[,1][2]

tau <- 6 #sliding window of one week
n_days_sim <- 100 #outbreak duration 
mu <- 15.3 #literature value for the mean of the continuous EVD serial interval distribution 
sd <- 9.3 #literature value of the standard deviation of the continuous EVD serial interval distribution 

phi <- 1 #phi = 1 when mean serial is not age-specific. Otherwise, three scenarios are tested: phi=0.75,0.5, 0.25.

m_sim <- seq(1,n_days_sim,by=1)
w_lit_sim <- discr_si(m_sim,mu,sd)
w_lit_sim <- w_lit_sim/sum(w_lit_sim) # normalised discretised serial interval distribution that assumes literature values for mean and standard deviation


w_corrected_sim <- pop_si_dist_corrected(phi,m=m_sim) #phi ≠ 1

#Setting underlying values of R_t

FUN<- function(x){
  return(-tanh(0.2*x-10)+1.5)
}

x<-seq(1,n_days_sim,by=1)
R_prior<- unlist(lapply(x,FUN)) 


#Daily incidence data by age simulated by the two-host model

I1_simu <- c()
I2_simu <- c()
B<-1000
incidence_simulation <- replicate(B, { 
  
  w_adults <- w_lit_sim
  w_children <- discr_si(m_sim,mu*phi,sd)
  w_children<- w_children/sum(w_children) 
  
  for (i in 1:(n_days_sim-1)){
    I1_simu[1] <- 10
    I2_simu[1] <- 10
    value1 <- rpois(1,(R_prior[i]/lambda)*(C_eff[1,1]*sum(I1_simu[i:1]*w_children[1:i])+C_eff[2,1]*sum(I2_simu[i:1]*w_adults[1:i])))
    value2 <- rpois(1,(R_prior[i]/lambda)*(C_eff[2,2]*sum(I2_simu[i:1]*w_adults[1:i])+C_eff[1,2]*sum(I1_simu[i:1]*w_children[1:i])))
    I1_simu <- c(I1_simu,value1)
    I2_simu <- c(I2_simu,value2)
  }
  
  return(list(I1_simu,I2_simu))
  
  
}
)

I1_simulated <- round(colMeans(do.call(rbind, incidence_simulation[1,])),digits=0) #child cases
I2_simulated <- round(colMeans(do.call(rbind, incidence_simulation[2,])),digits=0) #adult cases
cases_simulated <- I1_simulated+I2_simulated #total cases

#plotting daily incidence counts by age 

Children <- I1_simulated[-1]
Adults <- I2_simulated[-1]
df_incidence <- rbind(Children,Adults)
barplot(df_incidence,legend.text=TRUE,ylim=c(0,45)) # Figure 1a) in the report (k=1,phi=1)

estimate_one_host_sim <- one_host_model(n_days=n_days_sim,cases=cases_simulated,w=w_lit_sim)
estimate_one_host_sim_corrected <- one_host_model(n_days=n_days_sim,cases=cases_simulated,w=w_corrected_sim)
estimate_two_host_sim <- two_host_model(I1=I1_simulated,I2=I2_simulated,n_days=n_days_sim, cases=cases_simulated,phi,m=m_sim)

R_prior <-R_prior[(tau+2):length(cases_simulated)]
t<-seq(tau+2,length.out=length(estimate_one_host_sim$R_mean),by=1) #R_t is estimated from day tau+2

df_true_value <- data.frame(day=t,R_t=R_prior)
df_one_host_sim <- data.frame(day=t,R_t=estimate_one_host_sim$R_mean)
df_one_host_sim_corrected <- data.frame(day=t,R_t=estimate_one_host_sim_corrected$R_mean)
df_two_host_sim<-data.frame(day=t,R_t=estimate_two_host_sim$R_mean)
df_two_host_sim_children <- data.frame(day=t,R_t=estimate_two_host_sim$R_mean*(C_eff[1,1]+C_eff[1,2])/lambda)
df_two_host_sim_adults<-data.frame(day=t,R_t=estimate_two_host_sim$R_mean*(C_eff[2,1]+C_eff[2,2])/lambda)

upper95_one_host_sim <- qgamma(0.975,shape=estimate_one_host_sim$shape,scale=estimate_one_host_sim$scale) #lower bounds for the 95% intervals for the one-host model R_t estimates
lower95_one_host_sim<-qgamma(0.025,shape=estimate_one_host_sim$shape,scale=estimate_one_host_sim$scale) #upper bounds for the 95% intervals for the one-host model R_t estimates
upper95_one_host_sim_corrected <- qgamma(0.975,shape=estimate_one_host_sim_corrected$shape,scale=estimate_one_host_sim_corrected$scale) #the population serial interval distribution is corrected for the one-host model (phi ≠ 1)
lower95_one_host_sim_corrected<-qgamma(0.025,shape=estimate_one_host_sim_corrected$shape,scale=estimate_one_host_sim_corrected$scale)
upper95_two_host_sim<- qgamma(0.975,shape=estimate_two_host_sim$shape,scale=estimate_two_host_sim$scale)#lower bounds for the 95% intervals for the two-host model R_t estimates
lower95_two_host_sim<-qgamma(0.025,shape=estimate_two_host_sim$shape,scale=estimate_two_host_sim$scale)#lower bounds for the 95% intervals for the two-host model R_t estimates

#Generation of plots

p1 <- ggplot()+
  geom_line(data=df_true_value,aes(x=day,y=R_t,colour="True Value"))+
  geom_line(data=df_one_host_sim,aes(x=day,y=R_t,colour="Mean (One-host)"))+
  ggtitle("One-host model")+
  xlim(16,NA)+
  ylim(0,4)+
  geom_ribbon(data=df_one_host_sim,aes(ymin=lower95_one_host_sim,ymax=upper95_one_host_sim,x=day, y=R_t,fill="95% CrI"),alpha=0.3)+
  scale_colour_manual("",values=c("black","red"))+
  scale_fill_manual("",values="grey12") 

p2 <- ggplot()+
  geom_line(data=df_true_value,aes(x=day,y=R_t,colour="True Value"))+
  geom_line(data=df_two_host_sim,aes(x=day,y=R_t,colour="Mean"))+
  geom_line(data=df_two_host_sim_children,aes(x=day,y=R_t,colour="Mean (Children)"))+
  geom_line(data=df_two_host_sim_adults,aes(x=day,y=R_t,colour="Mean (Adults)"))+
  geom_ribbon(data=df_two_host_sim,aes(ymin=lower95_two_host_sim,ymax=upper95_two_host_sim,x=day, y=R_t,fill="95% CrI"),alpha=0.3)+
  ggtitle("Two-host model")+
  xlim(16,NA)+
  ylim(0,4)+
  scale_colour_manual("",values=c("black","green","blue","red"))+
  scale_fill_manual("",values="grey12")


multiplot(p1,p2,cols=1) #run for phi=1, k=1 (Figure 1 b) in the report, phi=0.75, 0.5, 0.25 and k=1 (Figure 2 in the report), k=1.25, 1.5, 1.75 and phi = 1 (Figure S1 in Supplementary File E)


p3 <- ggplot()+
  geom_line(data=df_true_value,aes(x=day,y=R_t,colour="True Value"))+
  geom_line(data=df_one_host_sim_corrected,aes(x=day,y=R_t,colour="Mean (One-host)"))+
  ggtitle("One-host model")+
  xlim(16,NA)+
  ylim(0,4)+
  geom_ribbon(data=df_one_host_sim_corrected,aes(ymin=lower95_one_host_sim_corrected,ymax=upper95_one_host_sim_corrected,x=day, y=R_t,fill="95% CrI"),alpha=0.3)+
  scale_colour_manual("",values=c("black","red"))+
  scale_fill_manual("",values="grey12")

multiplot(p3,p2,cols=1) #run for phi=0.75, 0.5, 0.25 and k=1 (Figure 3 in the report)




#Testing models for the real-world data ####

#Aggregating age groups to build a 2x2 contact matrix 

c11<- c()
c12<- c()
c22 <- c()
c21<-c()
for (i in 1:4){
  c11[i] <- P[i]*sum(M_all[i,1:4])/sum(P[1:4])
  c12[i] <- P[i]*sum(M_all[i,5:16])/sum(P[1:4])
}
for(i in 5:16){
  c21[i-4] <- P[i]*sum(M_all[i,1:4])/sum(P[5:16])
  c22[i-4] <- P[i]*sum(M_all[i,5:16])/sum(P[5:16])
}
C11 <- sum(c11)
C12 <- sum(c12)
C21 <- sum(c21)
C22 <- sum(c22)

C <- matrix(c(C11,C12,C21,C22),nrow=2,ncol=2) #standard contact matrix weighted to contain only the relevant age groups (< 20 and >=20)

p_da <- 0.7 # adult case fatality ratio for community cases 
p_dc <- 0.58 # child case fatality ratio for community cases 
p_f <- 0.571 #probability of unsafe burial given death in the community
C[1,2] <- C[1,2]*(1+5.67*p_f*p_dc) #scaling up the child-to-adult contacts
C[2,2] <- C[2,2]*(1+5.67*p_f*p_da) #scaling up the adult-to-adult contacts
C_eff <- C #effective contact matrix 
lambda<-spectralRadius(C_eff) # the dominant eigenvalue of the effective contact matrix 

#entries of the eigenvector corresponding to lambda (used further in the code)

v11 <- eigen(t(C_eff))$vectors[,1][1]
v12 <- eigen(t(C_eff))$vectors[,1][2]

#Load incidence data

real_world_data <- as.data.frame(read_xlsx("./equateur.xlsx")) #daily incidence data by age for the 2020 EVD outbreak in Equateur Province, DRC 

real_world_data_formatted <- real_world_data %>% group_by(DateOnset)%>%summarise(I1=length(which(Age<20)),I2=length(which(Age>=20))) #daily incidence data split into < 20 years and >= 20 years for the 2020 EVD outbreak in Equateur Province, DRC 


I1_real_world <- real_world_data_formatted$I1 #child cases (<20 years)
I2_real_world<- real_world_data_formatted$I2 #adult cases (>=20 years)
cases_real_world <- I1_real_world+I2_real_world #total cases

#plotting the daily incidence counts by age
Children <- I1_real_world
Adults <- I2_real_world
df_incidence <- rbind(Children,Adults)
barplot(df_incidence,legend.text=TRUE,xlab="Days",ylab="Total Cases",ylim=c(0,10)) # Figure 4 a) in the report (phi=1)

n_days_real <- length(cases_real_world) #duration of the outbreak

m_real <- seq(1,n_days_real,by=1) 
w_lit_real <- discr_si(m_real,mu,sd) 
w_lit_real <-w_lit_real/sum(w_lit_real)

estimate_one_host_real <- one_host_model(n_days=n_days_real,cases=cases_real_world,w=w_lit_real)
estimate_one_host_real_corrected <- one_host_model(n_days=n_days_real,cases=cases_real_world,w=pop_si_dist_corrected(phi=0.25,m=m_real)) 

estimate_two_host_real <- two_host_model(I1=I1_real_world,I2=I2_real_world,n_days=n_days_real,cases=cases_real_world,phi=1,m=m_real) #no difference between children and adults serial intervals
estimate_two_host_75_real <- two_host_model(I1=I1_real_world,I2=I2_real_world,n_days=n_days_real,cases=cases_real_world,phi=0.75,m=m_real) # mu_SI_child = 0.75mu_SI_adult
estimate_two_host_50_real <- two_host_model(I1=I1_real_world,I2=I2_real_world,n_days=n_days_real,cases=cases_real_world,phi=0.5,m=m_real) # mu_SI_child = 0.5mu_SI_adult
estimate_two_host_25_real <- two_host_model(I1=I1_real_world,I2=I2_real_world,n_days=n_days_real,cases=cases_real_world,phi=0.25,m=m_real) # mu_SI_child = 0.25mu_SI_adult

t<-seq(tau+2,length.out=length(estimate_two_host_real$R_mean),by=1)

upper95_one_host_real <- qgamma(0.975,shape=estimate_one_host_real$shape,scale=estimate_one_host_real$scale) #lower bounds for the 95% intervals for the one-host model R_t estimates
lower95_one_host_real<-qgamma(0.025,shape=estimate_one_host_real$shape,scale=estimate_one_host_real$scale) #lower bounds for the 95% intervals for the one-host model R_t estimates

upper95_one_host_real_corrected <- qgamma(0.975,shape=estimate_one_host_real_corrected$shape,scale=estimate_one_host_real_corrected$scale) # with corrected population serial interval distribution
lower95_one_host_real_corrected <-qgamma(0.025,shape=estimate_one_host_real_corrected$shape,scale=estimate_one_host_real_corrected$scale) 

upper95_two_host_real<- qgamma(0.975,shape=estimate_two_host_real$shape,scale=estimate_two_host_real$scale) #lower bounds for the 95% intervals for the two-host model R_t estimates
lower95_two_host_real<-qgamma(0.025,shape=estimate_two_host_real$shape,scale=estimate_two_host_real$scale) #upper bounds for the 95% intervals for the two-host model R_t estimates

df_one_host_real <- data.frame(day=t,R_t=estimate_one_host_real$R_mean)
df_one_host_real_corrected <- data.frame(day=t,R_t=estimate_one_host_real_corrected$R_mean)

df_two_host_real<-data.frame(day=t,R_t=estimate_two_host_real$R_mean) 
df_two_host_real_children <- data.frame(day=t,R_t=estimate_two_host_real$R_mean*(C_eff[1,1]+C_eff[1,2])/lambda) 
df_two_host_real_adults<-data.frame(day=t,R_t=estimate_two_host_real$R_mean*(C_eff[2,1]+C_eff[2,2])/lambda)

df_two_host_75_real <- data.frame(day=t,R_t=estimate_two_host_75_real$R_mean) 
df_two_host_50_real <-data.frame(day=t,R_t=estimate_two_host_50_real$R_mean) 
df_two_host_25_real <- data.frame(day=t,R_t=estimate_two_host_25_real$R_mean) 

# Generating plots

p1<-ggplot()+
  geom_line(data=df_one_host_real,aes(x=day,y=R_t,colour="Mean (One-host)"))+
  geom_line(data=df_two_host_real,aes(x=day,y=R_t,colour="Mean (Two-host)"))+
  ylim(0,4)+
  xlim(16,NA)+
  geom_ribbon(data=df_one_host_real,aes(ymin=lower95_one_host_real,ymax=upper95_one_host_real,x=day, y=R_t,fill="95% CrI"),alpha=0.3)+
  scale_colour_manual("",values=c("black","red"))+
  scale_fill_manual("",values="grey12")

p2 <- ggplot()+
  geom_line(data=df_two_host_real,aes(x=day,y=R_t,colour="Mean"))+
  geom_line(data=df_two_host_real_children,aes(x=day,y=R_t,colour="Mean (Children)"))+
  geom_line(data=df_two_host_real_adults,aes(x=day,y=R_t,colour="Mean (Adults)"))+
  geom_ribbon(data=df_two_host_real,aes(ymin=lower95_two_host_real,ymax=upper95_two_host_real,x=day, y=R_t,fill="95% CrI"),alpha=0.3)+
  ylim(0,4)+
  xlim(16,NA)+
  scale_colour_manual("",values=c("black","green","blue"))+
  scale_fill_manual("",values="grey12")


multiplot(p1,p2,cols=1) #Figure 4 b) in the report

p3 <- ggplot()+
  geom_line(data=df_one_host_real,aes(x=day,y=R_t,colour="Mean (One-host)"))+
  geom_line(data=df_two_host_75_real,aes(x=day,y=R_t,colour="Mean (Two-Host)   mu_child = 0.75mu_adult"))+
  geom_line(data=df_two_host_50_real,aes(x=day,y=R_t,colour="Mean (Two-host)  mu_child = 0.5mu_adult"))+
  geom_line(data=df_two_host_25_real,aes(x=day,y=R_t,colour="Mean (Two-host)  mu_child = 0.25mu_adult"))+
  ylim(0,4)+
  xlim(16,NA)+
  geom_ribbon(data=df_one_host_real,aes(ymin=lower95_one_host_real,ymax=upper95_one_host_real,x=day, y=R_t,fill="95% CrI"),alpha=0.3)+
  scale_colour_manual("",values=c("black","blue","red","green"))+
  scale_fill_manual("",values="grey12")


p4 <- ggplot()+
  geom_line(data=df_one_host_real_corrected,aes(x=day,y=R_t,colour="Mean (One-host)"))+
  geom_line(data=df_two_host_25_real,aes(x=day,y=R_t,colour="Mean (Two-host)  mu_child = 0.25mu_adult"))+
  ylim(0,4)+
  xlim(16,NA)+
  geom_ribbon(data=df_one_host_real,aes(ymin=lower95_one_host_real_corrected,ymax=upper95_one_host_real_corrected,x=day, y=R_t,fill="95% CrI"),alpha=0.3)+
  scale_colour_manual("",values=c("black","red"))+
  scale_fill_manual("",values="grey12")


multiplot(p3,p4,cols=1) #Figure 5 in the report
