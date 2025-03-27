#############################################################################
############### Sauer et al. Diet-SIR script            #####################
############### Code developed by Erin L. Sauer         #####################
#############################################################################

#### Diet x disease S-I-R model with reinfection ####
library(deSolve)
library(tidyverse)
library(ggpubr)

## parameters taken/modified from Fleming-Davies et al. 2018 ############
beta.protein <- 0.00275
beta.lipid <- beta.protein*2

nu <- 0.25 
psi <-  0.01

## parameters taken/modified from Cressler & Adelman 2024 ##########
delta <- 0.01 
lamda <- 0.025 

## parameters based on canary data ##############
a <- 1.7 #based on data from Sudnick et al. 2025

#gamma -> recovery rate; using median recovery times from Perrine et al. 2025
gamma.protein <- 1/21 
gamma.lipid <- 1/30 

#  mu -> inflammation induced mortality rate; using data from Perrine et al. 2025 
mu.protein <- -log(1-1/10)/35 
mu.lipid <- -log(1-3/11)/35 

################### S-I-R model function #########################
SIR <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{ 
    N <- S+I1+I2+R1+R2
    dS <- lamda*(S+R1+R2)-beta*S*(I1+I2)-delta*S + psi*R1
    dI1 <- beta*S*(I1+I2)-I1*(gamma+delta+mu+nu)
    dR1 <- I1*gamma - delta*R1 - (beta*R1*(I1+I2)) - psi*R1
    dI2 <- (beta*R1*(I1+I2))-I2*((gamma*a)+delta+(mu/2)+(nu/2))
    dR2 <- (I2*(gamma*a))-delta*R2
    return(list(c(dS, dI1, dR1, dI2, dR2)))
  })
}

############### model running ################################
times <- 0:100

params.protein <- c(beta=beta.protein, a=a, delta=delta, lamda=lamda, psi=psi, nu=nu,
                    gamma=gamma.protein, mu=mu.protein)
params.lipid <- c(beta=beta.lipid, a=a, delta=delta, lamda=lamda, psi=psi, nu=nu,
                  gamma=gamma.lipid, mu=mu.lipid)

S=100; I1=1; R1=0; I2=0; R2=0

initial_state <- c(S=S, I1=I1, R1=R1, I2=I2, R2=R2)
model.lipid <- ode(initial_state, times, SIR, params.lipid)
summary(model.lipid)

initial_state <- c(S=S, I1=I1, R1=R1, I2=I2, R2=)
model.protein <- ode(initial_state, times, SIR, params.protein)
summary(model.protein)

######### plotting  #####
#cumulative infection & total population size

model.protein <- as.data.frame(model.protein)
model.protein <- model.protein %>% mutate(diet="protein")

model.lipid <- as.data.frame(model.lipid)
model.lipid <- model.lipid %>% mutate(diet="lipid")

models <- rbind(model.lipid,model.protein)
models <- models %>% group_by(time) %>% mutate(I.total = I1 + I2)
models <- models %>% group_by(time) %>% 
  mutate(total.birds = S+R1+R2+I1+I2)
models <- models %>% group_by(time) %>% mutate(I.Prop = (I.total/total.birds))

#population size
PopSize <- ggplot(models, aes(x = time, y = total.birds)) + 
  geom_line(aes(color = diet)) +
  theme_bw()+
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.8,0.9),
        legend.text= element_text(size=20),
        legend.title= element_blank()
  )+
  ylab("Population size") + xlab("Days")
#proportion infected
PInf <- ggplot(models, aes(x = time, y = I.Prop)) + 
  geom_line(aes(color = diet)) +
  theme_bw()+
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.8,0.9),
        legend.text= element_text(size=20),
        legend.title= element_blank()
  )+
  ylab("Proportion infected") + xlab("Days")

ggarrange(PopSize, PInf, ncol=2, nrow=1, labels=c("A","B"),
          font.label = list(size=20,face="bold"))
