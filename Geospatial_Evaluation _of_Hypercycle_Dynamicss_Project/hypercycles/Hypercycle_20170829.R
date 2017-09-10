############################################################
# Title:			             Hypercycle Model1 for AgroEV
# Project Descriptor:	     AgroEv
# Author:			             bmarron
# Origin Date:		         29 Aug 2017
################################################################

#http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004853



 
library(deSolve)

###################
#Hypercycle_Test1
####################
#Define grid of time steps using seq(from, to, by="step")
time_steps1 <- seq(0, 200, by = 1)

#Define initial agent population (value in each compartment)
agent_start1 <- c(a1=1, a2=1, a3=1)

#Define dummy params
parms1 <- c(k=1)

#The function for the system of coupled nonlinear equations
hypercycle_test1 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
 da1.dt = (runif(1, -1, 1)*a1*k + runif(1, -1, 1)*a1*a3)
 da2.dt = (runif(1, -1, 1)*a2 + runif(1, -1, 1)*a2*a1)
 da3.dt = (runif(1, -1, 1)*a3 + runif(1, -1, 1)*a3*a2)
results = c(da1.dt, da2.dt, da3.dt)
list(results)
  })
}


# lsoda{deSolve}
set.seed(47)
sink("/home/bruce/Desktop/Hypercycle_Test1.txt")
(as.data.frame(lsoda(agent_start1, time_steps1, hypercycle_test1, parms1)))
sink()


#######################
#Hypercycle_Test2
#######################
#Define grid of time steps using seq(from, to, by="step")
time_steps2 <- seq(0, 200, by = 1)

#Define initial agent population (value in each compartment)
agent_start2 <- c(a1=1, a2=1, a3=1)

#Define dummy params
parms2 <- c(k1=runif(1, -1, 1), 
           k2=runif(1, -1, 1),
           k3=runif(1, -1, 1),
           k4=runif(1, -1, 1),
           k5=runif(1, -1, 1),
           k6=runif(1, -1, 1)
           )

#The function for the system of coupled nonlinear equations
hypercycle_test2 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    da1.dt = (k1*a1 + k2*a1*a3)
    da2.dt = (k3*a2 + k4*a2*a1)
    da3.dt = (k5*a3 + k6*a3*a2)
    results = c(da1.dt, da2.dt, da3.dt)
    list(results)
  })
}


# lsoda{deSolve}
set.seed(47)
sink("/home/bruce/Desktop/Hypercycle_Test2.txt")
(as.data.frame(lsoda(agent_start2, time_steps2, hypercycle_test2, parms2)))
sink()


#######################
#Hypercycle_Test3
# no change!
#######################
#Define grid of time steps using seq(from, to, by="step")
time_steps3 <- seq(0, 200, by = 1)

#Define initial agent population (value in each compartment)
agent_start3 <- c(a1=1, a2=1, a3=1)

#Define dummy params
parms3 <- c(k1=-.1, 
            k2=.1,
            k3=-.1,
            k4=.1,
            k5=-.1,
            k6=.1
)

#The function for the system of coupled nonlinear equations
hypercycle_test3 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    da1.dt = (k1*a1 + k2*a1*a3)
    da2.dt = (k3*a2 + k4*a2*a1)
    da3.dt = (k5*a3 + k6*a3*a2)
    results = c(da1.dt, da2.dt, da3.dt)
    list(results)
  })
}


# lsoda{deSolve}
set.seed(47)
sink("/home/bruce/Desktop/Hypercycle_Test3.txt")
(as.data.frame(lsoda(agent_start3, time_steps3, hypercycle_test3, parms3)))
sink()


#######################
#Hypercycle_Test4
#######################
#Define grid of time steps using seq(from, to, by="step")
time_steps4 <- seq(0, 10, by = 1)

#Define initial agent population (value in each compartment)
agent_start4 <- c(a1=1, a2=1, a3=1)

#Define dummy params
parms4 <- c(k1=runif(1, -1, 1), 
            k2=runif(1, -1, 1),
            k3=runif(1, -1, 1),
            k4=runif(1, -1, 1),
            k5=runif(1, -1, 1),
            k6=runif(1, -1, 1)
)

#The function for the system of coupled nonlinear equations
hypercycle_test4 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    da1.dt = (k1*a1 + k2*a1*a3)
    da2.dt = (k3*a2 + k4*a2*a1)
    da3.dt = (k5*a3 + k6*a3*a2)
    results = c(da1.dt, da2.dt, da3.dt)
    list(results)
  })
}


# lsoda{deSolve}
set.seed(47)
sink("/home/bruce/Desktop/Hypercycle_Test3.txt")
(as.data.frame(lsoda(agent_start4, time_steps4, hypercycle_test4, parms4)))
sink()


#######################
#Hypercycle_Test5
#######################
set.seed(47)

agent_start4 <- c(a1=1, a2=1, a3=1)

#The function for the system of coupled nonlinear equations
hypercycle_test4 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    da1.dt = (k1*a1 + k2*a1*a3)
    da2.dt = (k3*a2 + k4*a2*a1)
    da3.dt = (k5*a3 + k6*a3*a2)
    results = c(da1.dt, da2.dt, da3.dt)
    list(results)
  })
}

time_steps4 <- c()
time_steps4[1] <- 1
for (i in 2:2){
  time_steps4[i] <- i
  parms4 <- c(k1=runif(1, -1, 1), 
              k2=runif(1, -1, 1),
              k3=runif(1, -1, 1),
              k4=runif(1, -1, 1),
              k5=runif(1, -1, 1),
              k6=runif(1, -1, 1)
  )
  
  (as.data.frame(lsoda(agent_start4, time_steps4[i], hypercycle_test4, parms4)))
}



# lsoda{deSolve}
sink("/home/bruce/Desktop/Hypercycle_Test4.txt")
(as.data.frame(lsoda(agent_start4, time_steps4, hypercycle_test4, parms4)))
sink()

###########
time_steps4 <- c()
time_steps4[1] <- 1

for (i in 2:17){
  time_steps4[i] <- i
}
