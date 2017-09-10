############################################################
# Title:			             Hypercycle Model1 for AgroEV
# Project Descriptor:	     AgroEv
# Author:			             bmarron
# Origin Date:		         29 Aug 2017
################################################################

#http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004853



 
library(deSolve)


#######################
#Hypercycle_Test3.2
#######################
#Define grid of time steps using seq(from, to, by="step")
time_steps3.2 <- seq(0, 100, by = 1)

#Define initial agent population (value in each compartment)
agent_start3.2 <- c(a=1, b=1, c=1)

#Define dummy params
parms3.2 <- c(ka1=.1, 
            ka2=.1,
            ka3=.1,
            kb1=.1,
            kb2=.1,
            kb3=.1,
            kc1=.1,
            kc2=.1,
            kc3=.1
)

#The function for the system of coupled nonlinear equations
hypercycle_test3.2 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    da.dt = ((1-(a+b+c)/100)*ka1*a + (1-(a+b+c)/100)*ka2*a*c - ka3*a)
    db.dt = ((1-(a+b+c)/100)*kb1*b + (1-(a+b+c)/100)*kb2*b*a - ka3*b)
    dc.dt = ((1-(a+b+c)/100)*kc1*c + (1-(a+b+c)/100)*kc2*c*b - kc3*c)
    results = c(da.dt, db.dt, dc.dt)
    list(results)
  })
}


# lsoda{deSolve}
set.seed(47)
sink("/home/bruce/Desktop/Hypercycle_Test3.2.txt")
(as.data.frame(lsoda(agent_start3.2, time_steps3.2, hypercycle_test3.2, parms3.2)))
sink()


#######################
#Hypercycle_Test3.3
#######################
#Define grid of time steps using seq(from, to, by="step")
time_steps3.3 <- seq(0, 100, by = 1)

#Define initial agent population (value in each compartment)
agent_start3.3 <- c(a=1, b=1, c=1)

#Define dummy params
parms3.3 <- c(ka1=.1, 
              ka2=.5,
              ka3=.05,
              kb1=.1,
              kb2=.5,
              kb3=.05,
              kc1=.1,
              kc2=.5,
              kc3=.05
)

#The function for the system of coupled nonlinear equations
hypercycle_test3.3 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    da.dt = ((1-(a+b+c)/1000)*ka1*a + (1-(a+b+c)/1000)*ka2*a*c - ka3*a)
    db.dt = ((1-(a+b+c)/1000)*kb1*b + (1-(a+b+c)/1000)*kb2*b*a - ka3*b)
    dc.dt = ((1-(a+b+c)/1000)*kc1*c + (1-(a+b+c)/1000)*kc2*c*b - kc3*c)
    results = c(da.dt, db.dt, dc.dt)
    list(results)
  })
}


# lsoda{deSolve}
set.seed(47)
sink("/home/bruce/Desktop/Hypercycle_Test3.3.txt")
(as.data.frame(lsoda(agent_start3.3, time_steps3.3, hypercycle_test3.3, parms3.3)))
sink()



#######################
#Hypercycle_Test3.4
#######################
#Define grid of time steps using seq(from, to, by="step")
time_steps3.4 <- seq(0, 100, by = 1)

#Define initial agent population (value in each compartment)
agent_start3.4 <- c(a=1, b=1, c=1)

#Define dummy params
set.seed(47)
parms3.4 <- c(ka1=runif(1, -1, 1), 
              ka2=runif(1, -1, 1),
              ka3=runif(1, -1, 1),
              kb1=runif(1, -1, 1),
              kb2=runif(1, -1, 1),
              kb3=runif(1, -1, 1),
              kc1=runif(1, -1, 1),
              kc2=runif(1, -1, 1),
              kc3=runif(1, -1, 1)
)

#The function for the system of coupled nonlinear equations
hypercycle_test3.4 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    da.dt = ((1-(a+b+c)/1000)*ka1*a + (1-(a+b+c)/1000)*ka2*a*c - ka3*a)
    db.dt = ((1-(a+b+c)/1000)*kb1*b + (1-(a+b+c)/1000)*kb2*b*a - ka3*b)
    dc.dt = ((1-(a+b+c)/1000)*kc1*c + (1-(a+b+c)/1000)*kc2*c*b - kc3*c)
    results = c(da.dt, db.dt, dc.dt)
    list(results)
  })
}


# lsoda{deSolve}
sink("/home/bruce/Desktop/Hypercycle_Test3.4.txt")
(as.data.frame(lsoda(agent_start3.4, time_steps3.4, hypercycle_test3.4, parms3.4)))
sink()


#######################
#Hypercycle_Test3.5
#######################
#Define grid of time steps using seq(from, to, by="step")
time_steps3.5 <- seq(0, 100, by = 1)

#Define initial agent population (value in each compartment)
agent_start3.5 <- c(a=1, b=1, c=1)

#Define dummy params
set.seed(47)
parms3.5 <- c(ka1=runif(1, 0, 1), 
              ka2=runif(1, 0, 1),
              ka3=runif(1, 0, 1),
              kb1=runif(1, 0, 1),
              kb2=runif(1, 0, 1),
              kb3=runif(1, 0, 1),
              kc1=runif(1, 0, 1),
              kc2=runif(1, 0, 1),
              kc3=runif(1, 0, 1)
)

#The function for the system of coupled nonlinear equations
hypercycle_test3.5 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    da.dt = ((1-(a+b+c)/1000)*ka1*a + (1-(a+b+c)/1000)*ka2*a*c - ka3*a)
    db.dt = ((1-(a+b+c)/1000)*kb1*b + (1-(a+b+c)/1000)*kb2*b*a - ka3*b)
    dc.dt = ((1-(a+b+c)/1000)*kc1*c + (1-(a+b+c)/1000)*kc2*c*b - kc3*c)
    results = c(da.dt, db.dt, dc.dt)
    list(results)
  })
}


# lsoda{deSolve}
sink("/home/bruce/Desktop/Hypercycle_Test3.5.txt")
(as.data.frame(lsoda(agent_start3.5, time_steps3.5, hypercycle_test3.5, parms3.5)))
sink()


############################
# Processing Hypercycle 3.5
# (Linux)
#############################
'''
   a. remove all whitespaces in tmp.txt; pipe out
   b. change all whitespaces to ',' and send output to tmp2.txt
   c. delete first line of tmp1.txt

$ tr -s ' ' < Hypercycle_Test3.5.txt | tr ' ' ',' > tmp1.txt


   d. remove Field 1 in temp2.txt (a comma-separated file) and keep the rest
   e. send output to tmp.csv

$ cut -d',' -f1 --complement temp1.txt >> Hypercycle_3.5.csv

'''

#   f. import Hypercycle_3.5.csv into R
library(readr)
Hypercycle_3_5 <- read_csv("~/Desktop/Hypercycle_3.5.csv", col_names = FALSE)
View(Hypercycle_3_5)



