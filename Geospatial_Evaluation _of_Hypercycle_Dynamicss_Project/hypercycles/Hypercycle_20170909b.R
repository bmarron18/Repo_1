############################################################
# Title:			             Hypercycle Model1 for AgroEV
# Project Descriptor:	     AgroEv
# Author:			             bmarron
# Origin Date:		         09 Sept 2017
################################################################


##############
# Notes
#############
'''
 http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004853

 https://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution
 Normal dist has maximum entropy among all real-valued distributions supported on (-Inf, +inf)
 with a specified variance σ2 (a particular moment). Therefore, the assumption of normality 
 imposes the minimal prior structural constraint beyond this moment


 https://en.wikipedia.org/wiki/Weibull_distribution
 The Weibull distribution is the maximum entropy distribution for a non-negative real random variate 
 with a fixed expected value of xk equal to λk and a fixed expected value of ln(xk ) equal to 
 ln(λk ) - \gamma.

'''

##################
# Weibull Dist
##################

'''
https://en.wikipedia.org/wiki/Weibull_distribution
 The Weibull distribution is the maximum entropy distribution for a non-negative real random variate 
 with a fixed expected value of xk equal to λk and a fixed expected value of ln(xk ) equal to 
 ln(λk ) - \gamma.

k > 0 is the shape parameter and λ > 0 is the scale parameter 
'''


dweibull(x, shape, scale = 1, log = FALSE)
pweibull(q, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
qweibull(p, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
rweibull(n, shape, scale = 1)


library(Cairo)
CairoFonts(regular="DejaVu Sans:style=Regular")
CairoSVG("Weibull PDF.svg")
par(mar=c(3, 3, 1, 1))

x <- seq(0, 10, length.out=1000)
plot(x, dweibull(x, .5), type="l", col="blue", xlab="", ylab="", xlim=c(0, 4), ylim=c(0, 3), xaxs="i", yaxs="i")
lines(x, dweibull(x, 7), type="l", col="red")
lines(x, dweibull(x, 1.25), type="l", col="magenta")
lines(x, dweibull(x, 5), type="l", col="green")
legend("topright", legend=paste("\u03bb = 1, k =", c(.5, 1, 1.5, 5)), lwd=1, col=c("blue", "red", "magenta", "green"))





#######################
#Hypercycle_Test5.0
#######################

library(deSolve)

#Define grid of time steps using seq(from, to, by="step")
time_steps5.0 <- seq(0, 500, by = 1)

#Define initial agent population (value in each compartment)
agent_start5.0 <- c(a=1, b=1, c=1)

#Define dummy params
#Explore Weibull parameters ==> rweibull(n, shape, scale = 1)


#The function for the system of coupled nonlinear equations
hypercycle_test5.0 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    da.dt = ((1-(a+b+c)/1000)*ka1*a + (1-(a+b+c)/1000)*ka2*a*c - ka3*a)
    db.dt = ((1-(a+b+c)/1000)*kb1*b + (1-(a+b+c)/1000)*kb2*b*a - ka3*b)
    dc.dt = ((1-(a+b+c)/1000)*kc1*c + (1-(a+b+c)/1000)*kc2*c*b - kc3*c)
    results = c(da.dt, db.dt, dc.dt)
    list(results)
  })
}

# create the sets of params
# A set of ten (10) random pulls for each param
test5.0 <- list(ka1=rweibull(100, 1.25, scale = 1), 
              ka2=rweibull(100, 1.25, scale = 1),
              ka3=rweibull(100, 1.25, scale = 1),
              kb1=rweibull(100, 1.25, scale = 1),
              kb2=rweibull(100, 1.25, scale = 1),
              kb3=rweibull(100, 1.25, scale = 1),
              kc1=rweibull(100, 1.25, scale = 1),
              kc2=rweibull(100, 1.25, scale = 1),
              kc3=rweibull(100, 1.25, scale = 1)
)

# transfer the sets of parameters to a nested list
compiled_test5.0 <-c()
for (i in 1:100){
  compiled_test5.0[i]<-list(c(
                        ka1=test5.0[[1]][i], 
                        ka2=test5.0[[2]][i],
                        ka3=test5.0[[3]][i],
                        kb1=test5.0[[4]][i],
                        kb2=test5.0[[5]][i],
                        kb3=test5.0[[6]][i],
                        kc1=test5.0[[7]][i],
                        kc2=test5.0[[8]][i],
                        kc3=test5.0[[9]][i])
                       )
                 }


# https://stats.idre.ucla.edu/r/codefragments/multiple_filenames/
# do 10 runs of the hypercycle (one run for each param set); output of each run to its own file
for (i in 1:100) {
  a <- as.data.frame(lsoda(agent_start5.0, time_steps5.0, hypercycle_test5.0, compiled_test5.0[[i]]))
  mytime <- format(Sys.time(), "%b_%d_%H_%M_%S_%Y")
  myfile <- file.path(getwd(), paste0(mytime, "_", i, ".txt"))
  write.table(a, file = myfile, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = FALSE, append = FALSE)
}


