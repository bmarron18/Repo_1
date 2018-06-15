############################################################
# Title:			             Hypercycle_5.0-to-rasters
# Project Descriptor:	     AgroEv
# Author:			             bmarron
# Origin Date:		         17 Sept 2017
################################################################


##############
# websites
###############

'''
https://stats.idre.ucla.edu/r/codefragments/multiple_filenames/
https://stackoverflow.com/questions/4309217/cartesian-product-data-frame-in-r
https://stat.ethz.ch/pipermail/r-sig-geo/2012-September/016175.html
http://rstudio-pubs-static.s3.amazonaws.com/7993_6b081819ba184047802a508a7f3187cb.html
https://gis.stackexchange.com/questions/20018/how-can-i-convert-data-in-the-form-of-lat-lon-value-into-a-raster-file-using-r
http://neondataskills.org/R/Multi-Band-Rasters-In-R/
http://neondataskills.org/R/Raster-Times-Series-Data-In-R/
https://gis.stackexchange.com/questions/126811/combine-rasters-into-single-raster-with-multiple-field-values-in-r
https://artax.karlin.mff.cuni.cz/r-help/library/raster/html/writeRaster.html
https://gis.stackexchange.com/questions/207238/how-to-mix-two-raster-in-r-same-extent
http://idlcoyote.com/ip_tips/where3.html

'''


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%% Hypercycles-to-Time Series Raster Maps %%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'''
The following step-wise algorithm generates time series output from a degree p=2, n=3
hypercycle (p=degree of the growth functions; n=number of entities) and encapsulates 
this output in multi-band (bands=3) raster maps.

Eigen, Manfred, and Peter Schuster. 1979. 
The Hypercycle, a Principle of Natural Self-Organization. 
Berlin; New York: Springer-Verlag.
'''

# The Weibull dist is often used in survival analysis
# k = shape parameter, λ = scale parameter (often set to 1 for ease of scale)

# k < 1 ==> failure/mortality/linkage decreses thru time
# k = 1 ==> failure/mortality/linkage constant thru time (Weibull reduces to exponential)
# k > 1 ==> failure/mortality/linkage increases thru time 
# k ==> inf; Weibull reduces to Dirac delta

x <- seq(0, 10, length.out=1000)
plot(x, dweibull(x, .5, 1), type="l", col="blue", xlab="", ylab="", xlim=c(0, 4), ylim=c(0, 3), xaxs="i", yaxs="i", main="")
title(main="Weibull Distributions")
lines(x, dweibull(x, 1.25, 1), type="l", col="magenta")
lines(x, dweibull(x, 5, 1), type="l", col="green")
lines(x, dweibull(x, 7, 1), type="l", col="red")
lines(x, dweibull(x, 2000, 1), type="l", col="black")
legend("topright", legend=paste("\u03bb = 1, k =", c(.5, 1.25, 5, 7, 2000)), lwd=1, col=c("blue", "magenta", "green", "red", "black" ))



##################
#Libraries
##################
library(deSolve)
library(raster)
library(readr)
library(data.table)
library(rgdal)




#################################################
#Generate Monte Carlo output from Hypercycle
##################################################

library(deSolve)

#STEP1
#Define a set of 500 time steps in 1 tick increments
time_steps <- seq(0, 500, by = 1)

#STEP2
#Define the initial entity populations
entity_start <- c(a=1, b=1, c=1, d=1, e=1)

#STEP3
#Define the hypercycle: A function of coupled, nonlinear equations
#Hypercycle is a p=2 (# catalytic outputs), n=5 (# entities) hypercycle
#entities a, b, c, d, e
hypercycle <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    da.dt = ((1-(a+b+c+d+e)/1000)*ka1*a + (1-(a+b+c+d+e)/1000)*ka2*a*e - ka3*a)
    db.dt = ((1-(a+b+c+d+e)/1000)*kb1*b + (1-(a+b+c+d+e)/1000)*kb2*b*a - kb3*b)
    dc.dt = ((1-(a+b+c+d+e)/1000)*kc1*c + (1-(a+b+c+d+e)/1000)*kc2*c*b - kc3*c)
    dd.dt = ((1-(a+b+c+d+e)/1000)*kd1*d + (1-(a+b+c+d+e)/1000)*kc2*d*c - kd3*d)
    de.dt = ((1-(a+b+c+d+e)/1000)*ke1*e + (1-(a+b+c+d+e)/1000)*kc2*e*d - ke3*e)
    results = c(da.dt, db.dt, dc.dt, dd.dt, de.dt)
    list(results)
  })
}

#STEP4a
# Create 100 sets of unique hypercycle params. (initial conditions)
# Each set contains 15 parameter values obtained as random draws (Monte Carlo sampling)
# from the the Weibull dist (k=shape=1.25, λ=scale=1)

set.seed(47)
test_params <- list(
              ka1=rweibull(100, 1.25, scale = 1), 
              ka2=rweibull(100, 1.25, scale = 1),
              ka3=rweibull(100, 1.25, scale = 1),
              kb1=rweibull(100, 1.25, scale = 1),
              kb2=rweibull(100, 1.25, scale = 1),
              kb3=rweibull(100, 1.25, scale = 1),
              kc1=rweibull(100, 1.25, scale = 1),
              kc2=rweibull(100, 1.25, scale = 1),
              kc3=rweibull(100, 1.25, scale = 1),
              kd1=rweibull(100, 1.25, scale = 1),
              kd2=rweibull(100, 1.25, scale = 1),
              kd3=rweibull(100, 1.25, scale = 1),          
              ke1=rweibull(100, 1.25, scale = 1),
              ke2=rweibull(100, 1.25, scale = 1),
              ke3=rweibull(100, 1.25, scale = 1)             
              )

str(test_params)

#STEP4b
# transfer the 100 sets of parameters produced in STEP4a to a nested list
# so that each set (of 15 params each) can be accessed individually
compiled_test_params <-c()
for (i in 1:100){
  compiled_test_params[i]<-list(c(
                        ka1=test_params[[1]][i], 
                        ka2=test_params[[2]][i],
                        ka3=test_params[[3]][i],
                        kb1=test_params[[4]][i],
                        kb2=test_params[[5]][i],
                        kb3=test_params[[6]][i],
                        kc1=test_params[[7]][i],
                        kc2=test_params[[8]][i],
                        kc3=test_params[[9]][i],
                        kd1=test_params[[10]][i], 
                        kd2=test_params[[11]][i],
                        kd3=test_params[[12]][i],
                        ke1=test_params[[13]][i],
                        ke2=test_params[[14]][i],
                        ke3=test_params[[15]][i])
                       )
                 }


head(compiled_test_params)


#STEP5a
#make a directory to hold raw output from hypercycle runs (Linux)
mkdir RawData

#STEP5b
# Do 100 simulation RUNS of the hypercycle model (one run for each unique param set;
# ie one run for each compiled_test_params[[i]]) for 500 TICKS per RUN
# Send the output from each RUN to its own file
#
# Each simulation RUN will be a raster cell so this step will create simulation data 
# to populate 100 cells (a 10x10 grid) with Z=5 values (bands) per cell
# Each simulation TICK is a timestep in the time series of raster maps

for (i in 1:100) {
  h <- as.data.frame(lsoda(entity_start, time_steps, hypercycle, compiled_test_params[[i]]))
  mytime <- format(Sys.time(), "%b_%d_%H_%M_%S_%Y")
  myfile <- file.path(paste(getwd(),"RawData", sep="/"), paste0(mytime, "_", i, ".txt"))
  write.table(h, file = myfile, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = FALSE, append = FALSE)
}



#STEP6a
# Post-processing of hypercycle output data in directory RawData (Linux)
# sed ==> remove first row (line 1)
# cut ==> 1) remove column 1 (field 1) and keep the rest, 2) ouput to *.txt.flag files
# && ... remove *.txt files (leaves only the *.txt.flag files)

for i in *.txt; do 
sed -i '1d' $i;
cut -d "," --complement -f1 $i > $i.flag;
done && rm *.txt

#STEP6b 
# Post-processing of hypercycle output data (Linux)
# change *.txt.flag files to *.csv files
find -type f -name '*.txt.flag' | while read FILE ; do 
newfile="$(echo ${FILE} | sed -e 's/txt.flag/csv/')" ;
mv "${FILE}" "${newfile}" ;
done




############################################
# Create time-series raster data w/ bands
############################################

#STEP7
#load all * .csv files as a list
hypercycle_data_path <- "~/Desktop/NCSU/2017_DFER_06_GeospatialDataMining/Works_InProgress/Geospatial_Eval _of_Hypercycle_Dynamics/hypercycles/RawData" 
all_hypercycle <- list.files(hypercycle_data_path,
                             full.names = TRUE,
                             pattern = ".csv")

summary(all_hypercycle)

# view and check list
# note the full path, relative to our working directory, is included
all_hypercycle

#STEP8
# make 'n' directories to hold raster outputs (Linux);
# one directory per timestep (here 5)
# Make sure to set variable 't_steps' to same 'n'
for num in {1..200}; do
mkdir timebricks_t$num
done

#STEP9
# options for raster building:
# i counts timesteps for the same raster cell OR
# i counts raster cells for the same timestep   <== choose this

library(raster)
library(readr)
library(data.table)
library(rgdal)

# Every simulation RUN is now a *.csv file with 501 lines (Line 1 = time0)
# Every *.csv will hold the time series data for ONE raster cell
# Every row of data in a *.csv is a timestep (a simulation TICK)
# 
# j counts timesteps
# d is a Cartesian product that gives cell coordinates for a 10x10 grid
# i counts raster cells (for the same timestep); there are 100 raster cells
# p,q,r,s,t are five, individual rasters (each 10x10) that will each hold ONE
# entity attribute (Z1, Z2, Z3, Z4, Z5). The attributes will be transferred to raster cell 
# attributes (Z1, Z2, Z3, Z4, Z5) and then stacked/packed into a single, multi-band raster
# 
# u is a multi-band raster object that stacks/packs attributes Z1, Z2, Z3, Z4, Z5 from 
# individual rasters p, q, r, s, t and makes a single, 5-band raster object
#
# Each individual u object is a 5-band raster (*.tif) that contains the Z1,Z2,Z3,Z4,Z5 data 
# for a single cell at a single timestep. All of the 5-band rasters for a given timestep 
# are collected into 'timebrick' directories; one timebrick folder per timestep
# 
#set the number of timesteps to be processed into maps; one map per timestep
# variable 't_steps' to be used in STEP11 also

t_steps <- 200
for (j in 1:t_steps){
  vector1<- c(1:10)
  d <- CJ(vector1, vector1)
  for (i in 1:100){
    a<-read_csv(all_hypercycle[i], col_names = FALSE)
    p<-raster(nrow=10,ncol=10)
    p[d$V1[i], d$V2[i]]<-a[[1]][j]
    q<-raster(nrow=10,ncol=10)
    q[d$V1[i], d$V2[i]]<-a[[2]][j]
    r<-raster(nrow=10,ncol=10)
    r[d$V1[i], d$V2[i]]<-a[[3]][j]
    s<-raster(nrow=10,ncol=10)
    s[d$V1[i], d$V2[i]]<-a[[4]][j]
    t<-raster(nrow=10,ncol=10)
    t[d$V1[i], d$V2[i]]<-a[[5]][j]
    u = brick(p,q,r,s,t)
    
    mytime <- format(Sys.time(), "%b_%d_%Y_%H_%M")
    mypath <- file.path(paste(getwd(), paste0("timebricks", "_t", j), sep="/"))
    myfile <- paste0(mytime, "_cell", i, ".tif")
    b<-getwd()
    setwd(mypath)
    writeRaster(u, filename= myfile, format="GTiff", overwrite=TRUE)
#   writeRaster(u, filename= myfile, bandorder='BIL', format="raster", overwrite=TRUE)
    setwd(b)
  }
}


##################################################
# Create a single map per timestep 
##################################################


# STEP10
# make 'n' directories to hold raster time series outputs (Linux);
# one directory per time series map (here 5)
# No. directories should be same value as variable 't_steps'
for num in {1..200}; do
mkdir TimeSeries_t$num
done


#STEP11
# Stack the individual cell maps to create one, 5-band raster of 100 cells
# per time step

for (j in 1:t_steps){
 b<-getwd()
 mypath1 <- file.path(paste(b, paste0("timebricks", "_t", j), sep="/"))
 tif_maps <- list.files(mypath1, full.names=TRUE)

 x <- vector("list", 100)
 for(i in 1:100) {
   m <- stack(tif_maps[[i]])
   x[[i]] <- m
  }
 time_map <- do.call(merge, x)

 mytime <- format(Sys.time(), "%b_%d_%Y")
 #mypath2 <- file.path(paste(b, paste0("TimeSeries", "_t", j), sep="/"))
 mypath2 <- file.path(paste(b, paste0("TimeSeries"), sep="/"))
 myfile <- paste0(mytime, "_TimeSeriesMap", j, ".tif")
 setwd(mypath2)
 writeRaster(time_map, filename= myfile, format="GTiff", overwrite=TRUE)
 setwd(b)
}


##############
# Misc checks
##############
b<-getwd()

  # check timebrick outputs
setwd("/home/bruce/Desktop/NCSU/2017_DFER_06_GeospatialDataMining/Works_InProgress/Geospatial_Eval _of_Hypercycle_Dynamics/hypercycles/timebricks_t1/")
test<-raster("Sep_14_2017_13_56_cell1.tif")
test2<-raster("Sep_14_2017_13_56_cell2.tif")
test3<-raster("Sep_14_2017_13_56_cell100.tif")

test
test2
test3

plot(test)
plot(test2)
plot(test3)


  # check TimeSeries outputs
setwd("/home/bruce/Desktop/NCSU/2017_DFER_06_GeospatialDataMining/Works_InProgress/Geospatial_Eval _of_Hypercycle_Dynamics/hypercycles/TimeSeries_t1/")
test4<- raster("Sep_17_2017_16_36_TimeSeriesMap1.tif")
test4
plot(test4)

setwd("/home/bruce/Desktop/NCSU/2017_DFER_06_GeospatialDataMining/Works_InProgress/Geospatial_Eval _of_Hypercycle_Dynamics/hypercycles/TimeSeries_t5/")
test5<- raster("Sep_17_2017_16_37_TimeSeriesMap5.tif")
test5
plot(test5)

plot(test5, colNA=1)

  # reset wd
setwd(b)
-