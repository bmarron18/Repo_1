############################################################
# Title:			             Hypercycle-to-maps
# Project Descriptor:	     AgroEv
# Author:			             bmarron
# Origin Date:		         14 Sept 2017
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
#%%%%%%%%%%%% Hypercycles-to-Time Series Maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#######################
#Hypercycle_5.0
#######################

library(deSolve)

#STEP1
#Define a set of time steps using seq(from, to, by="step")
time_steps5.0 <- seq(0, 500, by = 1)

#STEP2
#Define initial agent population (value in each compartment)
agent_start5.0 <- c(a=1, b=1, c=1)

#STEP3
#Define the hypercycle: The function for the system of coupled nonlinear equations
#A ... hypercycle
hypercycle_test5.0 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    da.dt = ((1-(a+b+c)/1000)*ka1*a + (1-(a+b+c)/1000)*ka2*a*c - ka3*a)
    db.dt = ((1-(a+b+c)/1000)*kb1*b + (1-(a+b+c)/1000)*kb2*b*a - ka3*b)
    dc.dt = ((1-(a+b+c)/1000)*kc1*c + (1-(a+b+c)/1000)*kc2*c*b - kc3*c)
    results = c(da.dt, db.dt, dc.dt)
    list(results)
  })
}

#STEP4a
# Create the sets of params for the hypercycle from randomly sampling the Weibull dist
# This produces a set of one-hundred (100) random pulls for each param
set.seed(47)
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


#STEP4b
# transfer the sets of parameters produced in STEP4a to a nested list
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


#STEP5
# do 100 runs of the hypercycle model (one run for each unique param set) and send the 
# output of each run to its own file
# each simulation run is a raster cell
# each simulation tick is a timestep
for (i in 1:100) {
  h <- as.data.frame(lsoda(agent_start5.0, time_steps5.0, hypercycle_test5.0, compiled_test5.0[[i]]))
  mytime <- format(Sys.time(), "%b_%d_%H_%M_%S_%Y")
  myfile <- file.path(paste(getwd(),"data", sep="/"), paste0(mytime, "_", i, ".txt"))
  write.table(h, file = myfile, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = FALSE, append = FALSE)
}



#STEP6a
# Post-processing of hypercycle output data (Linux)
# do sed ...remove first row (line 1)
# cut -d ... 1) remove column 1 (field 1) and keep the rest, 2) ouput to .txt.flag files
# && ... remove .txt files (leaves only the .txt.flag files)

for i in *.txt;
do sed -i '1d' $i;
cut -d "," --complement -f1 $i > $i.flag;
done && rm *.txt

#STEP6b 
# Data post-processing of hypercycle output data (Linux)
# change .txt.flag files to .csv files
find -type f -name '*.txt.flag' | while read FILE ;
do newfile="$(echo ${FILE} |sed -e 's/txt.flag/csv/')" ;
mv "${FILE}" "${newfile}" ;
done




############################################
# Create time-series raster data w/ bands
############################################

#STEP7
#load all .csv files as a list
hypercycle_data_path <- "~/Desktop/NCSU/2017_DFER_06_GeospatialDataMining/Works_InProgress/Geospatial_Eval _of_Hypercycle_Dynamics/hypercycles/data" 
all_hypercycle <- list.files(hypercycle_data_path,
                             full.names = TRUE,
                             pattern = ".csv")

# view list - note the full path, relative to our working directory, is included
all_hypercycle

#STEP8
#make directories to hold raster outputs (Linux)
for num in {1..5}; do
mkdir timebricks_t$num
done

#STEP9
# options:
# i counts timesteps for the same raster cell OR
# i counts raster cells for the same timestep (choose this)

library(readr)
library(data.table)
library(rgdal)

# every sim run is a .csv
# every .csv is a cell
# every row of data in the .csv is a timestep
# 
# j counts the number of sim runs (=timesteps)
# d is a Cartesian product that gives cell coordinates
# i counts raster cells for the same timestep
# r,s,t are attributes (Z1, Z2, Z3)
# u is a multi-layer raster object

for (j in 1:5){
  vector1<- c(1:10)
  d <- CJ(vector1, vector1)
  for (i in 1:100){
    a<-read_csv(all_hypercycle[i], col_names = FALSE)
    r<-raster(nrow=10,ncol=10)
    r[d$V1[i], d$V2[i]]<-a[[1]][j]
    s<-raster(nrow=10,ncol=10)
    s[d$V1[i], d$V2[i]]<-a[[2]][j]
    t<-raster(nrow=10,ncol=10)
    t[d$V1[i], d$V2[i]]<-a[[3]][j]
    u = brick(r,s,t)
    
    mytime <- format(Sys.time(), "%b_%d_%Y_%H_%M")
    mypath <- file.path(paste(getwd(), paste0("timebricks", "_t", j), sep="/"))
    myfile <- paste0(mytime, "_cell", i, ".tif")
    b<-getwd()
    setwd(mypath)
    writeRaster(u, filename= myfile, format="GTiff", overwrite=TRUE)
    #    writeRaster(u, filename= myfile, bandorder='BIL', format="raster", overwrite=TRUE)
    setwd(b)
  }
}


##################################################
# merge data to create a single map per timestep 
##################################################


# STEP10
# must use 'stack' to import all data layers
mypath <- file.path(paste(getwd(), paste0("timebricks", "_t", 1), sep="/"))
tif_maps_t1 <- list.files(mypath, full.names=TRUE)

x <- vector("list", 100)
for(i in 1:100) {
  m <- stack(tif_maps_t1[i])
  x[[i]] <- m
}
t1 <- do.call(merge, x)
plot(t1)












-