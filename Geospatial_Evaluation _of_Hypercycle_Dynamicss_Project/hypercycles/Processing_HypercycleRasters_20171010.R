############################################################
# Title:			             Processing Hypercycle Data
# Project Descriptor:	     AgroEv
# Author:			             bmarron
# Origin Date:		         17 Sept 2017
################################################################


##############
# websites
###############

'''
https://geoscripting-wur.github.io/IntroToRaster/
https://gis.stackexchange.com/questions/60527/how-to-extract-values-from-rasters-at-location-of-points-in-r
https://cran.r-project.org/web/packages/spacetime/vignettes/jss816.pdf
https://informatique-mia.inra.fr/biosp/sites/informatique-mia.inra.fr.biosp-d7/files//analyzing-spatio-temporal.pdf
https://cran.r-project.org/web/views/SpatioTemporal.html
https://cran.r-project.org/web/views/Spatial.html
https://cran.r-project.org/web/views/TimeSeries.html

'''


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%% Hypercycles-to-Time Series Raster Maps %%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


##################
#Libraries
##################
library(deSolve)
library(raster)
library(readr)
library(data.table)
library(rgdal)
library(sp)


##########################
Plot single raster layers
###########################
a<-paste(getwd(), paste0("TimeSeries_t2"), sep="/")
m1<-raster(paste(a, paste0("Oct_07_2017_TimeSeriesMap2.tif"), sep="/"), band=1)
m2<-raster(paste(a, paste0("Oct_07_2017_TimeSeriesMap2.tif"), sep="/"), band=2)
m3<-raster(paste(a, paste0("Oct_07_2017_TimeSeriesMap2.tif"), sep="/"), band=3)
m4<-raster(paste(a, paste0("Oct_07_2017_TimeSeriesMap2.tif"), sep="/"), band=4)
m5<-raster(paste(a, paste0("Oct_07_2017_TimeSeriesMap2.tif"), sep="/"), band=5)


raster::plot(m1)
raster::plot(m2)
raster::plot(m3)
raster::plot(m4)
raster::plot(m5)


######################
Plot as all layers 
using brick()
#######################
a2<-paste(getwd(), paste0("TimeSeries_t2"), sep="/")
brickt2 <-brick(paste(a2, paste0("Oct_07_2017_TimeSeriesMap2.tif"), sep="/"))
raster::plot(brickt2, nc=2)

a3<-paste(getwd(), paste0("TimeSeries_t3"), sep="/")
brickm3 <-brick(paste(a3, paste0("Oct_07_2017_TimeSeriesMap3.tif"), sep="/"))
raster::plot(brickm3, nc=2)

a4<-paste(getwd(), paste0("TimeSeries_t4"), sep="/")
brickm4 <-brick(paste(a4, paste0("Oct_07_2017_TimeSeriesMap4.tif"), sep="/"))
raster::plot(brickm4, nc=2)

a47<-paste(getwd(), paste0("TimeSeries_t47"), sep="/")
brickm4 <-brick(paste(a47, paste0("Oct_07_2017_TimeSeriesMap47.tif"), sep="/"))
raster::plot(brickm4, nc=2)
                

##################
Plot raw data
#################

c2 <- read_csv("~/Desktop/NCSU/2017_DFER_06_GeospatialDataMining/Works_InProgress/Geospatial_Eval _of_Hypercycle_Dynamics/hypercycles/RawData/Oct_07_14_44_35_2017_2.csv", 
                                   col_names = FALSE)


[[2]]
ka1       ka2       ka3       kb1       kb2       kb3       kc1       kc2 
0.9869579 0.3913386 0.7354640 1.5386085 0.3991725 0.1381785 0.2179922 0.1880356 
kc3       kd1       kd2       kd3       ke1       ke2       ke3 
0.1737409 0.7273013 0.4034075 0.3022127 0.6373453 0.7845258 0.9946186 


par(mfrow=c(3,2))
time<-seq(1:501)
plot(c2$X1~time)
plot(c2$X2~time)
plot(c2$X3~time)
plot(c2$X4~time)
plot(c2$X5~time)

par(mfrow=c(1,1))

# ----------------------------------------------------------------------

c3 <- read_csv("~/Desktop/NCSU/2017_DFER_06_GeospatialDataMining/Works_InProgress/Geospatial_Eval _of_Hypercycle_Dynamics/hypercycles/RawData/Oct_07_14_44_36_2017_3.csv", 
               col_names = FALSE)

[[3]]
ka1        ka2        ka3        kb1        kb2        kb3        kc1 
0.35338272 0.08816453 1.10813510 0.26727043 0.04849966 0.28449084 0.72561371 
kc2        kc3        kd1        kd2        kd3        ke1        ke2 
2.30812321 0.53069768 0.10606029 0.34399637 0.40090048 0.03000603 0.19188286 
ke3 
0.01100272 

time<-seq(1:501)
plot(c3$X1~time)
plot(c3$X2~time)
plot(c3$X3~time)
plot(c3$X4~time)
plot(c3$X5~time)

#------------------------------------------------
c4 <- read_csv("~/Desktop/NCSU/2017_DFER_06_GeospatialDataMining/Works_InProgress/Geospatial_Eval _of_Hypercycle_Dynamics/hypercycles/RawData/Oct_07_14_44_37_2017_4.csv", 
               col_names = FALSE)

[[4]]
ka1        ka2        ka3        kb1        kb2        kb3        kc1 
0.27087559 1.01285406 0.67677843 1.23256957 1.37302803 0.51992960 0.01116306 
kc2        kc3        kd1        kd2        kd3        ke1        ke2 
1.06793810 0.37469391 0.92309504 0.86285628 1.38386400 0.76391577 0.16455464 
ke3 

time<-seq(1:501)
plot(c4$X1~time)
plot(c4$X2~time)
plot(c4$X3~time)
plot(c4$X4~time)
plot(c4$X5~time)            


#---------------------------------------------------------------
c47 <- read_csv("~/Desktop/NCSU/2017_DFER_06_GeospatialDataMining/Works_InProgress/Geospatial_Eval _of_Hypercycle_Dynamics/hypercycles/RawData/Oct_07_14_45_55_2017_74.csv", 
               col_names = FALSE)

[[47]]
ka1        ka2        ka3        kb1        kb2        kb3        kc1 
0.49280284 0.03370478 0.54509658 0.43538533 0.16694514 0.80047158 1.81412344 
kc2        kc3        kd1        kd2        kd3        ke1        ke2 
1.76697621 2.47574192 1.16754038 0.95275677 0.43154429 0.42320746 0.20582169 
ke3 
1.23272531


time<-seq(1:501)
plot(c47$X1~time)
plot(c47$X2~time)
plot(c47X3~time)
plot(c47$X4~time)
plot(c47$X5~time)



#########################
extract data from maps
(multi-band rasters)
#########################


a2<-paste(getwd(), paste0("TimeSeries_t2"), sep="/")
brickt2 <-brick(paste(a2, paste0("Oct_07_2017_TimeSeriesMap2.tif"), sep="/"))
brickt2.data <-data.frame(coordinates(brickt2), extract(brickt2, seq(1:100), layer=1, nl=1))
names(brickt2.data) <- c("x", "y", "a @ t2")
brickt2.data

# START HERE
re-run everything
put TimeSeries maps al in one folder
loop thru timeseries maps and pull out layer=1 data
cbind to 


#loop this
a<-paste(getwd(), paste0("TimeSeries"), sep="/")
mypath2 <- file.path(paste(b, paste0("timebricks", "_t", j), sep="/"))

myfile <- paste0(mytime, "_TimeSeriesMap", j, ".tif")

tif_maps <- list.files(mypath1, full.names=TRUE)
timebrick <-brick(paste(a2, paste0("Oct_07_2017_TimeSeriesMap2.tif"), sep="/"))
brickt2.data <-data.frame(coordinates(brickt2), extract(brickt2, seq(1:100), layer=1, nl=1))
names(brickt2.data) <- c("x", "y", "a @ t2")
brickt2.data



