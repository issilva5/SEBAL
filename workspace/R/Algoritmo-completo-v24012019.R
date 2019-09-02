######################################################################################
#                                                                                      #
#                         EU BRAZIL Cloud Connect                                      #
#                                                                                      #
#                                                                                      #
########################################################################################

options(echo=TRUE)
options(digits=22)
options(error=traceback)
rm(list=ls())

library(R.utils)
library(raster)
library(rgdal)
library(maptools)
library(ncdf4)
library(sp)
library(snow)
library(here)

args = commandArgs(trailingOnly=TRUE)
WD <- args[1]
setwd(WD) # Working Directory

# Changing raster tmpdir
rasterOptions(tmpdir=args[2])

# Load the source code in landsat.R to this code
source(here("workspace/R", "landsat2.R"))

# File that stores the Image Directories (TIFs, MTL, FMask)
dados <- read.csv(here("workspace/R", "dados.csv"), sep=";", stringsAsFactors=FALSE)

auxProc1 <- proc.time()

#################################### Constants ##########################################

k <- 0.41		# Von Karman
g <- 9.81		# Gravity
rho <- 1.15		# Air density
cp <- 1004		# Specific heat of air
Gsc <- 0.082		# Solar constant (0.0820 MJ m-2 min-1)
clusters <- 7		# Number of clusters used in image processing - some raster library methods are naturally coded to run in a clustered way

######################### Reading sensor parameters #####################################

p.s.TM1 <- read.csv(here("workspace/R", "parametros do sensor/parametrosdosensorTM1.csv"), sep=";", stringsAsFactors=FALSE)
p.s.TM2 <- read.csv(here("workspace/R", "parametros do sensor/parametrosdosensorTM2.csv"), sep=";", stringsAsFactors=FALSE)
p.s.ETM <- read.csv(here("workspace/R", "parametros do sensor/parametrosdosensorETM.csv"), sep=";", stringsAsFactors=FALSE)
p.s.LC <- read.csv(here("workspace/R", "parametros do sensor/parametrosdosensorLC.csv"), sep=";", stringsAsFactors=FALSE)

# Read relative distance from Sun to Earth
load(here("workspace/R", "d_sun_earth.RData"))

# Set projection and spatial resolution
WGS84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84"

######################### Image Information ######################################

fic.dir <- dados$File.images[1]  # Images Directory
print(fic.dir)
MTL <- read.table(dados$MTL[1], skip=0, nrows=140, sep="=", quote = "''", as.is=TRUE) # Reading MTL File

fic <- substr(MTL$V2[MTL$V1 == grep(pattern="LANDSAT_SCENE_ID", MTL$V1, value=T)], 3, 23)

n.sensor <- as.numeric(substr(fic, 3, 3)) # Sensor Number

if (n.sensor==8) MTL <- read.table(dados$MTL[1], skip=0, nrows=-1, sep="=", quote="''", as.is=TRUE, fill=TRUE) # Reading MTL File for Sensor number 8

WRSPR <- substr(fic, 4, 9)						#WRSPR
Ano <- as.numeric(substr(fic, 10, 13))			#Images year
Dia.juliano <- as.numeric(substr(fic, 14, 16))	#Julian Day

# Getting the sum elevation at the time of Image Capture
sun_elevation <- as.numeric(MTL$V2[MTL$V1 == grep(pattern="SUN_ELEVATION", MTL$V1, value=TRUE)])
costheta <- sin(sun_elevation*pi/180) # From SUN ELEVATION

# Setting the sensor parameter by the Sattelite sensor type and data
if (n.sensor==8) p.s <- p.s.LC
if (n.sensor==7) p.s <- p.s.ETM
if (Ano < 1992 & n.sensor==5) p.s <- p.s.TM1 
if (Ano > 1992 & n.sensor==5) p.s <- p.s.TM2

# Time image
acquired_date <- as.Date(MTL$V2[MTL$V1==grep(pattern="DATE_ACQUIRED", MTL$V1, value=TRUE)])
daysSince1970 <- as.numeric(acquired_date)
tdim <- ncdim_def("time", "days since 1970-1-1", daysSince1970, unlim=TRUE, create_dimvar=TRUE, "standard", "time")

# Reading image file
# The Images are of the type ".tif" that represents each spectral band captured by the satellite
# Depending on the sattelite the number of spectral bands captured are differents
fichs.imagens <- list.files(path=fic.dir, pattern="*.TIF", full.names=TRUE)

getBandsPath <- function(n.sensor){
  wanted_bands <- NULL
  if(n.sensor == 8) {
    wanted_bands <- c("B2", "B3", "B4", "B5", "B6", "B7", "B10")
  } else if(n.sensor == 7) {
    wanted_bands <- c("B1", "B2", "B3", "B4", "B5", "B6_VCID_1", "B6_VCID_2", "B7")
  } else if(n.sensor == 5) {
    wanted_bands <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7")
  }
  
  bands_path <- list()
  for (i in 1:length(wanted_bands)) {
    for (j in 1:length(fichs.imagens)) {
      if(regexpr(paste(wanted_bands[i], '.TIF', sep=""), fichs.imagens[[j]]) != -1) {
        bands_path[[i]] <- fichs.imagens[[j]]
      }
    }
  }

  return(bands_path)
}

#proc.time()

# Reading
fic.st <- stack(as.list(getBandsPath(n.sensor)))

#print(proc.time())

################################# Fmask ###########################################

# Identifier of clouds and shadows
# The FMask serves to identify if exist any cloud or shadow on the image, the existence of clouds or shadow disturbs the results.

n.fmask <- length(fichs.imagens)
Fmask <- raster(fichs.imagens[[n.fmask]])
fmask <- as.vector(Fmask)

if (n.sensor != 8) mask_filter <- 672 else 
  mask_filter <- 2720
for (i in 1:nlayers(fic.st)) {
  f <- fic.st[[i]][]
  f[fmask != mask_filter] <- NaN
  fic.st[[i]][] <- f 
}

#print(proc.time())

if (0.99<=(sum(is.na(values(fic.st)))/7)/(fic.st@ncols*fic.st@nrows)) { 
	print("Imagem incompativel para o processamento,mais de 99% nuvem e sombra de nuvem")
	quit("no", 1, FALSE)
}

#print(proc.time())

# Changing the projection of the images (UTM to GEO)
# This operation can be done in a parallel way by Clusters, projectRaster is implemented to naturally be executed by clusters
# The number of used clusters is given by the 'clusters' constant

#beginCluster(clusters)
	fic.st <- projectRaster(fic.st, crs=WGS84)
#endCluster()

#print(proc.time())

# Reading Bounding Box
# The Bounding Box area that is important and has less noise in the Image
fic.bounding.boxes <- paste(here("workspace/R", "wrs2_asc_desc/wrs2_asc_desc.shp"))
BoundingBoxes <- readShapePoly(fic.bounding.boxes, proj4string=CRS(WGS84))
BoundingBox <- BoundingBoxes[BoundingBoxes@data$WRSPR == WRSPR, ]

# Reading Elevation
# Read the File that stores the Elevation of the image area, this influence on some calculations
fic.elevation <- paste(here("workspace/R", "Elevation/srtm_29_14.tif"))
raster.elevation <- raster(fic.elevation)

print("Raster elevation before crop")
print(raster.elevation)

raster.elevation <- crop(raster.elevation, extent(BoundingBox))

print("Raster elevation after crop")
print(raster.elevation)

#print(proc.time())

# Setting the raster elevation resolution as equals to the Fmask raster resolution
raster.elevation.aux <- raster(raster.elevation)

print("Raster elevation aux before resolution exchange for images band")
print(raster.elevation.aux)

res(raster.elevation.aux) <- res(fic.st) # The raster elevation aux resolution is the same of raster fmask

print("Raster elevation aux after res")
print(raster.elevation.aux)

# Resample images
beginCluster(clusters) # ?? beginClusters, whereis endClusters??
	raster.elevation <- resample(raster.elevation, raster.elevation.aux, method="ngb")
endCluster()

print("Raster elevation after resampling")
print(raster.elevation)

print("Image band before resampling")
print(fic.st[[1]])

#print(proc.time())

#################### Resampling satellite bands images #####################################

# This block of code resample the image based on the Elevation of the terrain captured by the sat
# The Elevation of the terrain needs to be taken into account
# This block is already Clustered

# See if timeouts presented here will be the default or distinct between sites
# timeout before = 2177.062
# timeout now is 3600 (cause: Azure slowness)

image.rec <- NULL;
imageResample <- function() {
  beginCluster(clusters)
	image_resample <- resample(fic.st, raster.elevation, method="ngb")
  endCluster()
  return(image_resample)
}

res <- NULL;
tryCatch({
  res <- withTimeout({
    image.rec <- imageResample();
  }, timeout=3600);
}, TimeoutException=function(ex) {
  cat("Image resample timedout. Exiting with 124 code...\n");
  quit("no", 124, FALSE)
})
#names(fic.st) <- c("B2", "B3", "B4", "B5", "B6", "B7", "B10")
#output.path <- dados$Path.Output[1]
#writeRaster(fic.st, output.path, overwrite=TRUE, format="CDF", varname = fic, varunit = "daily", longname = fic, xname = "lon", yname = "lat", bylayer=TRUE, suffix="names")
############################################################################################
#quit()
print("Image band after resampling")
print(fic.st[[1]])

#print(proc.time())

# Reading file Station weather
fic.sw <- dados$File.Station.Weather[1]
table.sw <- (read.csv(fic.sw, sep=";", header=FALSE, stringsAsFactors=FALSE))

# Transmissivity 
tal <- 0.75+2*10^-5*raster.elevation

print(proc.time())
#print("Fim da fase 1 - pre-processamento")

################## Phase 1: Calculating the image energy balance ##################################

# This block calculate the energy balance of the image

output <- NULL;
outputLandsat <- function() {
  output <- landsat()
  return(output)
}

# timeout before = 2665.151
# timeout now is 7200 (cause: Azure slowness)
res <- NULL;
tryCatch({
  res <- withTimeout({
    output <- outputLandsat(); print("sai da landasat3.R");
  }, timeout=7200);
}, TimeoutException=function(ex) {
  cat("Output landsat timedout. Exiting with 124 code...\n");
  quit("no", 124, FALSE)
})

###########################################################################################

#auxProc12 <- proc.time()
#print(auxProc12)
################## Masking landsat rasters output #########################################

# This block mask the values in the landsat output rasters that has cloud cells and are inside the Bounding Box required
# This block is already Clustered

outputMask <- function() {
  output <- mask(output, BoundingBox)
  return(output)
}

# timeout before = 1716.853
# timeout now is 10800 (cause: Azure slowness)

res <- NULL;
tryCatch({
 res <- withTimeout({
    output <- outputMask();
  }, timeout=10800);
}, TimeoutException=function(ex) {
  cat("Output Fmask timedout. Exiting with 124 code...\n");
  quit("no", 124, FALSE)
})

##########################################################################################

#print(proc.time())

################## Write to files landsat output rasters #################################

# This block write landsat outputs rasters to files

output.path<-paste(dados$Path.Output[1], "/", fic, ".nc", sep="")
outputWriteRaster <- function() {
  names(output) <- c("Rn", "TS", "NDVI", "EVI", "LAI", "G", "alb")
  writeRaster(output, output.path, overwrite=TRUE, format="CDF", varname=fic, varunit="daily", longname=fic, xname="lon", yname="lat", bylayer=TRUE, suffix="names")
}

# timeout before = 1708.507
# timeout now is 10800 (cause: Azure slowness)

res <- NULL;
tryCatch({
  res <- withTimeout({
    outputWriteRaster();
  }, timeout=10800);
}, TimeoutException=function(ex) {
  cat("Output write raster timedout. Exiting with 124 code...\n");
  quit("no", 124, FALSE)
})

#proc.time()

#########################################################################################

variablesNames <- c("alb", "EVI", "G", "LAI", "NDVI", "Rn", "TS")
dimLatDef <- NULL
dimLonDef <- NULL

for(i in 1:length(variablesNames)){

	end_file_name <- paste("_", variablesNames[i], ".nc", sep="")
	#proc.time()

	#Opening old NetCDF
	var_output <- paste(dados$Path.Output[1], "/", fic, end_file_name, sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

	#proc.time()

	if(i == 1){
		# Getting lat and lon values from old NetCDF
		oldLat <- ncvar_get(nc, "lat", start=1, count=raster.elevation@nrows)
		oldLon <- ncvar_get(nc, "lon", start=1, count=raster.elevation@ncols)

		# Defining latitude and longitude dimensions
		dimLatDef <- ncdim_def("lat", "degrees", oldLat, unlim=FALSE, longname="latitude")
		dimLonDef <- ncdim_def("lon", "degrees", oldLon, unlim=FALSE, longname="longitude")
	}

	#proc.time()

	#New file
	file_output <- paste(dados$Path.Output[1], "/", fic, end_file_name, sep="")
	oldValues <- ncvar_get(nc, fic)
	newValues <- ncvar_def(variablesNames[i], "daily", list(dimLonDef, dimLatDef, tdim), longname=variablesNames[i], missval=NaN, prec="double")
	nc_close(nc)
	newNCDF4 <- nc_create(file_output, newValues)
	ncvar_put(newNCDF4, variablesNames[i], oldValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newNCDF4)

	#proc.time()
}

print(proc.time())
#print("Fim da fase 2 - Calculo do Rn e G")

################################################################################################

	#CONSTANTES
	# Weather station data
    x<-3                                    # Wind speed sensor Height (meters)
    hc<-0.2                                 # Vegetation height (meters)
    Lat<-table.sw$V4[1]     # Station Latitude
    Long<-table.sw$V5[1]    # Station Longitude

    # Surface roughness parameters in station
    azom <- -3              #Parameter for the Zom image
    bzom <- 6.47    #Parameter for the Zom image
    F_int <- 0.16   #internalization factor for Rs 24 calculation (default value)

print(trunc)


hotPixelSelection <- function(Rn, G, TS, NDVI){

	HO<-Rn[]-G[] # Read as a Vector
	
	x<-TS[][(NDVI[]>0.15 &!is.na(NDVI[]))  & (NDVI[]<0.20 &!is.na(NDVI[])) ] # Returns a vector
	x<-x[x>273.16]
	TS.c.hot<-sort(x)[round(0.95*length(x))] # Returns one value
	HO.c.hot<-HO[][(NDVI[]>0.15 &!is.na(NDVI[])) & (NDVI[]<0.20 &!is.na(NDVI[])) & TS[]==TS.c.hot] # Returns one value
	print("NDVI ente 0.15 e 0.20")
	print(length(x))
	print("Temperatura escolhida")
	print(TS.c.hot)
	print("HO dos que tem essa temp")
	print(HO.c.hot)
	if (length(HO.c.hot)==1){
		ll.hot<-which(TS[]==TS.c.hot & HO[]==HO.c.hot)
		xy.hot <- xyFromCell(TS, ll.hot)
		ll.hot.f<-cbind(as.vector(xy.hot[1,1]), as.vector(xy.hot[1,2]))
	  }else{
		HO.c.hot.min<-sort(HO.c.hot)[ceiling(0.25*length(HO.c.hot))]
		HO.c.hot.max<-sort(HO.c.hot)[ceiling(0.75*length(HO.c.hot))]
		ll.hot<-which(TS[]==TS.c.hot & HO[]>HO.c.hot.min & HO[]<HO.c.hot.max)
                xy.hot <- xyFromCell(TS, ll.hot)
		print("Localizacao dos pixels")
		print(xFromCol(TS, xy.hot[,1]))
		print(yFromRow(TS, xy.hot[,2]))
		print(xy.hot)
		NDVI.hot<-extract(NDVI,xy.hot, buffer=105)
		NDVI.hot.2<-NDVI.hot[!sapply(NDVI.hot, is.null)]
		print("Extract result")
		print(NDVI.hot.2)
		NDVI.hot.cv <- sapply(NDVI.hot.2,sd, na.rm=TRUE)/sapply(NDVI.hot.2, mean, na.rm=TRUE)
                print("CVs")
		print(NDVI.hot.cv)
		i.NDVI.hot.cv<-which.min(NDVI.hot.cv)
		ll.hot.f<-cbind(as.vector(xy.hot[i.NDVI.hot.cv,1]), as.vector(xy.hot[i.NDVI.hot.cv,2]))
	}

	return(ll.hot.f)
}

coldPixelSelection <- function(Rn, G, TS, NDVI){
		
	HO<-Rn[]-G[] # Read as a Vector
	
	x<-TS[][(NDVI[]<0 &!is.na(NDVI[])) & !is.na(HO)]
	x<-x[x>273.16]
	
	TS.c.cold<-sort(x)[round(0.5*length(x))]
	
	HO.c.cold<-HO[(NDVI[]<0 & !is.na(NDVI[])) & TS[]==TS.c.cold & !is.na(HO)]
	print("NDVI menor que 0")
	print(length(x))
	print("Temperatura escolhida")
	print(TS.c.cold)
	print("HOs dos que tem essa temp")
	print(HO.c.cold)

	if (length(HO.c.cold)==1){
		ll.cold<-which(TS[]==TS.c.cold & HO==HO.c.cold)
		xy.cold <- xyFromCell(TS, ll.cold)
		ll.cold.f<-cbind(as.vector(xy.cold[1,1]), as.vector(xy.cold[1,2]))
	}else{
		HO.c.cold.min<-sort(HO.c.cold)[ceiling(0.25*length(HO.c.cold))]
		HO.c.cold.max<-sort(HO.c.cold)[ceiling(0.75*length(HO.c.cold))]
		
		ll.cold<-which(TS[]==TS.c.cold & (HO>HO.c.cold.min &!is.na(HO)) & (HO<HO.c.cold.max & !is.na(HO)))
		print(ll.cold)
		print("Localizacao")
		xy.cold <- xyFromCell(TS, ll.cold)
		print(xFromCol(TS, xy.cold[,1]))
		print(yFromRow(TS, xy.cold[,2]))
		print(xy.cold)
		NDVI.cold<-extract(NDVI,xy.cold, buffer=105)
		print("Extract result")
		print(NDVI.cold)
		NDVI.cold.2<-NDVI.cold[!sapply(NDVI.cold, is.null)]
		
		# Maximum number of neighboring pixels with $NVDI < 0$
		t<-function(x){ sum(x<0,na.rm = TRUE)}
		n.neg.NDVI<-sapply(NDVI.cold.2,t)
		print(n.neg.NDVI)
		i.NDVI.cold<-which.max(n.neg.NDVI)
		
		ll.cold.f<-cbind(as.vector(xy.cold[i.NDVI.cold,1]), as.vector(xy.cold[i.NDVI.cold,2]))
	}

	return(ll.cold.f)
}

windVelocity200 <- function(){
	#VER A DECLARAÇÃO DAS CONSTANTES HC, K, X
	zom.est <- hc*0.12

	# Friction velocity at the station (ustar.est)
	ustar.est <- k*table.sw$V6[2]/log((x)/zom.est)
	
	# Velocity 200 meters
	u200 <- trunc(ustar.est/k*log(200/zom.est), 4)

	return(u200)
}

phase2 <- function() {
	####################### Selection of reference pixels ###################################
	
	# Getting the rasters output
	Rn<-output[[1]]
	TS<-output[[2]]
	NDVI<-output[[3]]
	EVI<-output[[4]]
	LAI<-output[[5]]
	G<-output[[6]]
	alb<-output[[7]]
	
	ll.hot.f <- hotPixelSelection(Rn, G, TS, NDVI)
	ll.cold.f <- coldPixelSelection(Rn, G, TS, NDVI)
	
	# Location of reference pixels (hot and cold)
	ll_ref<-rbind(ll.hot.f[1,],ll.cold.f[1,])
	colnames(ll_ref)<-c("long", "lat")
	rownames(ll_ref)<-c("hot","cold")
	#print(proc.time())
	print(ll_ref)
	####################################################################################
	
	# Weather station data
#	x<-3 					# Wind speed sensor Height (meters)
#	hc<-0.2 				# Vegetation height (meters)
#	Lat<-table.sw$V4[1] 	# Station Latitude
#	Long<-table.sw$V5[1] 	# Station Longitude
	
	# Surface roughness parameters in station
#	azom <- -3		#Parameter for the Zom image
#	bzom <- 6.47	#Parameter for the Zom image
#	F_int <- 0.16	#internalization factor for Rs 24 calculation (default value)
	
	#print(proc.time())
	
	# Velocity 200 meters
	u200 <- windVelocity200()
	
	# Zom for all pixels
	#zom<-exp(azom+bzom*NDVI[]) # Changed from Raster to Vector
	zom <- NDVI
	zom[] <- exp(azom+bzom*NDVI[])
	zom[] <- trunc(zom[], 4)
	print("Passou do zom")
	
	# Initial values
	ustar<-NDVI
	ustar[]<-k*u200/(log(200/zom[]))		# Friction velocity for all pixels #RASTER - VETOR 
	ustar[] <- trunc(ustar[], 4)
	
	ustar_before <- NDVI
	ustar_before[] <- k*u200/(log(200/zom[])) #FIXME:
	ustar_before[] <- trunc(ustar_before[], 4)
	print("Passou ustar")
	
	rah<-NDVI
	rah[]<-(log(2/0.1))/(ustar[]*k) 		# Aerodynamic resistance for all pixels #RASTER - VETOR
	rah[] <- trunc(rah[], 4)
	
	rah_before <- NDVI
	rah_before[] <-(log(2/0.1))/(ustar[]*k)  #FIXME:
	rah_before[] <- trunc(rah_before[], 4)
	print("Passou rah")

	base_ref<-stack(NDVI,TS,Rn,G,ustar,rah) # Raster
	nbase<-c("NDVI","TS","Rn","G")
	names(base_ref)<-c(nbase,"ustar","rah")
	
	value.pixels.ref<-extract(base_ref,ll_ref)
	rownames(value.pixels.ref)<-c("hot","cold")
	H.hot<-value.pixels.ref["hot","Rn"]-value.pixels.ref["hot","G"]
	value.pixel.rah<-value.pixels.ref["hot","rah"]
	print(value.pixels.ref)

	i<-1
	Erro<-TRUE
	print("Before rah correct")
	#print(proc.time())
	# Beginning of the cycle stability
	while(Erro){
	  rah.hot.0<-value.pixel.rah[i] # Value
	  print("rah.hot.0")
	  # Hot and cold pixels      
	  dt.hot<-H.hot*rah.hot.0/(rho*cp) # Value
	  b<-dt.hot/(value.pixels.ref["hot","TS"]-value.pixels.ref["cold","TS"]) # Value
	  a<- -b*(value.pixels.ref["cold","TS"]-273.15) # Value
	  print("H.hot")
	  print(H.hot)
	  print("dt.hot")
	  print(dt.hot)
	  print("a")
	  print(a)
	  print("b")
	  print(b)
    	  print("before H")
	  # All pixels
	  H<-rho*cp*(a+b*(TS[]-273.15))/rah[] # Changed from Raster to Vector
	  H <- trunc(H, 4)
	  L<- -1*((rho*cp*ustar[]^3*TS[])/(k*g*H)) # Changed from Raster to Vector
	  L <- trunc(L, 4)
	  y_0.1<-(1-16*0.1/L)^0.25 # Changed from Raster to Vector
	  y_0.1 <- trunc(y_0.1, 4)
	  y_2<-(1-16*2/L)^0.25 # Changed from Raster to Vector
	  y_2 <- trunc(y_2, 4)
	  x200<-(1-16*200/L)^0.25 # Changed from Raster to Vector
	  x200 <- trunc(x200, 4)
	  psi_0.1<-2*log((1+y_0.1^2)/2) # Changed from Raster to Vector
	  psi_0.1[L>0 &!is.na(L)]<--5*(0.1/L[L>0 &!is.na(L)]) # Changed from Raster to Vector
	  psi_0.1 <- trunc(psi_0.1, 4)
	  psi_2<-2*log((1+y_2^2)/2)  # Changed from Raster to Vector
	  psi_2[L>0 &!is.na(L) ]<--5*(2/L[L>0 &!is.na(L)]) # Changed from Raster to Vector
	  psi_2 <- trunc(psi_2, 4)
	  psi_200<-2*log((1+x200)/2)+log((1+x200^2)/2)-2*atan(x200)+0.5*pi # Changed from Raster to Vector
	  psi_200[L>0 &!is.na(L) ]<--5*(2/L[(L>0 &!is.na(L))]) # Changed from Raster to Vector
	  psi_200 <- trunc(psi_200, 4)
	  ustar<-k*u200/(log(200/zom[])-psi_200) # Changed from Raster to Vector # Friction velocity for all pixels
	  ustar <- trunc(ustar, 4)
	  rah<-NDVI
	  rah[]<-(log(2/0.1)-psi_2+psi_0.1)/(ustar*k) # Changed from Raster to Vector # Aerodynamic resistency for all pixels
	  rah[] <- trunc(rah[], 4)
	  rah.hot<-extract(rah,matrix(ll_ref["hot",],1,2)) # Value
	  value.pixel.rah<-c(value.pixel.rah,rah.hot) # Value
	  print("after rah")
	  print(i)
	  print(value.pixel.rah)
	  Erro<-(abs(1-rah.hot.0/rah.hot)>=0.05)
	  print(abs(1-rah.hot.0/rah.hot))
	  i<-i+1
	}
	L_final <- NDVI
	L_final[] <- L

	y_0.1_final <- NDVI
	y_0.1_final[] <- y_0.1

	y_2_final <- NDVI
	y_2_final[] <- y_2

	x200_final <- NDVI
	x200_final[] <- x200

	psi_0.1_final <- NDVI
	psi_0.1_final[] <- psi_0.1

	psi_2_final <- NDVI
	psi_2_final[] <- psi_2

	psi_200_final <- NDVI
	psi_200_final[] <- psi_200


	#print(proc.time())
	print("After rah correct")
	# End sensible heat flux (H)
	
	# Hot and cold pixels
	dt.hot<-H.hot*rah.hot/(rho*cp)                  
	b<-dt.hot/(value.pixels.ref["hot","TS"]-value.pixels.ref["cold","TS"]) 
	a<- -b*(value.pixels.ref["cold","TS"]-273.15)                          
	
	#print(proc.time())
	
	# All pixels
	H <-rho*cp*(a+b*(TS[]-273.15))/rah[] # Vector 
	H[(H>(Rn[]-G[]) &!is.na(H))]<-(Rn[]-G[])[(H>(Rn[]-G[]) &!is.na(H))] # Vector
	H <- trunc(H, 4)
	
	ustar_final <- NDVI
	ustar_final[] <- ustar
	print("Passou ustar final")

	#rah_final <- NDVI
	#rah_final[] <- rah
	#print("Passou rah final")

	H_final <- NDVI
	#H_final[] <-rho*cp*(a+b*(TS[]-273.15))/rah[] # Vector 
	#H_final[(H>(Rn[]-G[]) &!is.na(H))]<-(Rn[]-G[])[(H>(Rn[]-G[]) &!is.na(H))] # Vector
	H_final[] <- H
	print("Passou H final")

	#print(proc.time())

	# Instant latent heat flux (LE)
	LE<-Rn[]-G[]-H
	LE_final <- NDVI
	LE_final[] <- trunc(Rn[]-G[]-H, 4)
	print("Passou LE")
	# Upscalling temporal
	dr<-(1/d_sun_earth$dist[Dia.juliano])^2 		# Inverse square of the distance on Earth-SOL
	sigma<-0.409*sin(((2*pi/365)*Dia.juliano)-1.39) # Declination Solar (rad)
	phi<-(pi/180)*Lat 								# Solar latitude in degrees
	omegas<-acos(-tan(phi)*tan(sigma)) 				# Angle Time for sunsets (rad)
	Ra24h<-(((24*60/pi)*Gsc*dr)*(omegas*sin(phi)*
	        sin(sigma)+cos(phi)*cos(sigma)*sin(omegas)))*(1000000/86400)
	
	#print(proc.time())
	
	# Short wave radiation incident in 24 hours (Rs24h)
	Rs24h<-F_int*sqrt(max(table.sw$V7[])-min(table.sw$V7[]))*Ra24h
	Rs24h <- trunc(Rs24h, 4)
	
	FL<-110                                
	Rn24h_dB<-(1-alb[])*Rs24h-FL*Rs24h/Ra24h		# Method of Bruin #VETOR
	Rn24h_dB_final <- NDVI
	Rn24h_dB_final[] <-(1-alb[])*Rs24h-FL*Rs24h/Ra24h
	Rn24h_dB_final <- trunc(Rn24h_dB_final, 4)
	print("Passou Rn24h")
	# Evapotranspiration fraction Bastiaanssen
	EF<-NDVI
	EF[]<-LE/(Rn[]-G[])
	EF[] <- trunc(EF[], 4)
	
	# Sensible heat flux 24 hours (H24h), nao foi usado
	# H24h_dB<-(1-EF[])*Rn24h_dB
	
	# Latent Heat Flux 24 hours (LE24h)
	LE24h_dB<-EF[]*Rn24h_dB
	LE24h_dB_final <- NDVI
	LE24h_dB_final[] <- EF[]*Rn24h_dB
	LE24h_dB_final <- trunc(LE24h_dB_final, 4)
	print("Passou LE24h")
	
	# Evapotranspiration 24 hours (ET24h)
	ET24h_dB<-NDVI
	ET24h_dB[]<-LE24h_dB*86400/((2.501-0.00236* (max(table.sw$V7[])+min(table.sw$V7[]))/2)*10^6)
	ET24h_dB[] <- trunc(ET24h_dB[], 4)
	print("Antes de salvar")
	#print(proc.time())
	
	print(zom)
	print(ustar_before)
	print(rah_before)
	print(ustar_final)
	print(rah)
	print(H_final)
	print(LE_final)
	print(Rn24h_dB_final)
	print(LE24h_dB_final)
	print(EF)
	print(ET24h_dB)

	output.evapo<-stack(zom, ustar_before, rah_before, ustar_final, rah, H_final, LE_final, Rn24h_dB_final, LE24h_dB_final, EF, ET24h_dB, L_final, y_0.1_final, y_2_final, x200_final, psi_0.1_final, psi_2_final, psi_200_final)
	print(nlayers(output.evapo))
	print(output.evapo)
	output.names <- c('zom', 'ustar', 'Rah', 'ustar_after', 'Rah_after', 'H', 'LatentHF', 'Rn24h', 'LatentHF24h', 'EF', 'ET24h', 'L', 'y01', 'y2', 'x200', 'psi01', 'psi2', 'psi200')
	print(length(output.names))
	print(output.names)
	names(output.evapo) <- output.names
	writeRaster(output.evapo, output.path, overwrite=TRUE, format="CDF", varname=fic, varunit="daily", longname=fic, xname="lon", yname="lat", bylayer=TRUE, suffix="names")
	print("Depois do writeRaster")
	#print(proc.time())
	
	# Opening old EF NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_EF.nc",sep="")
	nc<-nc_open(var_output, write=TRUE,readunlim=FALSE,verbose=TRUE,auto_GMT=FALSE,suppress_dimvals=FALSE)
	
	# New EF file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_EF.nc",sep="")
	oldEFValues<-ncvar_get(nc,fic)
	newEFValues<-ncvar_def("EF","daily",list(dimLonDef,dimLatDef,tdim),longname="EF",missval=NaN,prec="double")
	nc_close(nc)
	newEFNCDF4<-nc_create(file_output,newEFValues)
	ncvar_put(newEFNCDF4,"EF",oldEFValues,start=c(1,1,1),count=c(raster.elevation@ncols,raster.elevation@nrows,1))
	nc_close(newEFNCDF4)
	
	#print(proc.time())
	
	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_ET24h.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New ET24h file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_ET24h.nc",sep="")
	oldET24hValues<-ncvar_get(nc,fic)
	newET24hValues<-ncvar_def("ET24h","daily", list(dimLonDef, dimLatDef, tdim), longname="ET24h", missval=NaN, prec="double")
	nc_close(nc)
	newET24hNCDF4<-nc_create(file_output,newET24hValues)
	ncvar_put(newET24hNCDF4, "ET24h", oldET24hValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newET24hNCDF4)
	
	#print(proc.time())
	#FIXME:
	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_zom.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New zom file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_zom.nc",sep="")
	oldzomValues<-ncvar_get(nc,fic)
	newzomValues<-ncvar_def("zom","daily", list(dimLonDef, dimLatDef, tdim), longname="zom", missval=NaN, prec="double")
	nc_close(nc)
	newzomNCDF4<-nc_create(file_output,newzomValues)
	ncvar_put(newzomNCDF4, "zom", oldzomValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newzomNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_ustar.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New ustar file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_ustar.nc",sep="")
	oldustarValues<-ncvar_get(nc,fic)
	newustarValues<-ncvar_def("ustar","daily", list(dimLonDef, dimLatDef, tdim), longname="ustar", missval=NaN, prec="double")
	nc_close(nc)
	newustarNCDF4<-nc_create(file_output,newustarValues)
	ncvar_put(newustarNCDF4, "ustar", oldustarValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newustarNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_Rah.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New Rah file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_Rah.nc",sep="")
	oldRahValues<-ncvar_get(nc,fic)
	newRahValues<-ncvar_def("Rah","daily", list(dimLonDef, dimLatDef, tdim), longname="Rah", missval=NaN, prec="double")
	nc_close(nc)
	newRahNCDF4<-nc_create(file_output,newRahValues)
	ncvar_put(newRahNCDF4, "Rah", oldRahValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newRahNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_ustar_after.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New ustar_after file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_ustar_after.nc",sep="")
	oldustar_afterValues<-ncvar_get(nc,fic)
	newustar_afterValues<-ncvar_def("ustar_after","daily", list(dimLonDef, dimLatDef, tdim), longname="ustar_after", missval=NaN, prec="double")
	nc_close(nc)
	newustar_afterNCDF4<-nc_create(file_output,newustar_afterValues)
	ncvar_put(newustar_afterNCDF4, "ustar_after", oldustar_afterValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newustar_afterNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_Rah_after.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New Rah_after file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_Rah_after.nc",sep="")
	oldRah_afterValues<-ncvar_get(nc,fic)
	newRah_afterValues<-ncvar_def("Rah_after","daily", list(dimLonDef, dimLatDef, tdim), longname="Rah_after", missval=NaN, prec="double")
	nc_close(nc)
	newRah_afterNCDF4<-nc_create(file_output,newRah_afterValues)
	ncvar_put(newRah_afterNCDF4, "Rah_after", oldRah_afterValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newRah_afterNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_H.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New H file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_H.nc",sep="")
	oldHValues<-ncvar_get(nc,fic)
	newHValues<-ncvar_def("H","daily", list(dimLonDef, dimLatDef, tdim), longname="H", missval=NaN, prec="double")
	nc_close(nc)
	newHNCDF4<-nc_create(file_output,newHValues)
	ncvar_put(newHNCDF4, "H", oldHValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newHNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_LatentHF.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New LatentHF file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_LatentHF.nc",sep="")
	oldLatentHFValues<-ncvar_get(nc,fic)
	newLatentHFValues<-ncvar_def("LatentHF","daily", list(dimLonDef, dimLatDef, tdim), longname="LatentHF", missval=NaN, prec="double")
	nc_close(nc)
	newLatentHFNCDF4<-nc_create(file_output,newLatentHFValues)
	ncvar_put(newLatentHFNCDF4, "LatentHF", oldLatentHFValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newLatentHFNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_Rn24h.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New Rn24h file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_Rn24h.nc",sep="")
	oldRn24hValues<-ncvar_get(nc,fic)
	newRn24hValues<-ncvar_def("Rn24h","daily", list(dimLonDef, dimLatDef, tdim), longname="Rn24h", missval=NaN, prec="double")
	nc_close(nc)
	newRn24hNCDF4<-nc_create(file_output,newRn24hValues)
	ncvar_put(newRn24hNCDF4, "Rn24h", oldRn24hValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newRn24hNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_LatentHF24h.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New LatentHF24h file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_LatentHF24h.nc",sep="")
	oldLatentHF24hValues<-ncvar_get(nc,fic)
	newLatentHF24hValues<-ncvar_def("LatentHF24h","daily", list(dimLonDef, dimLatDef, tdim), longname="LatentHF24h", missval=NaN, prec="double")
	nc_close(nc)
	newLatentHF24hNCDF4<-nc_create(file_output,newLatentHF24hValues)
	ncvar_put(newLatentHF24hNCDF4, "LatentHF24h", oldLatentHF24hValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newLatentHF24hNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_L.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New L file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_L.nc",sep="")
	oldLValues<-ncvar_get(nc,fic)
	newLValues<-ncvar_def("L","daily", list(dimLonDef, dimLatDef, tdim), longname="L", missval=NaN, prec="double")
	nc_close(nc)
	newLNCDF4<-nc_create(file_output,newLValues)
	ncvar_put(newLNCDF4, "L", oldLValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newLNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_y01.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New L file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_y01.nc",sep="")
	oldLValues<-ncvar_get(nc,fic)
	newLValues<-ncvar_def("y01","daily", list(dimLonDef, dimLatDef, tdim), longname="L", missval=NaN, prec="double")
	nc_close(nc)
	newLNCDF4<-nc_create(file_output,newLValues)
	ncvar_put(newLNCDF4, "y01", oldLValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newLNCDF4)
	
	#print(proc.time())
	
	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_y2.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New L file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_y2.nc",sep="")
	oldLValues<-ncvar_get(nc,fic)
	newLValues<-ncvar_def("y2","daily", list(dimLonDef, dimLatDef, tdim), longname="L", missval=NaN, prec="double")
	nc_close(nc)
	newLNCDF4<-nc_create(file_output,newLValues)
	ncvar_put(newLNCDF4, "y2", oldLValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newLNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_x200.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New L file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_x200.nc",sep="")
	oldLValues<-ncvar_get(nc,fic)
	newLValues<-ncvar_def("x200","daily", list(dimLonDef, dimLatDef, tdim), longname="L", missval=NaN, prec="double")
	nc_close(nc)
	newLNCDF4<-nc_create(file_output,newLValues)
	ncvar_put(newLNCDF4, "x200", oldLValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newLNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_psi01.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New L file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_psi01.nc",sep="")
	oldLValues<-ncvar_get(nc,fic)
	newLValues<-ncvar_def("psi01","daily", list(dimLonDef, dimLatDef, tdim), longname="L", missval=NaN, prec="double")
	nc_close(nc)
	newLNCDF4<-nc_create(file_output,newLValues)
	ncvar_put(newLNCDF4, "psi01", oldLValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newLNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_psi2.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New L file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_psi2.nc",sep="")
	oldLValues<-ncvar_get(nc,fic)
	newLValues<-ncvar_def("psi2","daily", list(dimLonDef, dimLatDef, tdim), longname="L", missval=NaN, prec="double")
	nc_close(nc)
	newLNCDF4<-nc_create(file_output,newLValues)
	ncvar_put(newLNCDF4, "psi2", oldLValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newLNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
	var_output<-paste(dados$Path.Output[1],"/",fic,"_psi200.nc",sep="")
	nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New L file name
	file_output<-paste(dados$Path.Output[1],"/",fic,"_psi200.nc",sep="")
	oldLValues<-ncvar_get(nc,fic)
	newLValues<-ncvar_def("psi200","daily", list(dimLonDef, dimLatDef, tdim), longname="L", missval=NaN, prec="double")
	nc_close(nc)
	newLNCDF4<-nc_create(file_output,newLValues)
	ncvar_put(newLNCDF4, "psi200", oldLValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
	nc_close(newLNCDF4)
	
	#print(proc.time())
}

tryCatch({
  res <- withTimeout({
    phase2();
  }, timeout=5400);
}, TimeoutException=function(ex) {
  cat("Image phase two processing timedout. Exiting with 124 code...\n");
  quit("no", 124, FALSE)
})

print(proc.time())
print("Fim da fase 3 - Calculo da evapotranspiracao")
