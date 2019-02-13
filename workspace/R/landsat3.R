library(raster)

radiance <- function(){
  print("Radiancia")
  rad <- NULL
  if(n.sensor < 8){
    rad <- list()
    if (n.sensor==5) r<-7 else r<-8
    for(i in 1:r){
      rad[[i]]<- image.rec[[i]][]*p.s$Grescale[i]+p.s$Brescale[i]
      rad_subset <- rad[[i]] < 0 & !is.na(rad[[i]])
      rad[[i]][rad_subset]<-0
    }
  } else{
    rad.mult<-as.numeric(MTL$V2[MTL$V1==grep(pattern ="RADIANCE_MULT_BAND_10 ", MTL$V1, value = TRUE)])
    rad.add<-as.numeric(MTL$V2[MTL$V1==grep(pattern ="RADIANCE_ADD_BAND_10 ", MTL$V1, value = TRUE)])
    rad <- image.rec[[7]][]*rad.mult+rad.add
  }
  return(rad)
}

reflectance <- function(rad){
  print("Refletancia")
  ref <- list()
  if(n.sensor < 8){
    if (n.sensor==5) r<-7 else r<-8
    for(i in 1:r){
      ref[[i]] <- pi*rad[[i]]*d_sun_earth$dist[Dia.juliano]^2/(p.s$ESUN[i]*costheta)
    }
  } else{
    for(i in 1:6){
      ref[[i]]<-(image.rec[[i]][]*0.00002-0.1)/costheta
    }
  }
  return(ref)
}

surfaceAlbedo <- function(ref){
  print("Albedo")
  print(n.sensor)

  r <- NULL
  if (n.sensor == 5) r <-7
  if (n.sensor == 7) r <-8
  if (n.sensor == 8) r <-6

  albedo <- ref[[1]]*p.s$wb[1] + ref[[2]]*p.s$wb[2] + ref[[3]]*p.s$wb[3] + ref[[4]]*p.s$wb[4] +
    ref[[5]]*p.s$wb[5] + ref[[r]]*p.s$wb[r]

  albedo <- (albedo-0.03)/tal[]^2
  return(albedo)
}

incomingShortwaveRadiation <- function(){
  print("Incoming short wave radiation")
  rs <-1367*costheta*tal[]/d_sun_earth$dist[Dia.juliano]^2
  return(rs)
}

ndvi <- function(ref){
  print("ndvi")
  NDVI <- (ref[[4]]-ref[[3]])/(ref[[4]]+ref[[3]])
  return(NDVI)
}

savi <- function(ref){
  print("savi")
  SAVI <- ((1+0.05)*(ref[[4]]-ref[[3]]))/(0.05+ref[[4]]+ref[[3]])
  return(SAVI)
}

lai <- function(SAVI){
  print("savi")
  LAI <- SAVI
  SAVI_subset1 <- SAVI > 0.687 & !is.na(SAVI)
  SAVI_subset2 <- SAVI <= 0.687 & !is.na(SAVI)
  SAVI_subset3 <- SAVI < 0.1 & !is.na(SAVI)
  LAI[SAVI_subset1] <- 6
  LAI[SAVI_subset2] <- -log((0.69 - SAVI[SAVI_subset2]) / 0.59) / 0.91
  LAI[SAVI_subset3] <- 0
}

evi <- function(ref){
  print("evi")
  EVI <- 2.5*((ref[[4]]-ref[[3]])/(ref[[4]]+(6*ref[[3]])-(7.5*ref[[1]])+1))
  return(EVI)
}

narrowBandSurfaceEmissivity <- function(LAI, NDVI){
  print("narrow band surface emissivity")
  enb <- 0.97 + 0.0033*LAI
  enb[NDVI < 0 | LAI > 2.99] <- 0.98
  return(enb)
}

broadBandSurfaceEmissivity <- function(LAI, NDVI){
  print("broad band surface emissivity")
  eo <- 0.95+0.01*LAI
  eo[NDVI < 0 | LAI > 2.99] <- 0.98
  return(eo)
}

surfaceTemperature <- function(enb, rad){
  print("temperatura da superficie")
  ts <- NULL
  print(n.sensor)
  if(n.sensor == 5){
    k1 <- 607.76
    k2 <- 1260.56
    ts <- k2/log((enb*k1/rad[[6]])+1)
  }

  if(n.sensor == 7){
    k1 <- 666.09
    k2 <- 1282.71
    ts <- k2/log((enb*k1/rad[[7]])+1)
  }

  if(n.sensor == 8){
    k1 <- 774.8853
    k2 <- 1321.0789
    ts <- k2/log((enb*k1/rad)+1)
  }

  return(ts)
}

outgoingLongwaveRadiation <- function(eo, ts){
  print("outgoing long wave radiation")
  rlsup <- eo*5.67*10^-8*ts^4
  return(rlsup)
}

atmosphericEmissivity <- function(){
  print("atmospheric emissivity")
  ea <- 0.85*(-1*log(tal[]))^0.09
  return(ea)
}

atmosphericOutgoingLongwaveRadiation <- function(ea){
  print("atmospheric outgoing long wave radiation")
  RLatm<-ea*5.67*10^-8*(table.sw$V7[2]+273.15)^4
  return(RLatm)
}

netRadiationFlux <- function(rs, albedo, RLatm, RLsup, eo){
  print("net radiation flux")
  Rn <- rs - rs*albedo + RLatm - RLsup - (1-eo)*RLatm
  Rn[Rn < 0] <- 0
  return(Rn)
}

soilHeatFlux <- function(NDVI, TS, albedo, Rn){
  print("soil heat flux")
  G_temp_1 <-(NDVI>=0)*(((TS-273.15)*(0.0038+0.0074*albedo)*(1-0.98*NDVI^4))*Rn)
  G_temp_2 <- (NDVI<0)*(0.5*Rn)
  G <- G_temp_1 + G_temp_2
  G[G < 0] <- 0
  return(G)
}

landsat<-function(){
  
  rad <- radiance()
  print(rad)

  ref <- reflectance(rad)
  print(ref)

  albedo <- surfaceAlbedo(ref)
  print(albedo)

  rs <- incomingShortwaveRadiation()
  print(rs)

  NDVI <- ndvi(ref)
  print(NDVI)

  SAVI <- savi(ref)
  print(SAVI)

  LAI <- lai(SAVI)
  print(LAI)

  EVI <- evi(ref)
  print(EVI)

  Enb <- narrowBandSurfaceEmissivity(LAI, NDVI)
  print(Enb)

  Eo <- broadBandSurfaceEmissivity(LAI, NDVI)
  print(Eo)

  TS <- surfaceTemperature(Enb, rad)
  print(TS)

  RLsup <- outgoingLongwaveRadiation(Eo, TS)
  print(RLSup)

  Ea <- atmosphericEmissivity()
  print(Ea)

  RLatm <- atmosphericOutgoingLongwaveRadiation(Ea)
  print(RLatm)

  Rn <- netRadiationFlux(rs, albedo, RLatm, RLsup, Eo)
  print(Rn)

  G <- soilHeatFlux(NDVI, TS, albedo, Rn)
  print(G)

  print(class(Rn))
  print(class(TS))
  print(class(NDVI))
  print(class(EVI))
  print(class(LAI))
  print(class(G))
  print(class(albedo))
  print(Rn)
  print(G)


  print("passou em todas")
  output <- stack(Rn,TS,NDVI,EVI,LAI,G,albedo)
  print("empilhado")
  return(output)
}
