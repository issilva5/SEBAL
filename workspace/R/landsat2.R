trunc <- function(obj, digits) {
    pot <- 10 ** digits
    obj <- as.integer(obj[] * pot)
    obj <- obj[] / pot
    return(obj)
}

landsat<-function(){
  if (n.sensor<8){
   # print("Radiancia")
    rad<-list()
    if (n.sensor==5) r<-7 else r<-8
    for(i in 1:r){
      rad[[i]]<- image.rec[[i]][]*p.s$Grescale[i]+p.s$Brescale[i]
      rad_subset <- rad[[i]] < 0 & !is.na(rad[[i]])
      rad[[i]][rad_subset]<-0
    }
  
  print("Radiancia")
  print(proc.time())
    
  print("Reflectancia")
  ref<-list()
  for(i in 1:r){
    ref[[i]]<-pi*rad[[i]]*d_sun_earth$dist[Dia.juliano]^2/
      (p.s$ESUN[i]*costheta)
  }
  print(proc.time())

  print("Albedo de superfcie")
  alb_temp<-ref[[1]]*p.s$wb[1]+ref[[2]]*p.s$wb[2]+ref[[3]]*p.s$wb[3]+ref[[4]]*p.s$wb[4]+
    ref[[5]]*p.s$wb[5]+ref[[r]]*p.s$wb[r]
  alb_temp<-(alb_temp-0.03)/tal[]^2
  print(proc.time())
  print(alb_temp)
  print("Radiao de onda curta incidente")
  Rs_temp<-1367*costheta*tal[]/d_sun_earth$dist[Dia.juliano]^2 # C�u claro
  print(proc.time())
  
  print("NDVI")
  #print("NDVI,SAVI,LAI e EVI")
  NDVI_temp<-(ref[[4]]-ref[[3]])/(ref[[4]]+ref[[3]])
  print(proc.time())
  print("SAVI")
  print(proc.time())
  SAVI_temp<- ((1+0.05)*(ref[[4]]-ref[[3]]))/(0.05+ref[[4]]+ref[[3]])
  print("LAI")
  LAI_temp <- SAVI_temp
  SAVI_subset1 <- SAVI_temp > 0.687 & !is.na(SAVI_temp)
  SAVI_subset2 <- SAVI_temp <= 0.687 & !is.na(SAVI_temp)
  SAVI_subset3 <- SAVI_temp < 0.1 & !is.na(SAVI_temp)
  LAI_temp[SAVI_subset1] <- 6
  LAI_temp[SAVI_subset2] <- -log(
    (0.69 - SAVI_temp[SAVI_subset2]) / 0.59
  ) / 0.91
  LAI_temp[SAVI_subset3] <- 0
  print(proc.time())
  print("EVI")
  EVI_temp<-2.5*((ref[[4]]-ref[[3]])/(ref[[4]]+(6*ref[[3]])-(7.5*ref[[1]])+1))
  print(proc.time())
  
  print("Emissividade Enb")
  Enb_temp <- 0.97+0.0033*LAI_temp
  Enb_temp[NDVI_temp<0 | LAI_temp>2.99]<-0.98
  print(proc.time())  

  print("Emissividade Eo")
  Eo_temp <- 0.95+0.01*LAI_temp
  Eo_temp[NDVI_temp<0 | LAI_temp>2.99]<-0.98
  
  print("Temperatura de Supercie em Kelvin (TS)")
  if (n.sensor==5) k1<-607.76 else k1<-666.09   #Constante Temperatura de superf�cie
  if (n.sensor==7) k2<-1260.56 else k2<-1282.71 #Constante Temperatura de superf�cie
  if (n.sensor==5) TS_temp<-k2/log((Enb_temp*k1/rad[[6]])+1) else TS_temp<-k2/log((Enb_temp*k1/rad[[7]])+1)
  print(proc.time())  

  print("Radicao de onda longa emitida pela superfcie (RLsup)")
  RLsup_temp<-Eo_temp*5.67*10^-8*TS_temp^4
  print(proc.time())

  print("Emissividade atmosfrica (Ea)")
  Ea_temp<-0.85*(-1*log(tal[]))^0.09 # C�u Claro
  print(proc.time())

  print("Radiao de onda longa emitida pela atmosfera (RLatm)")
  RLatm_temp<-Ea_temp*5.67*10^-8*(table.sw$V7[2]+273.15)^4
  print(proc.time())
  
  print("Saldo de radiao Instantnea (Rn)")
  Rn_temp<- Rs_temp-Rs_temp*alb_temp+RLatm_temp-RLsup_temp-(1-Eo_temp)*RLatm_temp
  Rn_temp[Rn_temp<0]<-0
  print(proc.time())

  print("Fluxo de Calor no Solo (G)")
  G_temp_1<-(NDVI_temp>=0)*(((TS_temp-273.15)*(0.0038+0.0074*alb_temp)*(1-0.98*NDVI_temp^4))*Rn_temp)
  G_temp_2 <- (NDVI_temp<0)*(0.5*Rn_temp)
  G_temp <- G_temp_1 + G_temp_2
  G_temp[G_temp<0]<-0
  print(proc.time())

  alb <- tal
  alb[] <- alb_temp
  Rn<-tal
  Rn[]<-Rn_temp
  TS<-tal 
  TS[]<-TS_temp
  NDVI<-tal 
  NDVI[]<-NDVI_temp
  EVI<-tal
  EVI[]<-EVI_temp
  LAI<-tal
  LAI[]<-LAI_temp
  G<-tal
  G[]<-G_temp
  
  output<-stack(Rn,TS,NDVI,EVI,LAI,G,alb)
  print(proc.time())
  } else {
    #print("Radi�ncia")
    rad.mult<-as.numeric(MTL$V2[MTL$V1==grep(pattern ="RADIANCE_MULT_BAND_10 ",
                                             MTL$V1, value = TRUE)])
    rad.add<-as.numeric(MTL$V2[MTL$V1==grep(pattern ="RADIANCE_ADD_BAND_10 ",
                                            MTL$V1, value = TRUE)])
    rad10<- image.rec[[7]][]*rad.mult+rad.add
    rad10 <- trunc(rad10, 4)
    print(proc.time())   
    #print("Reflect�ncia")
    ref<-list()
    for(i in 1:6){
      ref[[i]]<-(image.rec[[i]][]*0.00002-0.1)/costheta
      ref[[i]] <- trunc(ref[[i]], 4)
    }
    print(proc.time())
    print("Albedo")
    alb_temp<-ref[[1]]*p.s$wb[1]+ref[[2]]*p.s$wb[2]+ref[[3]]*p.s$wb[3]+ref[[4]]*p.s$wb[4]+
      ref[[5]]*p.s$wb[5]+ref[[6]]*p.s$wb[6]
    alb_temp<-(alb_temp-0.03)/tal[]^2
    alb_temp <- trunc(alb_temp, 4)
    #print(alb_temp)
    print(proc.time())

    print("Rs")
    Rs_temp<-1367*costheta*tal[]/d_sun_earth$dist[Dia.juliano]^2 # C�u claro
    Rs_temp <- trunc(Rs_temp, 4)
    #print(Rs_temp)
    print(proc.time())

    print("NDVI,SAVI,LAI e EVI")
    NDVI_temp<-(ref[[4]]-ref[[3]])/(ref[[4]]+ref[[3]])
    NDVI_temp <- trunc(NDVI_temp, 4)
    #print(NDVI_temp)
    print(proc.time())

    #NDVI<-tal 
    #NDVI[]<-NDVI_temp

    #print("NCol:")
    #print(ncol(NDVI))
    #print("NRow:")
    #print(nrow(NDVI))

    #xy <- xyFromCell(NDVI, 3886100) #500 * 7771 + 600
    #extxy <- extract(NDVI, xy, buffer=105)
    #print(xy)
    #print("Col:")
    #print(colFromX(NDVI, xy[1, 1]))
    #print("Line:")
    #print(rowFromY(NDVI, xy[1, 2]))
    #print("Ponto 1 - 3886100")
    #print(extxy)

    #xy <- xyFromCell(NDVI, 4000000)
    #extxy <- extract(NDVI, xy, buffer=105)
    #print(xy)
    #print("Col:")
    #print(colFromX(NDVI, xy[1, 1]))
    #print("Line:")
    #print(rowFromY(NDVI, xy[1, 2]))
    #print("Ponto 2 - 4000000")
    #print(extxy)

    #xy <- xyFromCell(NDVI, 2584720)
    #extxy <- extract(NDVI, xy, buffer=105)
    #print(xy)
    #print("Col:")
    #print(colFromX(NDVI, xy[1, 1]))
    #print("Line:")
    #print(rowFromY(NDVI, xy[1, 2]))
    #print("Ponto 3 - 2584720")
    #print(extxy)

    #quit()

    SAVI_temp<- ((1+0.05)*(ref[[4]]-ref[[3]]))/(0.05+ref[[4]]+ref[[3]])
    SAVI_temp <- trunc(SAVI_temp, 4)
    #print(SAVI_temp)
    print(proc.time())
    LAI_temp <- SAVI_temp
    SAVI_subset1 <- SAVI_temp > 0.687 & !is.na(SAVI_temp)
    SAVI_subset2 <- SAVI_temp <= 0.687 & !is.na(SAVI_temp)
    SAVI_subset3 <- SAVI_temp < 0.1 & !is.na(SAVI_temp)
    LAI_temp[SAVI_subset1] <- 6
    LAI_temp[SAVI_subset2] <- -log(
      (0.69 - SAVI_temp[SAVI_subset2]) / 0.59
    ) / 0.91
    LAI_temp[SAVI_subset3] <- 0
    LAI_temp <- trunc(LAI_temp, 4)
    #print(LAI_temp)
    print(proc.time())
    
    EVI_temp<-2.5*((ref[[4]]-ref[[3]])/(ref[[4]]+(6*ref[[3]])-(7.5*ref[[1]])+1))
    EVI_temp <- trunc(EVI_temp, 4)
    #print(EVI_temp)
    print(proc.time())
    
    print("Emissividade Enb")
    Enb_temp <- 0.97+0.0033*LAI_temp
    Enb_temp[NDVI_temp<0 | LAI_temp>2.99]<-0.98
    Enb_temp <- trunc(Enb_temp, 4)
    #print(Enb_temp)
    print(proc.time())

    print("Emissividade Eo")
    Eo_temp <- 0.95+0.01*LAI_temp
    Eo_temp[NDVI_temp<0 | LAI_temp>2.99]<-0.98
    Eo_temp <- trunc(Eo_temp, 4)
    #print(Eo_temp)
    print(proc.time())
    
    print("TS")
    k1<-774.8853      #Constante Temperatura de superf�cie
    k2<-1321.0789     #Constante Temperatura de superf�cie
    TS_temp<-k2/log((Enb_temp*k1/rad10)+1)
    TS_temp <- trunc(TS_temp, 4)
    #print(TS_temp)
    print(proc.time())
    print("RLsup")
    RLsup_temp<-Eo_temp*5.67*10^-8*TS_temp^4
    RLsup_temp <- trunc(RLsup_temp, 4)
    #print(RLsup_temp)
    print(proc.time())
    print("Ea")
    Ea_temp<-0.85*(-1*log(tal[]))^0.09 # C�u Claro
    Ea_temp <- trunc(Ea_temp, 4)
    #print(Ea_temp)
    print(proc.time())
    
    print("RLatm")
    RLatm_temp<-Ea_temp*5.67*10^-8*(table.sw$V7[2]+273.15)^4
    RLatm_temp <- trunc(RLatm_temp, 4)
    #print(RLatm_temp)
    print(proc.time())
    print("Rn")
    Rn_temp<- Rs_temp-Rs_temp*alb_temp+RLatm_temp-RLsup_temp-(1-Eo_temp)*RLatm_temp
    Rn_temp[Rn_temp<0]<-0
    Rn_temp <- trunc(Rn_temp, 4)
    #print(Rn_temp)
    print(proc.time())
    print("Fluxo de Calor no Solo (G)")
    G_temp_1<-(NDVI_temp>=0)*(((TS_temp-273.15)*(0.0038+0.0074*alb_temp)*(1-0.98*NDVI_temp^4))*Rn_temp)
    G_temp_2 <- (NDVI_temp<0)*(0.5*Rn_temp)
    G_temp <- G_temp_1 + G_temp_2
    G_temp[G_temp<0]<-0
    G_temp <- trunc(G_temp, 4)
    #print(G_temp)
    print(proc.time())
    alb <- tal
    alb[] <- trunc(alb_temp, 4)
    Rn<-tal
    Rn[]<-trunc(Rn_temp, 4)
    TS<-tal 
    TS[]<-trunc(TS_temp, 4)
    NDVI<-tal 
    NDVI[]<-trunc(NDVI_temp, 4)
    EVI<-tal
    EVI[]<-trunc(EVI_temp, 4)
    LAI<-tal
    LAI[]<-trunc(LAI_temp, 4)
    G<-tal
    G[]<-trunc(G_temp, 4)
    
    output<-stack(Rn,TS,NDVI,EVI,LAI,G,alb)
    #print("empilhado")
    print(proc.time())
    } 
  return(output)
}
