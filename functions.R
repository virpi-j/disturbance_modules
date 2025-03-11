
## ---------------------------------------------------------------------
## FUNCTIONS
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
## MAIN SCRIPT: uncRun for random segments, uncSeg for random values for segments
## ---------------------------------------------------------------------
runModelAdapt <- function(deltaID,sampleID=1, climScen=0, outType="dTabs",rcps = "CurrClim",
                     harvScen,harvInten,easyInit=FALSE, CO2fixed=0,
                     forceSaveInitSoil=F, cons10run = F,
                     procDrPeat=F,coeffPeat1=-240,coeffPeat2=70,
                     coefCH4 = 0.34,#g m-2 y-1
                     coefN20_1 = 0.23,coefN20_2 = 0.077,#g m-2 y-1
                     landClassUnman=NULL,compHarvX = 0,P0currclim=NA, fT0=NA,TminTmax = NA,
                     toRaster=F, disturbanceON = NA, ingrowth = F, clcut = 1){
  # outType determines the type of output:
  # dTabs -> standard run, mod outputs saved as data.tables 
  # testRun-> test run reports the mod out and initPrebas as objects
  # ststDeadW -> initialize the dead Wood volume;
  # uncRun -> reports the output table for the regional uncertainty run
  # uncSeg -> reports the list of output table for the segment uncertainty run
  # cons10run -> flag for conservation areas 10% run
  print(paste("Prebas version", vPREBAS))
  a0 <- Sys.time()
  if(outType=="testRun_region"){
    outType <- "testRun"
    sampleID <- deltaID
  }
  set.seed(1)  
  # print(date())
  nSegs <- dim(ops[[sampleID]])[1]
  path_to_inputs <- "/scratch/project_2000994/PREBASruns/finRuns/"
  print(paste("start sample ID",sampleID,"dim =",dim(ops[[sampleID]])[1]))
  if(rcps!="CurrClim" & climScen<0)  print(paste("start delta ID",deltaID,": deltaT=", deltaTP[1,deltaID]," deltaP=", deltaTP[2,deltaID]))
  
  initilizeSoil=T ###flag for soil initialization 
  procInSample=F
  ####in the protection scenarios consider buffer to protection areas
  ####if cons10run == TRUE run the model considering 10% area is conservation area according to zonation results
  if(harvScen %in% c("protect","protectNoAdH","protectTapio") & cons10run==FALSE ){
    # sampleX$cons[sampleX$Wbuffer==1] <- 1
    load(paste0(path_to_inputs,"input/maakunta/maakunta_",r_no,"_IDsBuffer.rdata"))
    xDat <- buffDat
    procInSample = T
    initilizeSoil = F
  }
  if(cons10run){
    load(paste0("path_to_inputs,input/maakunta/maakunta_",r_no,"_IDsCons10.rdata"))
    xDat <- cons10Dat
    procInSample = T
    initilizeSoil = F
  }
  if(procInSample){  
    if(identical(landClassX,1:3)) load(paste0("../initSoilC/station",station_id,"_LandClass1to3.rdata"))
    if(identical(landClassX,1:2)) load(paste0("../initSoilC/station",station_id,"_LandClass1to2.rdata"))
    if(identical(landClassX,1)) load(paste0("../initSoilC/station",station_id,"_LandClass1.rdata"))
    setnames(xDat,"nPix","N")
    xDat[,area:=N*16^2/10000]
    setkey(ops[[sampleID]],maakuntaID)
    setkey(xDat,maakuntaID)
    maakX <- ops[[sampleID]]$maakuntaID[which(ops[[sampleID]]$maakuntaID %in% xDat$maakuntaID)]
    posX <- which(ops[[sampleID]]$maakuntaID %in% xDat$maakuntaID)
    ops[[sampleID]][maakuntaID %in% maakX]$N <- xDat[maakuntaID %in% maakX]$N
    ops[[sampleID]][maakuntaID %in% maakX]$area <- xDat[maakuntaID %in% maakX]$area
    
    selX <- xDat[!maakuntaID %in% maakX &
                   oldMaakID %in% maakX]
    ops[[sampleID]][,oldMaakID:=maakuntaID]
    
    selX$newCons <- NULL
    selX$Wbuffer <- NULL
    ops[[sampleID]]$Wbuffer <- NULL
    
    sampleX <- rbind(ops[[sampleID]],selX)
    sampleX$segID <- sampleX$maakuntaID
    initSoilC <- abind(initSoilC,initSoilC[posX,,,],along=1)
    
    ###remove N==0 -> all seggment within the buffer
    x0 <- which(sampleX$N==0)    
    sampleX <- sampleX[-x0]
    initSoilC <- initSoilC[-x0,,,]
    # data.all <- rbind(data.all[!maakuntaID %in% xDat$maakuntaID],xDat)
  }else{
    sampleX <- ops[[sampleID]]
  }
  if(nrow(sampleX)<9000 & climScen < 0) station_id <- paste0(station_id,"s")
  if(outType %in% c("uncRun","uncSeg")){
    if(!exists("area_tot")) area_tot <- sum(data.all$area) # ha
    sampleX[,area := 16^2/10000] 
    cA <- 1/nrow(sampleX) #area_tot/nrow(sampleX) 
    harvestLims <- as.numeric(harvestLimsr[sampleID,])
    HarvLimMaak[,1]<-harvestLims[1]*HarvLimMaak[,1]
    HarvLimMaak[,2]<-harvestLims[2]*HarvLimMaak[,2]
    if(outType=="uncRun"){
      coeffPeat1 <- EC1[sampleID]
      coeffPeat2 <- EC2[sampleID]
    }
    if(uncRCP>0) {rcps <- paste0(climMod[climModids[sampleID]],rcpx[uncRCP])}
    else {rcps <- "CurrClim"}
    print(paste0("Climate model deltaID=",deltaID,": ",rcps))
    #print(paste("sampleID",sampleID,"harvestLims ="))
    #print(HarvLimMaak[1,] * sum(sampleX$area)/sum(data.all$area))
  } else if(!"area"%in%colnames(sampleX)){
    sampleX[,area := N*16^2/10000] 
  }
  #print(sum(sampleX$area))
  if(!exists("area_tot")) area_tot <- sum(data.all$area)
  sampleX[,id:=climID]
  HarvLimX <- harvestLims * sum(sampleX$area)/area_tot
  nSample = nrow(sampleX)#200#nrow(data.all)
  
  # leave unmaned land classes in landClassUnman
  if(!is.null(landClassUnman)) sampleX[landclass %in% landClassUnman]$cons=1
  
  ## ---------------------------------------------------------
  i = 0
  rcpfile = rcps
  print(paste("Clim:",rcpfile))
  if(rcpfile=="CurrClim"){
    load(paste(climatepath_orig, "CurrClim",".rdata", sep=""))
    #####process data considering only current climate###
    # dat <- dat[rday %in% 1:10958] #uncomment to select some years (10958 needs to be modified)
    maxRday <- max(dat$rday)
    #xday <- c(dat$rday,(dat$rday+maxRday),(dat$rday+maxRday*2))
    xday <- c(dat$rday,(dat$rday+maxRday),(dat$rday+maxRday*2),
              (dat$rday+maxRday*3),(dat$rday+maxRday*4),(dat$rday+maxRday*5))
    dat = rbind(dat,dat,dat,dat,dat,dat)
    #dat <- dat[rep(1:nrow(dat),4),]
    dat[,rday:=xday]
    rm(list = "xday")
  } else if(rcpfile=="CurrClim_fmi"){
      datname <- load(paste0(workdir, fmi_vars_PREBAS_file))#"fmi_vars_PREBAS.rdata"))
      assign("dat",get(datname))
      lookupname <- load(paste0(workdir, climID_lookup_file))#"climID_lookup.rdata"))
      assign("lookup",get(lookupname))
      rm(list=c("datname","lookupname"))
      gc()
      #####process data considering only current climate###
      dat[,rday := as.numeric(dat$time)-min(as.numeric(dat$time))+1]
      # dat <- dat[rday %in% 1:10958] #uncomment to select some years (10958 needs to be modified)
      maxRday <- max(dat$rday)
      colnames(dat)[which(colnames(dat)=="precip")] <- "Precip"
      colnames(dat)[which(colnames(dat)=="tair")] <- "TAir"
      colnames(dat)[which(colnames(dat)=="par")] <- "PAR"
      colnames(dat)[which(colnames(dat)=="vpd")] <- "VPD"
      colnames(dat)[which(colnames(dat)=="co2")] <- "CO2"
      dat$CO2[which(is.na(dat$CO2))] <- max(na.omit(dat$CO2))
      
      # TminTmax array, repeat Tmin and Tmax for all climIDs
      print("Setup Tmin Tmax values...")
      
      tmp <-  t( dcast(dat[, list(id, rday, tmin)], rday ~ id,
                       value.var="tmin")[, -1])
      TminTmax <- array(0,c(length(unique(dat$id)),ncol(tmp),2))
      TminTmax[,,1] <- t( dcast(dat[, list(id, rday, tmin)], rday ~ id,
                         value.var="tmin")[, -1])
      TminTmax[,,2] <- t( dcast(dat[, list(id, rday, tmax)], rday ~ id,
                                   value.var="tmax")[, -1])
      print("done.")
    } else if(climScen>=0){
      climatepath = "/scratch/project_2000994/RCP/"
      load(paste(climatepath, rcpfile,".rdata", sep=""))
      missingIDs <- setdiff(unique(sampleX$id), unique(dat$id))
      if(length(missingIDs)>0){
        coords <- fread("/scratch/project_2000994/RCP/coordinates")
        for(i in 1:length(missingIDs)){
          idX <- order((coords$x - coords$x[missingIDs[i]])^2 + (coords$y - coords$y[missingIDs[i]])^2)
          idX <- idX[idX%in%unique(dat$id)][1]
          nn<-which(sampleX$id==missingIDs[i]) 
          sampleX[nn,climID:=idX]
          sampleX[nn,id:=idX]
          print(paste("SampleID",sampleID,"clim ids: ",missingIDs[i], "was replaced with",idX))
        }
      }
    } else {
      climatepath_adapt <- "/scratch/project_2000994/PREBASruns/adaptFirst/tempData/"
      #read.csv2(file=paste0(climatepath,rcpsFile),sep = ",")
      dat2 <- read.csv(paste0(climatepath_adapt, rcpfile)) 
      dat2 <- dat2[which(dat2$Year2>=startingYear & 
                           dat2$deltaT==deltaTP[1,deltaID] & 
                           dat2$Pchange==deltaTP[2,deltaID]),]
      climIDs <- unique(sampleX$climID)
      # TminTmax array, repeat Tmin and Tmax for all climIDs
      dat1 <- read.csv(paste0(climatepath_adapt, str_replace(rcpfile, "v1","tmin_and_tmax")))
      dat1 <- dat1[which(dat1$Year2>=startingYear & 
                           dat1$deltaT==deltaTP[1,deltaID] & 
                           dat1$Pchange==deltaTP[2,deltaID]),]
      TminTmax <- array(0,c(length(climIDs),dim(dat1)[1],2))
      TminTmax[,,1] <- t(array(dat1$Tmin_perturbed,c(dim(dat1)[1],length(climIDs))))
      TminTmax[,,2] <- t(array(dat1$Tmax_perturbed,c(dim(dat1)[1],length(climIDs))))
      print("Run with Tmin Tmax values.")
      # CO2 array
      CO2<-as.numeric(sub(",",".",CO2_RCPyears[match(dat2$Year2,CO2_RCPyears$year),(Co2Col+1)]))
      if(CO2fixed==0){
        dat2 <- data.table(id=sampleX$climID[1],rday=1:nrow(dat2),
                           #PAR=-0.894+1.8*dat2$GLOB,
                           PAR=1.8*dat2$GLOB/1000,
                           TAir=dat2$Tmean_constant,#detrended,
                           VPD=dat2$VPdef_constant,#detrended,
                           Precip=dat2$Pre_constant,
                           CO2=CO2)
      } else {
        dat2 <- data.table(id=sampleX$climID[1],
                           rday=1:nrow(dat2),
                           #PAR=-0.894+1.8*dat2$GLOB,
                           PAR=1.8*dat2$GLOB/1000,
                           TAir=dat2$Tmean_seasonal,#detrended,
                           VPD=dat2$VPdef_seasonal,#detrended,
                           Precip=dat2$Pre_seasonal,
                           CO2=CO2)
      }
      nr <- length(climIDs)
      clim <- list(PAR = t(replicate(nr,dat2$PAR)),
                   TAir = t(replicate(nr,dat2$TAir)),
                 VPD = t(replicate(nr,dat2$VPD)),
                 Precip = t(replicate(nr,dat2$Precip)),
                 CO2 = t(replicate(nr,dat2$CO2)),
                 id = climIDs)
    rownames(clim$PAR)<-climIDs     
    colnames(clim$PAR)<-1:ncol(clim$PAR)    
    rm(list="dat2")
  }
  
  gc()
  ## Prepare the same initial state for all harvest scenarios that are simulated in a loop below
  data.sample = sample_data.f(sampleX, nSample)
  if(rcpfile%in%c("CurrClim")) data.sample$id <- data.sample$CurrClimID
  if(rcpfile%in%c("CurrClim_fmi")) data.sample$id <- data.sample$climID <- data.sample$CurrClimID <- lookup$climID[((sampleID-1)*nSample+1):(sampleID*nSample)] # dat-tiedostossa uudet clim-coordinaatit
  
  areas <- data.sample$area
  totAreaSample <- sum(data.sample$area)
  #print(sum(areas))
  
  if(rcpfile%in%c("CurrClim","CurrClim_fmi") | climScen >0) clim <-prep.climate.f(dat, data.sample, startingYear, nYears, rcps = rcpfile)
  rm("dat")
  gc()
  
  Region = nfiareas[ID==r_no, Region]
  if(!is.na(disturbanceON[1]) & !rcpfile%in%c("CurrClim","CurrClim_fmi")) ingrowth <- T # in disturbances ingrowth is always T
  ## Second, continue now starting from soil SS
  print(paste("Ingrowth = ",ingrowth))
  print(disturbanceON)
  print("Start initializing...")
  initPrebas = create_prebas_input_adapt.f(r_no, clim, data.sample, nYears = nYears,harv=harvScen,
                                           startingYear = startingYear,domSPrun=domSPrun,
                                            HcFactorX=HcFactor, 
                                           climScen=climScen, sampleX=sampleX, #clcut = clcut, 
                                           P0currclim=P0currclim, fT0=fT0, ingrowth = ingrowth,
                                           TminTmax = TminTmax, disturbanceON = disturbanceON)
  print("...done.")
  rm(list=c("data.sample"))
  gc()
  
  opsna <- which(is.na(initPrebas$multiInitVar))
  initPrebas$multiInitVar[opsna] <- 0.
  
  initPrebas$pPRELES <- pPREL
  initPrebas$pCROBAS <- pCROB
  
  ##### if the mortality model flag is 13 uses 
  ##### mortMod=1 (reineke) for managed forests
  ##### mortMod=3 (reineke + empirical model) for unmanaged forests
  if(mortMod==13){
    initPrebas$mortMod <- c(1,3)#rep(1,dim(initPrebas$multiOut)[1])
  }
  
  ##here mix years for weather inputs for Curr Climate
  if(rcpfile%in%c("CurrClim","CurrClim_fmi") & nYears>9){
    set.seed(10)
    resampleYear <- sample(1:nYears,nYears)
    resampleYear[(2015:2024)-2015+1] <- (2015:2024)-2015+1
    initPrebas$ETSy <- initPrebas$ETSy[,resampleYear]
    initPrebas$P0y <- initPrebas$P0y[,resampleYear,]
    initPrebas$weather <- initPrebas$weather[,resampleYear,,]
    initPrebas$weatherYasso <- initPrebas$weatherYasso[,resampleYear,]
  }
  
  ### for adapt and protect scenario Replanting schemes 
  ### do not replant pine in sitetypes 1 and 2
  ### do not replant spruce in sitetypes higher than 3
  ### ensure minimum 20% birch at replanting
  if(harvScen %in% c("adapt","protect","protectNoAdH","protectTapio",
                     "adaptNoAdH","adaptTapio")){
    sitesXs <- which(initPrebas$siteInfo[,3]>3)
    jj <- which(initPrebas$initCLcutRatio[sitesXs,2]>0.)
    initPrebas$initCLcutRatio[sitesXs[jj],2] <- 0.
    
    sitesXs <- which(initPrebas$siteInfo[,3]<3)
    jj <- which(initPrebas$initCLcutRatio[sitesXs,1]>0.)
    initPrebas$initCLcutRatio[sitesXs[jj],1] <- 0.
    
    xx <- 1/rowSums(initPrebas$initCLcutRatio)
    initPrebas$initCLcutRatio <- sweep(initPrebas$initCLcutRatio,MARGIN = 1,xx, `*`)
    jj <- which(is.na(rowSums(initPrebas$initCLcutRatio)))
    initPrebas$initCLcutRatio[jj,] <- 0.
    initPrebas$initCLcutRatio[jj,3] <- 1
    jj <- which(initPrebas$initCLcutRatio[,3]<0.2)
    xx <- 1-(0.2 - initPrebas$initCLcutRatio[jj,3])
    initPrebas$initCLcutRatio[jj,1:2] <- 
      sweep(initPrebas$initCLcutRatio[jj,1:2],MARGIN=1,xx, `*`)
    initPrebas$initCLcutRatio[jj,3] <- 0.2
  }
  
  
  #i = i + 1

  ## Assign harvesting quota for the region based on volume (in NFI startingYear) and MELA
  if(regSets!="maakunta"){
    Region = nfiareas[ID==r_no, Region]
    if(harvScen=="NoHarv"){
      if(!exists("clcut") & !is.na(disturbanceON[1])) clcut <- -1
      initPrebas$ClCut = initPrebas$defaultThin = rep(clcut,nSample)
      HarvLim1 = 0
      harvInten = "NoHarv"
    }else if(harvScen=="Tapio"){
      HarvLim1 = 0
    }else{
      HarvLim0 = nfiareas[ID==r_no, VOL_fraction]*rem[Scenario == harvScen & Area == Region, "1990-2013"]
      HarvLim0  = (totAreaSample/1000) / nfiareas[ID == r_no, AREA] * 1e3 *HarvLim0
      HarvLim = nfiareas[ID==r_no, VOL_fraction]*rem[Scenario == harvScen & Area == Region, "2015-2024"]
      HarvLim  = (totAreaSample/1000) / nfiareas[ID == r_no, AREA] * 1e3 *HarvLim
      HarvLim1 <- rep(as.numeric(HarvLim),10)
      HarvLim = nfiareas[ID==r_no, VOL_fraction]*rem[Scenario == harvScen & Area == Region, "2025-2034"]
      HarvLim  = (totAreaSample/1000) / nfiareas[ID == r_no, AREA] * 1e3 *HarvLim
      HarvLim1 <- c(HarvLim1,rep(as.numeric(HarvLim),10))
      HarvLim = nfiareas[ID==r_no, VOL_fraction]*rem[Scenario == harvScen & Area == Region, "2035-2044"]
      HarvLim  = (totAreaSample/1000) / nfiareas[ID == r_no, AREA] * 1e3 *HarvLim
      HarvLim1 <- c(HarvLim1,rep(as.numeric(HarvLim),10))
      HarvLim = nfiareas[ID==r_no, VOL_fraction]*rem[Scenario == harvScen & Area == Region, "2045-2054"]
      HarvLim  = (totAreaSample/1000) / nfiareas[ID == r_no, AREA] * 1e3 *HarvLim
      HarvLim1 <- c(HarvLim1,rep(as.numeric(HarvLim),10))
      HarvLim = nfiareas[ID==r_no, VOL_fraction]*rem[Scenario == harvScen & Area == Region, "2055-2064"]
      HarvLim  = (totAreaSample/1000) / nfiareas[ID == r_no, AREA] * 1e3 *HarvLim
      HarvLim1 <- c(HarvLim1,rep(as.numeric(HarvLim),44))
    }
    ## In the model, harvests are always per hectar units. If 1000 pixels (nSample)
    ## are simulated it corresponds to 1000 hectars, although pixels are only 16x16 m2.
    ## Therefore, we need to apply the areal fraction of removals scenarios
    ## nfiareas are in 1000 ha, model takes Harvlim in m3, while removals from Mela are 1000 m3
    #      HarvLim  = (nSample/1000) / nfiareas[ID == r_no, AREA] * 1e3 *HarvLim
    if(year1harv==1){
      HarvLim1 <- HarvLimX
      if(harvInten == "Low"){ HarvLim1 <- HarvLimX * 0.6}
      if(harvInten == "MaxSust"){HarvLim1 <- HarvLimX * 1.2}
      if(harvScen == "NoHarv"){
        HarvLim1 <- HarvLimX * 0.
        initPrebas$ClCut = initPrebas$defaultThin = rep(clcut,nSample)
        harvInten = harvScen
      }
    }else{
      roundWood <- HarvLim1 * roundTotWoodRatio
      enWood <- HarvLim1 - roundWood
      HarvLim1 <- cbind(roundWood,enWood)
    }
  }else{
    HarvLim1 <- HarvLimMaak*1000*sum(areas)/area_tot
    
    # If simulation time period is longer than harvest limit data, add lines
    if(nrow(HarvLimMaak)<nYears) HarvLim1<-rbind(HarvLim1,HarvLim1[rep(nrow(HarvLim1),nYears-nrow(HarvLim1)),])
  
    if(harvInten == "Low"){ HarvLim1 <- HarvLim1 * 0.6}
    if(harvInten == "MaxSust"){HarvLim1 <- HarvLim1 * 1.2}
    if(harvScen == "NoHarv"){
      HarvLim1 <- HarvLim1 * 0.
      initPrebas$ClCut = initPrebas$defaultThin = rep(0,nSample)
      harvInten = harvScen
    }
  }          
  
  ###calculate clearcutting area for the sample
  #if(!is.na(cutArX)){
  print("calculating clearcutting areas")
  clcutArX <- clcutAr * sum(areas)/area_tot
  if(length(clcutArX)<nYears) clcutArX<-c(clcutArX,clcutArX[rep(length(clcutArX),nYears-length(clcutArX))])
  clcutArX <- cbind(clcutArX[1:nYears],0.)
  tendX <- tendingAr * sum(areas)/area_tot
  if(length(tendX)<nYears) tendX<-c(tendX,tendX[rep(length(tendX),nYears-length(tendX))])
  tendX <- cbind(tendX[1:nYears],0.)
  fThinX <- firstThinAr * sum(areas)/area_tot
  if(length(fThinX)<nYears) fThinX<-c(fThinX,fThinX[rep(length(fThinX),nYears-length(fThinX))])
  fThinX <- cbind(fThinX[1:nYears],0.)
  cutArX <- cbind(clcutArX,tendX)
  cutArX <- cbind(cutArX,fThinX)
  if(harvInten == "Low"){ cutArX <- cutArX * 0.6}
  if(harvInten == "MaxSust"){cutArX <- cutArX * 1.2}
  if(harvScen == "NoHarv"){cutArX <- cutArX * 0.}
  
  ###run PREBAS
  if(initilizeSoil){
    if(!(harvScen =="Base" & harvInten == "Base" & rcpfile%in%c("CurrClim","CurrClim_fmi"))){
        if(!harvScen %in% c("protect","protectNoAdH","protectTapio")){
          pathToInitiSoilC <- "/scratch/project_2000994/PREBASruns/adaptFirst/initSoilC/"
          if(identical(landClassX,1:3)) load(paste0(pathToInitiSoilC,"station",station_id,"_LandClass1to3.rdata"))
          if(identical(landClassX,1:2)) load(paste0(pathToInitiSoilC,"station",station_id,"_LandClass1to2.rdata"))
          if(identical(landClassX,1)) load(paste0(pathToInitiSoilC,"station",station_id,"_LandClass1.rdata"))
        }
    }
  }
  initPrebas$yassoRun <- rep(1,initPrebas$nSites)
  if(exists("initSoilC")){ 
    nl_initPrebasSoilC <- dim(initPrebas$soilC[,1,,,])[4]
    nl_initSoilC <- dim(initSoilC)[4]
    #if(dim(initPrebas$soilC[,1,,,])[4]!=dim(initSoilC)[4]){
    #  initPrebas$soilC[,1,,,] <- initSoilC
    #} else {
      initPrebas$soilC[,1,,,1:nl_initSoilC] <- initSoilC
    #}
  }
  
  print(paste0("harvest scenario ", harvScen))
  print(paste0("harvest intensity ", harvInten))
  HarvLimX <- HarvLim1[1:nYears,]
  
  if(harvScen %in% c("adapt","adaptNoAdH","adaptTapio")){
    if(harvScen=="adaptNoAdH"){
      compHarvX=0.
    }
    ###set parameters to decrease rotation length of 25% (start)
    load(paste0(path_to_inputs,"input/",regSets,"/pClCut_adapt/ClCutplots_maak",r_no,".rdata"))
    ClcutX <- updatePclcut(initPrebas,pClCut)
    initPrebas$inDclct <- ClcutX$inDclct
    initPrebas$inAclct <- ClcutX$inAclct
    initPrebas$thinInt <- rep(thinIntX,initPrebas$nSites)
    ###set parameters to decrease rotation length of 25% (end)
    if(harvScen=="adaptTapio"){
      region <- regionPrebas(initPrebas,compHarv=compHarvX,
                             fertThin = fertThin,nYearsFert = nYearsFert)
    }else{
      region <- regionPrebas(initPrebas, HarvLim = as.numeric(HarvLimX),
                             cutAreas = cutArX,compHarv=compHarvX,
                             fertThin = fertThin,nYearsFert = nYearsFert)
    }
  }else if(harvScen %in% c("Mitigation","MitigationNoAdH","MitigationTapio")){
    if(harvScen=="MitigationNoAdH"){
      compHarvX=0.
    }
    HarvLimX[,2]=0.
    initPrebas$energyCut <- rep(0,length(initPrebas$energyCut))
    ###set parameters to increase rotation length of 25% (start)
    load(paste0("path_to_inputs,input/",regSets,"/pClCut_mitigation/ClCutplots_maak",r_no,".rdata"))
    ClcutX <- updatePclcut(initPrebas,pClCut)
    initPrebas$inDclct <- ClcutX$inDclct
    initPrebas$inAclct <- ClcutX$inAclct
    initPrebas$thinInt <- rep(thinIntX,initPrebas$nSites)
    ###set parameters to increase rotation length of 25% (end)
    if(harvScen=="MitigationTapio"){
      region <- regionPrebas(initPrebas,compHarv=compHarvX)
    }else{
      region <- regionPrebas(initPrebas, HarvLim = as.numeric(HarvLimX),
                             cutAreas =cutArX,compHarv=compHarvX,
                             ageHarvPrior = ageHarvPriorX)
    }
  }else if(harvScen %in% c("protect","protectNoAdH","protectTapio")){
    if(harvScen=="protectNoAdH"){
      compHarvX=0.
    }
    ####no energy cuts
    HarvLimX[,2]=0.
    initPrebas$energyCut <- rep(0,length(initPrebas$energyCut))
    
    ###set parameters to increase rotation length of 25% (start)
    load(paste0(path_to_inputs,"input/",regSets,"/pClCut_mitigation/ClCutplots_maak",r_no,".rdata"))
    ClcutX <- updatePclcut(initPrebas,pClCut)
    initPrebas$inDclct <- ClcutX$inDclct
    initPrebas$inAclct <- ClcutX$inAclct
    initPrebas$thinInt <- rep(thinIntX,initPrebas$nSites)
    ###set parameters to increase rotation length of 25% (end)
    
    if(harvScen=="protectTapio"){
      region <- regionPrebas(initPrebas,
                             compHarv=compHarvX,oldLayer = 1)
    }else{
      region <- regionPrebas(initPrebas, HarvLim = as.numeric(HarvLimX),
                             cutAreas =cutArX,compHarv=compHarvX,
                             ageHarvPrior = ageHarvPriorX,
                             oldLayer = 1)
    }
  }else{
    if(harvScen=="baseTapio"){
      #save(initPrebas,compHarvX,file="baseTapioInputs.rdata")
      #print("baseTapio input data saved.")
      region <- regionPrebas(initPrebas,compHarv=compHarvX)
    }else{
      if(climScen>10){ save(HarvLimX, cutArX,compHarvX,file =paste0("testDataRegion",restrictionSwitch,".rdata"))
        print("region run data saved")
        }
      ##Don't pass minDharvX if NA
      if (is.na(minDharvX)) {
        region <- regionPrebas(initPrebas, HarvLim = as.numeric(HarvLimX),
                               cutAreas =cutArX,compHarv=compHarvX)
      } else {
        #print("save regionPrebas input")
        #save(initPrebas, HarvLimX, minDharvX,cutArX,compHarvX,file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/testRun",rcps,".rdata"))
        print("start regionPrebas...")
        region <- regionPrebas(initPrebas, HarvLim = as.numeric(HarvLimX),
                               minDharv = minDharvX,cutAreas =cutArX,
                               compHarv=compHarvX)
        print("... done.")
        
      }
      
    }
  }
  print(paste("runModel",deltaID,"completed")) ##################################################
  ###############################################################################################
  
  ##calculate steady state carbon from prebas litter 
  if(harvScen=="Base" & harvInten =="Base" & initilizeSoil & rcpfile%in%c("CurrClim","CurrClim_fmi")){
    initSoilC <- stXX_GV(region, 1)
    print(paste("initSoilC:",sampleID))
    #if(outType!="testRun" | forceSaveInitSoil)
      if(!outType %in% c("uncRun","uncSeg")){
        pathToInitiSoilC <- "/scratch/project_2000994/PREBASruns/adaptFirst/initSoilC/"
        if(identical(landClassX,1:3)) save(initSoilC,file=paste0(pathToInitiSoilC,"station",station_id,"_LandClass1to3.rdata"))
        if(identical(landClassX,1:2)) save(initSoilC,file=paste0(pathToInitiSoilC,"station",station_id,"_LandClass1to2.rdata"))
        if(identical(landClassX,1)) save(initSoilC,file=paste0(pathToInitiSoilC,"station",station_id,"_LandClass1.rdata"))
        print("initSoilC saved.")
      }
    #}
    ###run yasso (starting from steady state) using PREBAS litter
    # region <- yassoPREBASin(region,initSoilC)
    initPrebas$yassoRun <- rep(1,initPrebas$nSites)
    initPrebas$soilC[,1,,,] <- initSoilC
    if (is.na(minDharvX)) {
      region <- regionPrebas(initPrebas, HarvLim = as.numeric(HarvLimX),
                             cutAreas =cutArX,compHarv=compHarvX)
    } else {
      region <- regionPrebas(initPrebas, HarvLim = as.numeric(HarvLimX),
                             minDharv = minDharvX,cutAreas =cutArX,
                             compHarv=compHarvX)
    }
    # out <- region$multiOut[,,,,1]
  }
  print(paste("all runs done",deltaID))
  
  #####process drained Peat
  if(procDrPeat){
    siteDrPeat1 <- which(sampleX$pseudoptyp==400 & region$siteInfo[,3]<4)
    siteDrPeat2 <- which(sampleX$pseudoptyp==400 & region$siteInfo[,3]>=4)
    
    ###CH4 <- N20
    # converts coeef to ha
    coefCH4 = coefCH4/1000*10000 #g m-2 y-1 -> kg ha-1
    coefN20_1 = coefN20_1/1000*10000 #g m-2 y-1 -> kg ha-1
    coefN20_2 = coefN20_2/1000*10000 #g m-2 y-1 -> kg ha-1
    region$CH4emisDrPeat_kgyear[1:nSegs] <- 0
    region$CH4emisDrPeat_kgyear[siteDrPeat1] = coefCH4
    region$CH4emisDrPeat_kgyear[siteDrPeat2] = coefCH4
    #region$CH4emisDrPeat_kgyear = coefCH4*region$areas[siteDrPeat1] +
    #  coefCH4*region$areas[siteDrPeat2]
    region$N2OemisDrPeat_kgyear[1:nSegs] <- 0
    region$N2OemisDrPeat_kgyear[siteDrPeat1] = coefN20_1
    region$N2OemisDrPeat_kgyear[siteDrPeat2] = coefN20_2
    #region$N2OemisDrPeat_kgyear = coefN20_1*region$areas[siteDrPeat1] +
    #  coefN20_2*region$areas[siteDrPeat2]
    
    region$multiOut[siteDrPeat1,,46,,1] = 0.
    region$multiOut[siteDrPeat1,,46,,1] = region$multiOut[siteDrPeat1,,18,,1] - 
      region$multiOut[siteDrPeat1,,26,,1]/10 - region$multiOut[siteDrPeat1,,27,,1]/10 - 
      region$multiOut[siteDrPeat1,,28,,1]/10 - region$multiOut[siteDrPeat1,,29,,1]/10
    region$multiOut[siteDrPeat1,,46,1,1] = region$multiOut[siteDrPeat1,,46,1,1] + 
      coeffPeat1 +  region$GVout[siteDrPeat1,,5]
    
    region$multiOut[siteDrPeat2,,46,,1] = 0.
    region$multiOut[siteDrPeat2,,46,,1] = region$multiOut[siteDrPeat2,,18,,1] - 
      region$multiOut[siteDrPeat2,,26,,1]/10 - region$multiOut[siteDrPeat2,,27,,1]/10 - 
      region$multiOut[siteDrPeat2,,28,,1]/10 - region$multiOut[siteDrPeat2,,29,,1]/10
    region$multiOut[siteDrPeat2,,46,1,1] = region$multiOut[siteDrPeat2,,46,1,1] + 
      coeffPeat2 +  region$GVout[siteDrPeat2,,5]
    print("Drained peatlands processed.")
  } else {
    print("No peatland post-procession")
  }
  #####start initialize deadWood volume
  ## identify managed and unmanaged forests
  manFor <-  which(sampleX$cons==0)
  unmanFor <- which(sampleX$cons==1)
  if(outType=="ststDeadW" | (harvScen =="Base" & harvInten == "Base" & rcpfile%in%c("CurrClim","CurrClim_fmi"))){
    yearsDeadW <- 1:nYears
    manDeadW <- initDeadW(region,manFor,yearsDeadW)
    print(paste("dim manDeadW:",dim(manDeadW$ssDeadW)))
    print(paste("mean:",mean(manDeadW$deadWV)))
    if(length(unmanFor)>1){
      unmanDeadW <- initDeadW(region,unmanFor,yearsDeadW)
      print(paste("dim unmanDeadW:",dim(unmanDeadW$deadWV)))
      print(paste("mean:",mean(unmanDeadW$deadWV)))
    } else {
      unmanDeadW <- data.frame()
    }
    #if(exists("station_id")){ 
      pathToInitDeadW <- "/scratch/project_2000994/PREBASruns/adaptFirst/initDeadWVss/"
      save(unmanDeadW,manDeadW,file=paste0(pathToInitDeadW,"station",
                                         station_id,"_deadWV_mortMod",mortMod,".rdata"))
      print("deadWood volume at steady state saved")
    #}
  }else{
    pathToInitDeadW <- "/scratch/project_2000994/PREBASruns/adaptFirst/initDeadWVss/"
    load(paste0(pathToInitDeadW,"station",
                station_id,"_deadWV_mortMod",mortMod,".rdata"))
    DeadWInit <- matrix(0,nrow = nYears, ncol = dim(manDeadW$ssDeadW)[2])
    if(nYears == dim(manDeadW$ssDeadW)[1]){
      DeadWInit[1:nrow(manDeadW$ssDeadW),] <- manDeadW$ssDeadW
    } else {
      nnYears <- min(nYears,dim(manDeadW$ssDeadW)[1])
      DeadWInit[1:nnYears,] <- manDeadW$ssDeadW[1:nnYears,]
    }
    nl_deadWInit <- dim(DeadWInit)[2]
    nl_multiOut <- dim(region$multiOut)[4]
    region$multiOut[manFor,,8,1:min(nl_deadWInit,nl_multiOut),1] <- 
      region$multiOut[manFor,,8,1:min(nl_deadWInit,nl_multiOut),1] + 
      aperm(replicate(length(manFor),DeadWInit),c(3,1:2))
#    region$multiOut[manFor,,8,1:3,1] <- region$multiOut[manFor,,8,1:3,1] + 
#      aperm(replicate(length(manFor),(manDeadW$ssDeadW[1:nYears,])),c(3,1:2))
    if(length(unmanFor)>1){
      DeadWInit <- matrix(0,nrow = nYears, ncol = dim(unmanDeadW$ssDeadW)[2])
      nl_deadWInit <- dim(DeadWInit)[2]
      if(nYears == dim(manDeadW$ssDeadW)[1]){
        DeadWInit[1:nrow(unmanDeadW$ssDeadW),] <- unmanDeadW$ssDeadW
      } else {
        nnYears <- min(nYears,dim(unmanDeadW$ssDeadW)[1])
        DeadWInit[1:nnYears,] <- unmanDeadW$ssDeadW[1:nnYears,]
      }
      region$multiOut[unmanFor,,8,1:min(nl_deadWInit,nl_multiOut),1] <- 
        region$multiOut[unmanFor,,8,1:min(nl_deadWInit,nl_multiOut),1] + 
        aperm(replicate(length(unmanFor),DeadWInit),c(3,1:2))
#      region$multiOut[unmanFor,,8,1:3,1] <- region$multiOut[unmanFor,,8,1:3,1] + 
#        aperm(replicate(length(unmanFor),(unmanDeadW$ssDeadW[1:nYears,])),c(3,1:2))
    }
    print("deadWood volume update processed.")
  }
  ####end initialize deadWood Volume
  # Return outputs
  if(outType=="testRun") return(list(region = region,initPrebas=initPrebas, clim=clim, TminTmax=TminTmax)) 
  if(outType=="dTabs"){
    print("Calculate outputs...")
    output <- runModOutAdapt(sampleID,deltaID,sampleX,region,r_no,harvScen,harvInten,climScen, rcpfile,areas,
              colsOut1,colsOut2,colsOut3,varSel,sampleForPlots,toRaster=toRaster)
    print(output[c(1,12,13,15,26:nrow(output)),c(1,2,6,ncol(output))])
    print("all outs calculated")
    #print(output)
    print(paste("Time",Sys.time()-a0))
    
    return(output)
  } 
  if(outType=="uncRun"){
    # results for all pixels
    uncTab <- UncOutProc(varSel=varSel,#c(46,39,30,37), 
                         funX=funX,#rep("sum",4),
                         modOut=region,sampleID=sampleID,
                         finPeats=finPeats,sampleX=sampleX)#,
    yy <- which(sampleX$peatID == 100) # mineral soils
    uncTab <- cbind(uncTab,UncOutProc(varSel=varSel,#c(46,39,30,37), 
                                      funX=funX,#rep("sum",4),
                                      modOut=region,sampleID=sampleID,
                                      finPeats=finPeats,sampleX=sampleX,vname="min",
                                      evalSegs=yy))#,
    yy <- which(sampleX$peatID == 400) # drained peatlands
    uncTab <- cbind(uncTab,UncOutProc(varSel=varSel,#c(46,39,30,37), 
                                      funX=funX,#rep("sum",4),
                                      modOut=region,sampleID=sampleID,
                                      finPeats=finPeats,sampleX=sampleX,
                                      vname="drPeat",
                                      evalSegs=yy))#,
    yy <- which(sampleX$consArea == 1) # conservation areas
    uncTab <- cbind(uncTab,UncOutProc(varSel=varSel,#c(46,39,30,37), 
                                      funX=funX,#rep("sum",4),
                                      modOut=region,sampleID=sampleID,
                                      finPeats=finPeats,sampleX=sampleX,
                                      vname="cons",
                                      evalSegs=yy))#,
    yy <- which(sampleX$consArea == 0) # managed & poorly productive forest
    uncTab <- cbind(uncTab,UncOutProc(varSel=varSel,#c(46,39,30,37), 
                                      funX=funX,#rep("sum",4),
                                      modOut=region,sampleID=sampleID,
                                      finPeats=finPeats,sampleX=sampleX,
                                      vname="man_pprod",
                                      evalSegs=yy))#,
    return(uncTab)
  } 
  if(outType=="uncSeg"){
    uncSegTab <- UncOutProcSeg(varSel=varSel, funX=funX,
                               modOut=region,sampleX,colsOut1,colsOut2,colsOut3)
    return(uncSegTab)
  }
  # rm(list=c("region","initPrebas")); gc()
  # rm(list=setdiff(ls(), c(toMem,"toMem")))
  # rm(out); gc()
  # }###harvest loop
  # } ###region loop
  # }rcps loop
  print(paste("end deltaID",deltaID))
  rm(list=setdiff(ls(), c(toMem,"toMem"))); gc()
  
  #print(uncRun)
  # }
}

runModOutAdapt <- function(sampleID,deltaID,sampleX,modOut,r_no,harvScen,harvInten,climScen,rcpfile,areas,
                      colsOut1,colsOut2,colsOut3,varSel,sampleForPlots,toRaster){#},SBBbp,PI,pSBB){
  ####create pdf for test plots 
  marginX= 1:2#(length(dim(out$annual[,,varSel,]))-1)
  nas <- data.table()
  output <- data.frame()
  ij<-1
  for (ij in 1:length(varSel)) {
    # print(varSel[ij])
    if(funX[ij]=="baWmean"){
      outX <- data.table(segID=sampleX$segID,baWmean(modOut,varSel[ij]))
    }
    if(funX[ij]=="sum"){
      outX <- data.table(segID=sampleX$segID,apply(modOut$multiOut[,,varSel[ij],,1],marginX,sum))
    }
    ####test plot
    #print(outX)
    #if(sampleID==sampleForPlots){testPlot(outX,varNames[varSel[ij]],areas)}
    pX <- calculatePerCols(outX = outX)
#    pX <- calculatePerColsAllRows(outX = outX)
    vnam <- varOuts[ij]#varNames[varSel[ij]] 
    if(toRaster){
      if(vnam=="GPPTot/1000") vnam<-"GPPTot_1000"
      assign(vnam,pX)
      save(list=vnam,
           file=paste0(path_output,"weatherStation",station_id,"/",
                       vnam,
                       "_harscen",harvScen,
                       "_harInten",harvInten,"_",
                       rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
      rm(list=vnam); gc()
    }    
    ##check for NAs
    nax <- data.table(segID=unique(which(is.na(pX),arr.ind=T)[,1]))
    if(nrow(nax)>0){
      nax$var <- varNames[varSel[ij]]
      nax$sampleID <- sampleID
      nas <- rbind(nas,nax)
    } 
    pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
    pX <- c(var = varNames[varSel[ij]], pX)
    #pX[1] <- varNames[varSel[ij]]
    #names(pX)[1] <- "var"
    output <- rbind(output, pX)
    colnames(output) <- names(pX)
    #print(output)
  }
  # save NAs
  #  if(nrow(nas)>0){
  #    save(nas,file=paste0("NAs/NAs_forCent_",r_no,
  #                         "_","sampleID",sampleID,
  #                         "_harscen",harvScen,
  #                         "_harInten",harvInten,"_",
  #                         rcpfile,".rdata"))        
  #  }
  ####process and save special variales
  print(paste("start special vars",deltaID))
  output <- specialVarProcAdapt(sampleX,modOut,r_no,harvScen,harvInten,rcpfile,sampleID,
                 areas,sampleForPlots,output,toRaster=toRaster)#,SBBbp,PI,pSBB)
  return(output)
}



#sample_data.f = function(data.all, nSample) {
sample_data.f = function(sampleX, nSample) {
  cloudpixels = sampleX[, sum(ba==32766)]
  nonforest = sampleX[, sum(ba==32767)]
  forest = sampleX[, sum(ba< 32766)]
  AREA = (forest + cloudpixels) * 16 * 16 * 1000 #m2
  AREA_1000ha = AREA / 10000 / 1000
  
  ## REMOVE CLOUD COVERED, AND WHERE cons = NA (...? why)
  sampleX = sampleX[ba < 32766]
  sampleX = sampleX[!is.na(cons)]
  gc()
  ## REDUCE SAMPLE FOR TESTING ---------------------------------------
  smp = floor(seq(1, dim(sampleX)[1], len=nSample))
  data.sample = sampleX[smp,]
  
  # summary(data.sample[, 3:11])
  colnams <- c("regID",  "N",      "ba",     "age",    "dbh",    "pine",   "spruce", "birch" )
  #for (col in colnames(data.sample)[c(3, 5:11)]){ 
  for (col in colnames(data.sample)[match(colnams, colnames(sampleX))]){ 
    #print(col)
    set(data.sample, j=col, value=as.double(data.sample[[col]]))
  }  
  #if("y"%in%colnames(data.sample))set(data.sample, j="y", value=as.double(data.sample[[col]]))
  #if("x"%in%colnames(data.sample))set(data.sample, j="x", value=as.double(data.sample[[col]]))
  #if("lat"%in%colnames(data.sample))set(data.sample, j="lat", value=as.double(data.sample[[col]]))
  
    
  ## -----------------------------------------------------------------
  
  
  ## AVOID ZERO CASES
  
  data.sample$dbh = as.double(data.sample$dbh)
  
  data.sample[pine == 0 & spruce == 0 & decid ==0 & fert ==1, decid:=1  ]
  data.sample[pine == 0 & spruce == 0 & decid ==0 & fert <= 3 & fert > 1, spruce:=1  ]
  data.sample[pine == 0 & spruce == 0 & decid ==0 & fert >= 4, pine:=1  ]
  siteX <- union(which(data.sample$ba <=0.041),which(data.sample$h<= 15))
  siteX <- union(siteX,which(data.sample$dbh<=0.5))
  data.sample$nTree <- data.sample$ba/(pi/4*( data.sample$dbh/100)^2)
  siteNN <- which(data.sample$nTree>5000)
  siteX <- union(siteX,siteNN)
  data.sample[siteX,h:=15]
  data.sample[siteX,dbh:=0.5]
  data.sample[siteX,ba:=0.0431969]
  data.sample
}



# StartingYear = climate data that detrermines simulation period must have year greater than this.
create_prebas_input_adapt.f = function(r_no, clim, data.sample, nYears, harv,
                                 startingYear=0, domSPrun=0, clcut=1,
                                 HcFactorX=HcFactor, climScen=climScen, ingrowth=F,
                                 sampleX=sampleX, P0currclim=NA, fT0=NA, TminTmax=NA,
                                 disturbanceON = F) { # dat = climscendataset
  #domSPrun=0 initialize model for mixed forests according to data inputs 
  #domSPrun=1 initialize model only for dominant species 
  nSites <- nrow(data.sample)
  ###site Info matrix. nrow = nSites, cols: 1 = siteID; 2 = climID; 3=site type;
  ###4 = nLayers; 5 = nSpecies;
  ###6=SWinit;   7 = CWinit; 8 = SOGinit; 9 = Sinit
  
  siteInfo <- matrix(c(NA,NA,NA,160,0,0,20,3,3,413,0.45,0.118),nSites,12,byrow = T)
  #siteInfo <- matrix(c(NA,NA,NA,3,3,160,0,0,20),nSites,9,byrow = T)
  siteInfo[,1] <- data.sample$segID
  siteInfo[,2] <- as.numeric(data.sample[,id])
  siteInfo[,3] <- data.sample[,fert]
  
  ####### Wind disturbance module from Jonathan
  sid <- NA
  if("wind"%in%disturbanceON){# & !rcps%in%c("CurrClim","CurrClim_fmi")){
    ###
    # EXTRACT WIND SPEEDS 
    # for prebas wind disturbance module
    print("EXTRACT WIND SPEEDS")
    library(sf)
    library(mapview)
    library(stars)
    
    #load("/scratch/project_2000994/PREBASruns/adaptFirst/rasters/ops_Jyvaskyla.rdata") # load ops object / sampled segments
    #idcoord<-  ops[[1]][,c("segID", "x", "y")] # unlist, reduce to coordinates (likely centroids)
    idcoord <-  data.sample[,c("segID", "x", "y")] # unlist, reduce to coordinates (likely centroids)
    #names(idcoord) <- c("segID", "x", "y")
    
    segids <- st_as_sf(idcoord, coords = c("x","y")) # as sf (vector format)
    st_crs(segids) <- "EPSG:3067"
    
    #mapview(segids) # check visually 
    #length(unique(segids$segID)) #check number of unique ids/samples
    
    
    #### WIND SPEED ####
    
    # read wind speed dataset
    ws_full <- read_stars("/appl/data/geo/ilmatiede/wind_speed/Wind_10y_return_level.tif") # directly available in puhti!
    
    # extract wspeed at plot/centroid locations
    wspeed = st_extract(ws_full, segids) 
    
    segids$wspeed <- wspeed$Wind_10y_return_level.tif 
    
    #mapview(segids, zcol="wspeed")
    
    #st_write(segids, dsn="/scratch/project_2000994/PREBASruns/adaptFirst/rasters/segids_wspeed_jyvaskyla.gpkg", delete_dsn = TRUE)
    
    segids_dt<- data.frame(segids)
    
    segids_dt[2] <- NULL
    #write.csv(segids_dt, dsn="/scratch/project_2000994/PREBASruns/adaptFirst/rasters/segids_wspeed_jyvaskyla.csv")
    ###
    #Xy <- st_read("/scratch/project_2000994/PREBASruns/adaptFirst/rasters/segids_wdist.addvars_jyvaskyla.gpkg")
    #ni <- which(Xy$z %in%  data.sample$segID)
    ni <- which(segids_dt$segID %in%  data.sample$segID)
    ni <- ni[match(data.sample$segID,segids_dt$segID[ni])]
    
    sid <- matrix(0, nSites,10) #create input matrix
    # identical demo inputs for all sites
    sid[,1] <- segids_dt$wspeed[ni] # localised 10a return max wind speed (Ven??l??inen et al. 2017). Average 12.2. For test purposes, this can be set to e.g. 50 to trigger disturbances more frequently...
    sid[,2] <- sample(1:30, nSites, replace=T) # init for time since thinning (e.g. sampling 1:40)
    sid[,3] <- 1 #Xy$soiltype[ni] # soiltype (0 = mineral, coarse; 1 = mineral, fine; 2 = organic)
    sid[,4] <- 0 # shallowsoil (0 = F, >30cm, 1 = T, <30cm)
    
    # salvage logging/mgmt reaction parameters
    sid[,5] <- 10 # salvlog threshhold, m3/ha; if total damaged volume exceeds this, site is considered for salvage logging (removal of directly damaged timber)
    sid[,6] <- 1 # salvlog share, 0-1; share of sites over salvlog threshold where salvage logging is applied
    sid[,7] <- 1 # pharvtrees for salvage logging (share of directly damaged vol to be collected) !!NOTE: 0.1 still going to harvest residues after this, i.e. 'harvest as in regular thin' = 1!)
    sid[,8] <- 20 # mgmtreact threshold, m3/ha: if total damaged volume exceeds this, site is considered for prioritisation in randomised Tapio mgmt
    sid[,9] <- 1 # mgmtreact share, 0-1; share of sites over mgmtreact threshold where prioritisation is applied. Note: if mgmt react/prioritisation is applied, salvage logging is conducted as well.
    sid[,10] <- 1 # sevdistccshare: share of sites with reldamvol>0.5 or sevclass 3 where force clearcut is applied
  }
  ##
  # litterSize <- matrix(0,3,3)
  # litterSize[1,1:2] <- 30
  # litterSize[1,3] <- 10
  # litterSize[2,] <- 2
  
  ###Initialise model
  # initVardension nSites,variables, nLayers
  # variables: 1 = species; 2 = Age; 3 = H; 4=dbh; 5 = ba; 6 = Hc
  initVar <- array(NA, dim=c(nSites,7,3))
  data.sample[,baP:= (ba * pine/(pine+spruce+decid))]
  data.sample[,baSP:= (ba * spruce/(pine+spruce+decid))]
  data.sample[,baB:= (ba * decid/(pine+spruce+decid))]
  data.sample[,dbhP:= dbh]
  data.sample[,dbhSP:= dbh]
  data.sample[,h:= h/10]
  data.sample[,hP:= h]
  data.sample[,hSP:= h]
  
  data.sample[,N:=ba/(pi*(dbh/2)^2/10000)]
  
  areas <- data.sample$area
  
  initVar[,1,] <- as.numeric(rep(1:3,each=nSites))
  initVar[,2,] <- round(as.numeric(data.sample[,age]))
  initVar[,3,] <- as.numeric(data.sample[,h])
  # initVar[,3,][which(initVar[,3,]<1.5)] <- 1.5  ####if H < 1.5 set to 1.5
  initVar[,4,] <- as.numeric(data.sample[,dbh])
  
  if(domSPrun==1){
    ##initialize model only for dominant species##
    initVar[,5,] = 0.
    ix = unlist(data.sample[, which.max(c(pine, spruce, decid)), by=1:nrow(data.sample)] [, 2])
    for(jx in 1:nSites) initVar[jx,5,ix[jx]] = as.numeric(data.sample[, ba])[jx]
  } else{
    ###initialize model for mixed forest runs
    initVar[,5,1] <- as.numeric(data.sample[,(ba * pine/(pine+spruce+decid))])
    initVar[,5,2] <- as.numeric(data.sample[,(ba * spruce/(pine+spruce+decid))])
    initVar[,5,3] <- as.numeric(data.sample[,(ba * decid/(pine+spruce+decid))])
    
    if(TRUE){ #### if true will vary H and D of pine and spruce using siteType
      
      ###increase spruceP dbh 10% for spruceP sitetype 1:2
      minDelta <- 0.75
      data.sample[pine>0. & spruce >0. & fert<2.5,X:=pmax(minDelta,(ba-1.1*baSP-baB)/baP)]
      data.sample[pine>0. & spruce >0. & fert<2.5,dbhSP:=1.1*dbh]
      data.sample[pine>0. & spruce >0. & fert<2.5 & X==minDelta,dbhSP:=dbh*(ba-minDelta* baP-baB)/baSP]
      data.sample[pine>0. & spruce >0. & fert<2.5,dbhP:=X*dbh]
      data.sample[pine>0. & spruce >0. & fert<2.5 & dbhP<0.5,dbhSP:=pmax(0.5,((ba-(0.5/dbh)*baP-baB)/baSP))]
      data.sample[pine>0. & spruce >0. & fert<2.5 & dbhP<0.5,dbhP:=0.5]
      
      # data.sample[pine>0. & spruce >0. & fert<2.5 & baSP <= baP,dbhSP:=dbh * (ba - 0.9*baP - baB)/baSP]
      # data.sample[pine>0. & spruce >0. & fert<2.5 & baSP <= baP,dbhP:=pmax(0.9*dbh,0.3)]
      
      ####increase spruce h 10% for spruce sitetype 1:2
      data.sample[pine>0. & spruce >0. & fert<2.5, X:=pmax(minDelta,(ba-1.1*baSP-baB)/baP)]
      data.sample[pine>0. & spruce >0. & fert<2.5,hSP:=1.1*h]
      data.sample[pine>0. & spruce >0. & fert<2.5 & X==minDelta,hSP:=h*(ba-minDelta* baP-baB)/baSP]
      data.sample[pine>0. & spruce >0. & fert<2.5, hP:=X*h]
      data.sample[pine>0. & spruce >0. & fert<2.5 & hSP<1.5,hSP:=1.5]
      data.sample[pine>0. & spruce >0. & fert<2.5 & hP<1.5,hP:=1.5]
      
      # data.sample[pine>0. & spruce >0. & fert<2.5 & baSP <= baP,hSP:=h * (ba - 0.9*baP - baB)/baSP]
      # data.sample[pine>0. & spruce >0. & fert<2.5 & baSP <= baP,hP:=pmax(0.9*h,1.3)]
      #  
      ####increase spruce dbh 5% for spruce sitetype 3
      data.sample[pine>0. & spruce >0. & fert==3, X:=pmax(minDelta,(ba-1.05*baSP-baB)/baP)]
      data.sample[pine>0. & spruce >0. & fert==3, dbhP:=X*dbh]   
      data.sample[pine>0. & spruce >0. & fert==3, dbhSP:=1.05*dbh]
      data.sample[pine>0. & spruce >0. & fert==3 & X==minDelta,dbhSP:=dbh*(ba-minDelta* baP-baB)/baSP]
      data.sample[pine>0. & spruce >0. & fert==3 & dbhP<0.5,dbhSP:=pmax(1.5,((ba-(0.5/dbh)*baP-baB)/baSP)*dbh)]
      data.sample[pine>0. & spruce >0. & fert==3 & dbhP<0.5,dbhP:=0.5]
      
      # data.sample[pine>0. & spruce >0. & fert==3 & baSP <= baP,dbhSP:=pmin(25,(dbh * (ba - 0.95*baP - baB)/baSP))]
      # data.sample[pine>0. & spruce >0. & fert==3 & baSP <= baP,dbhP:=pmax(0.95*dbh,0.3)]
      
      ####increase spruce h 5% for spruce sitetype 3
      data.sample[pine>0. & spruce >0. & fert==3, X:=pmax(minDelta,(ba-1.05*baSP-baB)/baP)]
      data.sample[pine>0. & spruce >0. & fert==3, hP:=X*h]
      data.sample[pine>0. & spruce >0. & fert==3, hSP:=1.05*h]
      data.sample[pine>0. & spruce >0. & fert==3 & X==minDelta,hSP:=h*(ba-minDelta* baP-baB)/baSP]
      data.sample[pine>0. & spruce >0. & fert==3 & hSP<1.5, hSP:=1.5]
      data.sample[pine>0. & spruce >0. & fert==3 & hP<1.5, hP:=1.5]
      
      # data.sample[pine>0. & spruce >0. & fert==3 & baSP <= baP,hSP:=pmin(30.,(h * (ba - 0.95*baP - baB)/baSP))]
      # data.sample[pine>0. & spruce >0. & fert==3 & baSP <= baP,hP:=pmax(0.95*h,1.3)]
      
      ####increase pine dbh 10% for sitetype >= 4
      data.sample[pine>0. & spruce >0. & fert>3.5, X:=pmax(minDelta,(ba-1.1*baP-baB)/baSP)]
      data.sample[pine>0. & spruce >0. & fert>3.5, dbhSP:=X*dbh]
      data.sample[pine>0. & spruce >0. & fert>3.5, dbhP:=1.1*dbh]
      data.sample[pine>0. & spruce >0. & fert>3.5 & X==minDelta,dbhP:=dbh*(ba-minDelta*baSP-baB)/baP]
      data.sample[pine>0. & spruce >0. & fert>3.5 & dbhSP<0.5,dbhP:=pmax(1.5,((ba-(0.5/dbh)*baSP-baB)/baP)*dbh)]
      data.sample[pine>0. & spruce >0. & fert>3.5 & dbhSP<0.5,dbhSP:=0.5]
      # data.sample[pine>0. & spruce >0. & fert>3.5 & baP <= baSP,dbhP:=dbh * (ba - 0.9*baSP - baB)/baP]
      # data.sample[pine>0. & spruce >0. & fert>3.5 & baP <= baSP,dbhSP:=pmax(0.9*dbh,0.3)]
      ####increase pine h 10% for sitetype >= 4
      data.sample[pine>0. & spruce >0. & fert>3.5, X:=pmax(minDelta,(ba-1.1*baP-baB)/baSP)]
      data.sample[pine>0. & spruce >0. & fert>3.5,hSP:=X*h]
      data.sample[pine>0. & spruce >0. & fert>3.5,hP:=1.1*h]
      data.sample[pine>0. & spruce >0. & fert>3.5 & X==minDelta,hP:=h*(ba-minDelta*baSP-baB)/baP]
      data.sample[pine>0. & spruce >0. & fert>3.5 & hP<1.5,hP:=1.5]
      data.sample[pine>0. & spruce >0. & fert>3.5 & hSP<1.5,hSP:=1.5]
      # data.sample[pine>0. & spruce >0. & fert>3.5 & baP <= baSP,hP:=h * (ba - 0.9*baSP - baB)/baP]
      # data.sample[pine>0. & spruce >0. & fert>3.5 & baP <= baSP,hSP:=pmax(0.9*h,1.3)]
      
      initVar[,3,1] <- as.numeric(data.sample[,hP])
      initVar[,3,2] <- as.numeric(data.sample[,hSP])
      initVar[,4,1] <- as.numeric(data.sample[,dbhP])
      initVar[,4,2] <- as.numeric(data.sample[,dbhSP])
      
    }
    
  }
  
  # initVar[,6,] <- as.numeric(data.sample[,hc])
  
  if(harv %in% c("adapt","protect","protectNoAdH","protectTapio",
                 "adaptNoAdH","adaptTapio")){
    ####always the 3 species layers in this two scenarios
    ###check which BA ==0. and set to 0 the rest of the variable
    NoPine <- which(initVar[,5,1]==0.)
    NoSpruce <- which(initVar[,5,2]==0.)
    NoDecid <- which(initVar[,5,3]==0.)
    
    # siteInfo[NoPine,8] <- siteInfo[NoPine,8] - 1
    # siteInfo[NoSpruce,8] <- siteInfo[NoSpruce,8] - 1
    # siteInfo[NoDecid,8] <- siteInfo[NoDecid,8] - 1
    
    initVar[NoPine,3:6,1] <- 0.
    initVar[NoSpruce,3:6,2] <- 0.
    initVar[NoDecid,3:6,3] <- 0.
    # initVar[NoSpruce,,2] <- initVar[NoSpruce,,3]
    # initVar[NoPine,,1:2] <- initVar[NoPine,,2:3]
    
    # nLay1 <- which(siteInfo[,8]==1)
    # nLay2 <- which(siteInfo[,8]==2)
    # initVar[nLay1,c(1,3:6),2:3] <- 0
    # initVar[nLay2,c(1,3:6),3] <- 0
  }else{
    NoPine <- which(initVar[,5,1]==0.)
    NoSpruce <- which(initVar[,5,2]==0.)
    NoDecid <- which(initVar[,5,3]==0.)
    
    siteInfo[NoPine,8] <- siteInfo[NoPine,8] - 1
    siteInfo[NoSpruce,8] <- siteInfo[NoSpruce,8] - 1
    siteInfo[NoDecid,8] <- siteInfo[NoDecid,8] - 1
    
    initVar[NoPine,3:6,1] <- 0.
    initVar[NoSpruce,3:6,2] <- 0.
    initVar[NoDecid,3:6,3] <- 0.
    initVar[NoSpruce,,2] <- initVar[NoSpruce,,3]
    initVar[NoPine,,1:2] <- initVar[NoPine,,2:3]
    
    nLay1 <- which(siteInfo[,8]==1)
    nLay2 <- which(siteInfo[,8]==2)
    initVar[nLay1,3:6,2:3] <- 0
    initVar[nLay2,3:6,3] <- 0
  }
  
  siteInfo[, 2]  = match(as.numeric(siteInfo[, 2]), as.numeric(rownames(clim[[1]])))
  # siteInfo[, 2]  = match(siteInfo[,2], unique(dat$id))
  
  defaultThin=as.numeric(1-data.sample[, cons])
  energyCut <- ClCut <- as.numeric(1-data.sample[, cons])
  ## Set to match climate data years
  if(!exists("ftTapioParX")) ftTapioParX = ftTapio
  if(!exists("tTapioParX")) tTapioParX = tTapio
  #initVar[,6,] <- aaply(initVar,1,findHcNAs,pHcM)[,6,]*HcFactorX
  initVar[,6,] <- aaply(initVar,1,findHcNAs,pHcM,pCrobasX,HcModVx)[,6,]*HcFactorX
  set_thin_PROJ6_warnings(TRUE)
  xy <- sampleX[,c("segID","x","y")]
  coordinates(xy) <- c("x","y")
  proj4string(xy) <- crsX
  #cord = SpatialPoints(xy, proj4string=CRS("+init=EPSG:3067"))
  location<-as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
  lat <- location$y
  #print(paste("check crobas:",pCrobasX[55,3]))
  #print(pCrobasX)
  if((nYears*365)>ncol(clim$PAR)){
    clim$PAR <- cbind(clim$PAR,clim$PAR,clim$PAR,clim$PAR)
    clim$TAir <- cbind(clim$TAir,clim$TAir,clim$TAir,clim$PAR)
    clim$VPD <- cbind(clim$VPD,clim$VPD,clim$VPD,clim$PAR)
    clim$Precip <- cbind(clim$Precip,clim$Precip,clim$Precip,clim$PAR)
    clim$CO2 <- cbind(clim$CO2,clim$CO2,clim$CO2,clim$PAR)
    if(rcps=="CurrClim_fmi"){
      Tmm <- array(0,c(dim(TminTmax)[1],dim(TminTmax)[2]*4,2))
      Tmm[,,1] <- cbind(TminTmax[,,1],TminTmax[,,1],TminTmax[,,1],TminTmax[,,1])
      Tmm[,,2] <- cbind(TminTmax[,,2],TminTmax[,,2],TminTmax[,,2],TminTmax[,,2])
      assign("TminTmax", Tmm)  
      rm(list="Tmm")
      gc()
    }
    print(paste("length of clim:",ncol(clim$PAR),"versus dimension",nYears*365))
  }
  
  if(!is.na(P0currclim[1])){
    print("initialization with N module")
    #print(P0currclim)
    #save(nYears,nSites,siteInfo,lat,pCrobasX,parsCN_new_alfar,restrictionSwitch,                                
    #      defaultThin,
    #       ClCut, 
    #       areas,
    #       energyCut, 
    #       ftTapioParX,
    #       tTapioParX,
    #       initVar,
    #       clim,
    #       mortMod,
    #       P0currclim, fT0, file=paste0("testDataInit",restrictionSwitch,".rdata"))
    # print("data saved")
    if(length(clim$id)!=length(P0currclim)){
      P0currclim <- as.vector(mean(P0currclim)*array(1,c(1,length(clim$id))))
      fT0 <- as.vector(mean(fT0)*array(1,c(1,length(clim$id))))
    }
    
    initPrebas <- InitMultiSite(nYearsMS = rep(nYears,nSites),siteInfo=siteInfo,
                                latitude = lat,
                                pCROBAS = pCrobasX,
                                pCN_alfar = parsCN_new_alfar,
                                alpharNcalc = T,
                                ingrowth = ingrowth,
                                alpharVersion = restrictionSwitch,                                
                                ECMmod = 1,
                                defaultThin = defaultThin,
                                ClCut = ClCut, 
                                areas =areas,
                                energyCut = energyCut, 
                                ftTapioPar = ftTapioParX,
                                tTapioPar = tTapioParX,
                                multiInitVar = as.array(initVar),
                                PAR = clim$PAR[, 1:(nYears*365)],
                                TAir=clim$TAir[, 1:(nYears*365)],
                                VPD=clim$VPD[, 1:(nYears*365)],
                                Precip=clim$Precip[, 1:(nYears*365)],
                                CO2=clim$CO2[, 1:(nYears*365)],
                                yassoRun = 1,
                                mortMod = mortMod,
                                p0currClim = P0currclim, fT0AvgCurrClim = fT0,
                                TminTmax=TminTmax[,1:(nYears*365),],
                                disturbanceON = disturbanceON)
  } else {
    if(FALSE){save(nYears,nSites,siteInfo,sid,lat,pCrobasX,defaultThin,ClCut,areas,ingrowth,energyCut,
         ftTapioParX,tTapioParX,initVar,clim,mortMod,TminTmax,disturbanceON, 
         file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/testDataInit_master.rdata"))
        print("data saved")}
    print("run initPrebas")
    initPrebas <- InitMultiSite(nYearsMS = rep(nYears,nSites),siteInfo=siteInfo,
                                siteInfoDist = sid,
                                latitude = lat,
                                pCROBAS = pCrobasX,
                                ECMmod = 1,
                                defaultThin = defaultThin,
                                ClCut = ClCut, 
                                ingrowth = ingrowth,
                                areas =areas,
                                energyCut = energyCut, 
                                ftTapioPar = ftTapioParX,
                                tTapioPar = tTapioParX,
                                multiInitVar = as.array(initVar),
                                PAR = clim$PAR[, 1:(nYears*365)],
                                TAir=clim$TAir[, 1:(nYears*365)],
                                VPD=clim$VPD[, 1:(nYears*365)],
                                Precip=clim$Precip[, 1:(nYears*365)],
                                CO2=clim$CO2[, 1:(nYears*365)],
                                yassoRun = 1,
                                mortMod = mortMod,
                                TminTmax = TminTmax, 
                                disturbanceON = disturbanceON)
      
  }

}

yasso.mean.climate.f = function(dat, data.sample, startingYear, nYears){
  dat = dat[id %in% data.sample[, unique(id)]]
  dat[, DOY:=rep(1:365, len=dim(dat)[1])]
  dat[, Year:=rep(1980:2099, each=365)]
  #dat[, Year:= as.numeric(format(pvm, "%Y"))]
  dat = dat[Year >= startingYear & Year <= startingYear+nYears]
  dat[, pvm:= as.Date(paste(Year, '-01-01', sep="")) - 1 + DOY ]
  #dat[, DOY:= as.numeric(format(pvm, "%j"))]
  dat[, Mon:= as.numeric(format(pvm, "%m"))]
  #dat[DOY==366, DOY:=365]
  Tmean = dat[, mean(TAir), by = Year]
  Tsum = dat[, sum(ifelse(TAir>5, TAir-5, 0)), by=.(id, Year)][, mean(V1), by=Year]
  PAR = dat[, mean(PAR), by = Year]
  VPD = dat[, mean(VPD), by = Year]
  CO2 = dat[, mean(CO2), by = Year]
  Precip = dat[, sum(Precip), by = .(id, Year)][, mean(V1), by=Year]
  Tampl = dat[, .(mean(TAir)), by = .(id, Year, Mon)][, (max(V1)-min(V1))/2, by=Year]
  
  out = cbind(Tmean, Precip[, -1], Tampl[, -1], CO2[, -1], PAR[, -1], VPD[, -1], Tsum[, -1])
  colnames(out) = c('Year','Tmean','Precip','Tampl', 'CO2', "PAR", "VPD", "Tsum5")
  out
}


prep.climate.f = function(dat, data.sample, startingYear, nYears, rcps){
  print("Formulate climate data...")
  dat = dat[id %in% data.sample[, unique(id)]]
 # print(colnames(dat))
  
  if(rcps == "CurrClim_fmi"){
    colnames(dat)[which(colnames(dat)=="time")] <- "pvm"
  } else {
    dat[, pvm:= as.Date('1991-01-01') - 1 + rday ]
  }
  dat[, DOY:= as.numeric(format(pvm, "%j"))]
  dat[, Year:= as.numeric(format(pvm, "%Y"))]
  dat = dat[Year >= startingYear]
  dat[DOY==366, DOY:=365]
  # }
  id = dat[,unique(id)]
  PARtran = t( dcast(dat[, list(id, rday, PAR)], rday ~ id,
                     value.var="PAR")[, -1])
  TAirtran = t( dcast(dat[, list(id, rday, TAir)], rday ~ id,
                      value.var="TAir")[, -1])
  VPDtran = t( dcast(dat[, list(id, rday, VPD)], rday ~ id,
                     value.var="VPD")[, -1])
  Preciptran = t( dcast(dat[, list(id, rday, Precip)], rday ~ id,
                        value.var="Precip")[, -1])
  CO2tran = t( dcast(dat[, list(id, rday, CO2)], rday ~ id,
                     value.var="CO2")[, -1])
  print("...done.")
  
  list(PAR=PARtran, TAir=TAirtran, VPD=VPDtran, 
       Precip=Preciptran, CO2=CO2tran,id=id)
}


# simSummary.f = function(region=region, r_no, nYears, startingYear, rcpfile, harvScen) {
#   
#   out = region[['multiOut']]
#   VOL = out[, , 30, , 1]
#   VOL = apply(VOL, c(1,2), sum)
#   ## SO THIS IS NOW MEAN VOL PER HA OF nSample SIMULATED SAMPLES (by YEAR):
#   VOL = apply(VOL, 2, mean)
#   ## Multiply by area (tha)
#   VOL_INAREA = VOL * nfiareas[ID==r_no, AREA] * 1000 / 1000000 ## mill m3
#   ## at the beginning 207.7 mill m3, vrt 189.9 according to NFI (for region 7 = Keski-Suomi)
#   
#   Vmort = out[, , 42, , 1]
#   Vmort = apply(Vmort, c(1,2), sum)
#   ## SO THIS IS NOW MEAN VOL PER HA OF nSample SIMULATED SAMPLES (by YEAR):
#   Vmort = apply(Vmort, 2, mean)
#   Vmort_INAREA = Vmort * nfiareas[ID==r_no, AREA] * 1000 / 1000000 ## mill m3
#   
#   
#   ## WHY THIS IS NOT THE SAME AS har?
#   Vharvested = out[, , 37, , 1]
#   Vharvested = apply(Vharvested, c(1,2), sum)
#   ## SO THIS IS NOW MEAN VOL PER HA OF nSample SIMULATED SAMPLES (by YEAR):
#   Vharvested = apply(Vharvested, 2, mean)
#   Vharvested_INAREA = Vharvested * nfiareas[ID==r_no, AREA] * 1000 / 1000000 ## mill m3
#   
#   grossgrowth = out[, , 43, , 1]
#   grossgrowth = apply(grossgrowth, c(1,2), sum)
#   ## SO THIS IS NOW MEAN VOL PER HA OF nSample SIMULATED SAMPLES (by YEAR):
#   grossgrowth = apply(grossgrowth, 2, mean)
#   grossgrowth_INAREA = grossgrowth * nfiareas[ID==r_no, AREA] * 1000 / 1000000 ## mill m3
#   
#   dbh = out[, , 12, , 1]
#   dbh = apply(dbh, c(1,2), mean)
#   ## SO THIS IS NOW MEAN VOL PER HA OF nSample SIMULATED SAMPLES (by YEAR):
#   dbh = apply(dbh, 2, mean)
#   
#   age = out[, , 7, , 1]
#   age = apply(age, c(1,2), mean)
#   ## SO THIS IS NOW MEAN VOL PER HA OF nSample SIMULATED SAMPLES (by YEAR):
#   age = apply(age, 2, mean)
#   
#   gpp = out[, , 10, , 1]
#   gpp = apply(gpp, c(1,2), sum)
#   ## SO THIS IS NOW MEAN VOL PER HA OF nSample SIMULATED SAMPLES (by YEAR):
#   gpp = apply(gpp, 2, mean)
#   #npp_INAREA = npp * nfiareas[ID==7, AREA] * 1000 / 1000000 ## mill m3
#   
#   
#   npp = out[, , 18, , 1]
#   npp = apply(npp, c(1,2), sum)
#   ## SO THIS IS NOW MEAN VOL PER HA OF nSample SIMULATED SAMPLES (by YEAR):
#   npp = apply(npp, 2, mean)
#   #npp_INAREA = npp * nfiareas[ID==7, AREA] * 1000 / 1000000 ## mill m3
#   
#   
#   nep = out[, , 46, , 1]
#   nep = apply(nep, c(1,2), sum, na.rm=TRUE)
#   ## SO THIS IS NOW MEAN VOL PER HA OF nSample SIMULATED SAMPLES (by YEAR):
#   nep = apply(nep, 2, mean)
#   
#   
#   B_tree = out[, , 35, , 1]
#   B_tree = apply(B_tree, c(1,2), sum)
#   ## SO THIS IS NOW MEAN VOL PER HA OF nSample SIMULATED SAMPLES (by YEAR):
#   B_tree = apply(B_tree, 2, mean)
#   
#   lproj = out[, , 21, , 1]
#   lproj = apply(lproj, c(1,2), mean)
#   ## SO THIS IS NOW MEAN VOL PER HA OF nSample SIMULATED SAMPLES (by YEAR):
#   lproj = apply(lproj, 2, mean)
#   data.table(r_no, rcpfile, harvScen, year=startingYear + (1:nYears),
#              VOL, VOL_INAREA, Vharvested, Vmort, Vmort_INAREA,
#              grossgrowth_INAREA, dbh, age, gpp, npp, nep, B_tree, lproj)
# }
# 
# # this function create maps in tif format from raw data.
# createTif <- function(climate, management, yearOut, variable, species, startingYear){
#   simYear <- yearOut - startingYear
#   
#   files <- intersect(list.files(path= "output/", pattern = climate), list.files(path= "output/",pattern = management))
#   
#   outX <- data.table()
#   ops <- split(data.all, sample(1:115, nrow(data.all), replace=T))
#   
#   for(i in 1:length(files)){
#     sampleID <- paste0("sample",i,".")
#     
#     fileX <- files[grep(sampleID,files,fixed = T)]
#     
#     load(paste0("output/",fileX))
#     
#     out <- data.table(out$annual[,simYear,variable,])
#     
#     set.seed(1)
#     sampleX <- ops[[i]]
#     sampleX[,area := N*16^2/10000]
#     sampleX[,id:=climID]
#     
#     outX <- rbind(outX,cbind(sampleX$segID,out))
#     print(i)
#   }
#   
#   
#   
#   setnames(outX,c("segID","pine","spruce","birch"))
#   
#   outX[, tot := rowSums(.SD), .SDcols = c("pine","spruce","birch")]
#   
#   
#   outXY <- merge(kokeIDsTab,outX,all = T)
#   
#   ###create raster 
#   rastX <- rasterFromXYZ(outXY[,c("x","y",species),with=F])
#   crs(rastX) <- crs(kokeShp)
#   
#   rastName <- paste0("outRast/",climate,"_",management,"_var",varNames[variable],
#                      "_spec",species,"_year",yearOut,".tif")
#   writeRaster(rastX,filename = rastName)
# }
# 
# # this function create maps in tif format from data.tables selecting one year or the average of a time priod if yearOut is a vector of years
# createTifFromDT <- function(climate, management, yearOut, variable, species, startingYear){
#   simYear <- yearOut - startingYear
#   fileDT=paste0("outputDT/",varNames[variable],"_",management,"_",climate,".rdata")  
#   load(fileDT)
#   
#   outX <- t(get(varNames[variable]))
#   if (length(simYear)==1) outX <- outX[simYear,]
#   if (length(simYear)>1) outX <- colMeans(outX[simYear,],na.rm = T)
#   
#   segID <- areas <-numeric(0)
#   set.seed(1)
#   ops <- split(data.all, sample(1:115, nrow(data.all), replace=T))
#   for(i in 1:115){
#     # set.seed(1)
#     sampleX <- ops[[i]]
#     sampleX[,area := N*16^2/10000]
#     sampleX[,id:=climID]
#     segID <- c(segID,sampleX$segID)
#     areas <- c(areas,sampleX$area)
#     # print(i)
#   }
#   outX <- data.table(cbind(segID,areas,outX))
#   
#   setnames(outX,c("segID","areas",varNames[variable]))
#   
#   # outX[, tot := rowSums(.SD), .SDcols = c("pine","spruce","birch")]
#   
#   outXY <- merge(kokeIDsTab,outX,all = T)
#   
#   ###create raster 
#   rastX <- rasterFromXYZ(outXY[,c("x","y",varNames[variable]),with=F])
#   crs(rastX) <- crs(kokeShp)
#   
#   rastName <- paste0("outRast/",climate,"_",management,"_var",varNames[variable],
#                      "_spec",species,"_year",min(yearOut),"_",max(yearOut),".tif")
#   writeRaster(rastX,filename = rastName,overwrite=T)
# }
# 
# 
# 
# ##function to compile all data and create data.table 
# createDT <- function(climate, management,variable, species, startingYear){
#   
#   files <- intersect(list.files(path= "output/", pattern = climate), list.files(path= "output/",pattern = management))
#   
#   for (ij in variable) assign(varNames[ij],data.table())
#   VenergyWood <- WenergyWood <- data.table()
#   
#   # segID <- areas <-numeric(0)
#   
#   for(i in 1:length(files)){
#     sampleID <- paste0("sample",i,".")
#     
#     fileX <- files[grep(sampleID,files,fixed = T)]
#     
#     load(paste0("output/",fileX))
#     
#     ###sum harvests
#     if(i==1){
#       harvest <- out$harvest
#     }else{
#       harvest <- harvest+out$harvest  
#     }
#     
#     VenergyWood <- rbind(VenergyWood,apply(out$energyWood[,,,1],1:2,sum))
#     WenergyWood <- rbind(WenergyWood,apply(out$energyWood[,,,2],1:2,sum))
#     
#     marginX= 1:2#(length(dim(out$annual[,,variable,]))-1)
#     for (ij in variable) {
#      varIndx <- match(varNames[ij],varNames[varSel])  
#      assign(varNames[ij],data.table(rbind(eval(parse(text = varNames[ij])),
#                  apply(out$annual[,,varIndx,],marginX,sum))))
#     }    
#     print(i)
#   }
#   
#   ###proc and save total harvests
#   totHarvest <- data.table(harvest)
#   setnames(totHarvest,c("roundWood","energyWood"))
#   save(totHarvest,file=paste0("outputDT/","totHarvest","_",management,"_",climate,".rdata"))
#   
#   save(VenergyWood,file=paste0("outputDT/","VenergyWood","_",management,"_",climate,".rdata"))
#   save(WenergyWood,file=paste0("outputDT/","WenergyWood","_",management,"_",climate,".rdata"))
#   
#   
#   for(ij in variable) save(list=varNames[ij],file=paste0("outputDT/",varNames[ij],"_",management,"_",climate,".rdata"))
# }
# 
# ##function to compile all data and create data.table by species
# createDTbySp <- function(climate, management,variable, species, startingYear){
#   
#   files <- intersect(list.files(path= "output/", pattern = climate), list.files(path= "output/",pattern = management))
#   
#   for (ij in variable){
#     assign(paste0(varNames[ij],1),data.table())
#     assign(paste0(varNames[ij],2),data.table())
#     assign(paste0(varNames[ij],3),data.table())
#   }
#   # segID <- areas <-numeric(0)
#   
#   for(i in 1:length(files)){
#     sampleID <- paste0("sample",i,".")
#     
#     fileX <- files[grep(sampleID,files,fixed = T)]
#     
#     load(paste0("output/",fileX))
#     
#     marginX= 1:2#(length(dim(out$annual[,,variable,]))-1)
#     for (ij in variable){
#       assign(paste0(varNames[ij],1),
#              data.table(rbind(eval(parse(text = paste0(varNames[ij],1))),
#                               out$annual[,,ij,1])))
#       assign(paste0(varNames[ij],2),
#              data.table(rbind(eval(parse(text = paste0(varNames[ij],2))),
#                               out$annual[,,ij,2])))
#       assign(paste0(varNames[ij],3),
#              data.table(rbind(eval(parse(text = paste0(varNames[ij],3))),
#                               out$annual[,,ij,3])))
#     } 
#     
#     print(i)
#   }
#   
#   for(ij in variable){
#     save(list=c(paste0(varNames[ij],1),paste0(varNames[ij],2),paste0(varNames[ij],3)),
#          file=paste0("outputDT/",varNames[ij],"_",management,"_",climate,"_bySpecies.rdata"))
#   } 
# }
# 
# 
# # this function compute the annual totals of the region from data.tables 
# aTOTfromDT <- function(yearOut, variable, species="tot", startingYear){
#   simYear <- yearOut - startingYear
#   segID <- areas <-numeric(0)
#   set.seed(1)
#   ops <- split(data.all, sample(1:115, nrow(data.all), replace=T))
#   for(i in 1:115){
#     sampleX <- ops[[i]]
#     sampleX[,area := N*16^2/10000]
#     sampleX[,id:=climID]
#     segID <- c(segID,sampleX$segID)
#     areas <- c(areas,sampleX$area)
#     # print(i)
#   }
#   files <- list.files("outputDT/",pattern=paste0(varNames[variable],"_"))
#   if (species=="tot") files <- files[-grep("bySpecies",files)]
#   allOut <- data.table()
#   for(i in 1:length(files)){
#     load(paste0("outputDT/",files[i]))
#     
#     dats <- strsplit(files[i], "[_.]+")
#     if(variable %in% 32:33) dats[[1]] <- dats[[1]][-2]
#     climate=dats[[1]][3]
#     management=dats[[1]][2]
#     outX <- get(varNames[variable])*areas
#     outX <- colMeans(outX,na.rm=T)
#     outX <- data.table(cbind(outX,climate,management))
#     outX[,year:=yearOut]
#     allOut <- rbind(allOut,outX)
#   }
#   
#   setnames(allOut,c(varNames[variable],"climate","management","year"))
#   allOut$climate <- factor(allOut$climate)
#   allOut$management <- factor(allOut$management)
#   allOut[,1] <- as.numeric(unlist(allOut[,1]))
#   
#   fwrite(allOut,file= paste0("plots/",varNames[variable],"_DT.txt"))  
#   
#   p <- ggplot(data=allOut, 
#               aes_string(x="year", y=varNames[variable])) +
#     # scale_shape_manual(values=1:nlevels(countryTot$harvScenario)) +
#     labs(title = varNames[variable])+
#     # geom_smooth() +
#     xlab("Year") +
#     ylab("") +
#     geom_point(aes(colour=management, shape = climate,group=interaction(management, climate))) +
#     geom_line(aes(colour=management, group=interaction(management, climate)))
#   
#   pdf(file=paste0("plots/",varNames[variable], ".pdf"))
#   print(p)
#   dev.off()
# }
# 
# ###compute total biomass from DTs
# Wtot <- function(manClim){
#   files <- paste0("outputDT/",varNames[c(24:25,31:33)],manClim)
#   for(i in 1:5) load(files[i])
#   Wtot <- Wstem + W_croot + wf_STKG + Wbranch + WfineRoots
#   save(Wtot,file = paste0("outputDT/","Wtot",manClim))
# }


# createPlotfromDT <- function(path, variable){
#   DT <- fread(paste0(path,varNames[variable],"_DT.txt"))
#   p <- ggplot(data=DT, 
#               aes_string(x="year", y=varNames[variable])) +
#     # scale_shape_manual(values=1:nlevels(countryTot$harvScenario)) +
#     labs(title = varNames[variable])+
#     # geom_smooth() +
#     xlab("Year") +
#     ylab("") +
#     geom_point(aes(colour=management, shape = climate,group=interaction(management, climate))) +
#     geom_line(aes(colour=management, group=interaction(management, climate)))
#   
#   png(file=paste0("plots/",varNames[variable], ".png"),width = 500,height = 500)
#   print(p)
#   dev.off()
# }

calMean <- function(varX,hscenX,areas){
  load(paste0("outputDT/",varX,"_",hscenX,"_CurrClim.rdata"))
  varAreas <- get(varX)*areas
  # Vareas <- Vareas[-siteX]
  totX <- colSums(varAreas,na.rm = T)
  meanX <- totX/sum(areas)#co
  return(meanX)
}

calculatePerCols <- function(outX){ #perStarts,perEnds,startingYear,
  iper <- 1
  for(iper in 1:length(perStarts)){      
    per <- perStarts[iper]:perEnds[iper]
    simYear = per - startingYear# + 1
    colsOut = c(paste("V", simYear, sep=""))
#    outX <- outX*area
    p <- outX[, .(per = rowMeans(.SD,na.rm=T)), .SDcols = colsOut, by = segID] 
    colnames(p)[2] <- paste0("per",iper)
    if(iper==1) {
      pX <- data.table(p)
  #    colnames(pX)[1] <- "var"
    } else {
      pX <- cbind(pX, p[,2])
    }
  }
  return(pX)
}

calculatePerColsAllRows <- function(outX){ #perStarts,perEnds,startingYear,
  for(iper in 1:length(perStarts)){      
    per <- perStarts[iper]:perEnds[iper]
    simYear = per - startingYear# + 1
    colsOut = c(paste("V", simYear, sep=""))
    p <- cbind(outX[,"segID"],rowMeans(outX[,..colsOut])) 
    colnames(p)[2] <- paste0("per",iper)
    if(iper==1) {
      pX <- data.table(p)
  #    colnames(pX)[1] <- "var"
    } else {
      pX <- cbind(pX, p[,2])
    }
  }
  return(pX)
}

specialVarProcAdapt <- function(sampleX,region,r_no,harvScen,harvInten,rcpfile,sampleID,
                           areas,sampleForPlots,output, toRaster){#},SBBbp,PI,pSBB){
  nYears <-  max(region$nYears)
  nSites <-  max(region$nSites)
  ####process and save special variables: 
  ###dominant Species
  outX <- domFun(region,varX="species")  
  ####test plot
  #if(sampleID==sampleForPlots){testPlot(outX,"domSpecies",areas)}
  ###take the most frequent species in the periods
  pX <- calculatePerCols(outX = outX)
  varNam <- "domSpecies"
  assign(varNam,pX)
  if(toRaster){
      save(list=varNam,
         file=paste0(path_output,"weatherStation",station_id,"/",
                     varNam,
                     "_harscen",harvScen,
                     "_harInten",harvInten,"_",
                     rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  rm(list=varNam); gc()
  #print(output)
  
  # rm(domSpecies); gc()
  ###age dominant species
  outX <- domFun(region,varX="age")
  pX <- calculatePerCols(outX = outX)
  varNam <- "domAge"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  rm(list=varNam); gc()
  #print(output)

  ### pine Volume Vpine
  outX <- vSpFun(region,SpID=1)
  #outX <- vDecFun(region)
  pX <- calculatePerCols(outX = outX)
  varNam <- "Vpine"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  #print(output)
  
  ### spruce Volume Vspruce
  outX <- vSpFun(region,SpID=2)
  pX <- calculatePerCols(outX = outX)
  varNam <-  "Vspruce"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  #print(output)
  
  ### deciduous Volume Vdec
  outX <- vSpFun(region,SpID=3)
  pX <- calculatePerCols(outX = outX)
  varNam <-  "Vdec"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  #print(output)
  

  ####WenergyWood
  outX <- data.table(segID=sampleX$segID,apply(region$multiEnergyWood[,,,2],1:2,sum))
  pX <- calculatePerCols(outX = outX)
  varNam <-  "Wenergywood"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  #print(output)

  ####VenergyWood
  outX <- data.table(segID=sampleX$segID,apply(region$multiEnergyWood[,,,1],1:2,sum))
  pX <- calculatePerCols(outX = outX)
  varNam <-  "Venergywood"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  #print(output)

  ####GVgpp
  outX <- data.table(segID=sampleX$segID,region$GVout[,,3])
  pX <- calculatePerCols(outX = outX)
  varNam <-  "GVgpp"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  #print(output)
  
  ####GVw
  outX <- data.table(segID=sampleX$segID,region$GVout[,,4])
  pX <- calculatePerCols(outX = outX)
  varNam <-  "GVw"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  #print(output)
  
  ####Wtot
  outX <- data.table(segID=sampleX$segID,apply(region$multiOut[,,c(24,25,31,32,33),,1],1:2,sum))
  pX <- calculatePerCols(outX = outX)
  varNam <-  "Wtot"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  #print(output)
  
  NUP <- T
  if(exists("parsCN_alfar")){
    #### alphar
    outX <- data.table(segID=sampleX$segID,region$multiOut[,,3,1,2])
    pX <- calculatePerCols(outX = outX)
    varNam <-  "alphar"
    assign(varNam,pX)
    if(toRaster){
      save(list=varNam,
         file=paste0(path_output,"weatherStation",station_id,"/",
                     varNam,
                     "_harscen",harvScen,
                     "_harInten",harvInten,"_",
                     rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
    }
    pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
    pX <- c(var = varNam, pX)
    output <- rbind(output, pX)
    colnames(output) <- names(pX)
    #print(output)
    
    ####### Nup
#    outX <- data.table(segID=sampleX$segID,region$multiOut[,,55,1,2])
    marginX <- 1:2
      outX <- data.table(segID=sampleX$segID,apply(region$multiOut[,,55,,2],marginX,sum))
    pX <- calculatePerCols(outX = outX)
    varNam <- "Nup"
    assign(varNam,pX)
    if(toRaster){
      save(list=varNam,
         file=paste0(path_output,"weatherStation",station_id,"/",
                     varNam,
                     "_harscen",harvScen,
                     "_harInten",harvInten,"_",
                     rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
    }
    pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
    pX <- c(var = varNam, pX)
    output <- rbind(output, pX)
    colnames(output) <- names(pX)
    #print(output)
    
    ### Ndem
#    outX <- data.table(segID=sampleX$segID,region$multiOut[,,56,1,2])
      outX <- data.table(segID=sampleX$segID,apply(region$multiOut[,,56,,2],marginX,sum))
    pX <- calculatePerCols(outX = outX)
    varNam <- "Ndem"
    assign(varNam,pX)
    if(toRaster){
      save(list=varNam,
         file=paste0(path_output,"weatherStation",station_id,"/",
                     varNam,
                     "_harscen",harvScen,
                     "_harInten",harvInten,"_",
                     rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
    }
    pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
    pX <- c(var = varNam, pX)
    output <- rbind(output, pX)
    colnames(output) <- names(pX)
    #print(output)
    
    ### "Umax"
    #outX <- data.table(segID=sampleX$segID,region$multiOut[,,57,1,2])
    tmp <- region$multiOut[,,57,,1]
    region$multiOut[,,57,,1] <- region$multiOut[,,57,,2] 
    outX <- data.table(segID=sampleX$segID,baWmean(region,57))
    region$multiOut[,,57,,1] <- tmp
    pX <- calculatePerCols(outX = outX)
    varNam <- "Umax"
    assign(varNam,pX)
    if(toRaster){
      save(list=varNam,
         file=paste0(path_output,"weatherStation",station_id,"/",
                     varNam,
                     "_harscen",harvScen,
                     "_harInten",harvInten,"_",
                     rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
    }
    pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
    pX <- c(var = varNam, pX)
    output <- rbind(output, pX)
    colnames(output) <- names(pX)
    #print(output)
    

    if(FALSE){
    ####### "Gf"
    outX <- data.table(segID=sampleX$segID,apply(region$multiOut[,,55,,1],marginX,sum))
    pX <- calculatePerCols(outX = outX)
    varNam <- "Gf"
    assign(varNam,pX)
    if(toRaster){
      save(list=varNam,
         file=paste0(path_output,"weatherStation",station_id,"/",
                     varNam,
                     "_harscen",harvScen,
                     "_harInten",harvInten,"_",
                     rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
    }
    pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
    pX <- c(var = varNam, pX)
    output <- rbind(output, pX)
    colnames(output) <- names(pX)
    #print(output)
    
    ### "Gr"
    outX <- data.table(segID=sampleX$segID,apply(region$multiOut[,,56,,1],marginX,sum))
#    outX <- data.table(segID=sampleX$segID,region$multiOut[,,56,1,1])
    pX <- calculatePerCols(outX = outX)
    varNam <- "Gr"
    assign(varNam,pX)
    if(toRaster){
      save(list=varNam,
         file=paste0(path_output,"weatherStation",station_id,"/",
                     varNam,
                     "_harscen",harvScen,
                     "_harInten",harvInten,"_",
                     rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
    }
    pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
    pX <- c(var = varNam, pX)
    output <- rbind(output, pX)
    colnames(output) <- names(pX)
    #print(output)
    
    ### "Gw"
    outX <- data.table(segID=sampleX$segID,apply(region$multiOut[,,57,,1],marginX,sum))
#    outX <- data.table(segID=sampleX$segID,region$multiOut[,,57,1,1])
    pX <- calculatePerCols(outX = outX)
    varNam <- "Gw"
    assign(varNam,pX)
    if(toRaster){
      save(list=varNam,
         file=paste0(path_output,"weatherStation",station_id,"/",
                     varNam,
                     "_harscen",harvScen,
                     "_harInten",harvInten,"_",
                     rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
    }
    pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
    pX <- c(var = varNam, pX)
    output <- rbind(output, pX)
    colnames(output) <- names(pX)
    #print(output)
    }
  }
  
  #### SBBprob
  outX <- data.table(segID=sampleX$segID,region$multiOut[,,45,1,2])
  pX <- calculatePerCols(outX = outX)
  varNam <-  "SBBprob"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  #print(output)
  
  #### SMI
  outX <- data.table(segID=sampleX$segID,region$multiOut[,,46,1,2])
  pX <- calculatePerCols(outX = outX)
  varNam <-  "SMI"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  
  ## BB intensity
  #xSMI <- region$multiOut[,,46,1,2]
  #BA <- apply(region$multiOut[,,13,,1],1:2,sum)
  ##  BAspruce <- BASpFun(modOut = region,SpID = 2)[,-1]
  #xBAspruceFract <- BASpFun(modOut = region,SpID = 2)[,-1]/BA
  #xBAspruceFract[BA==0] <- 0
  #SHI = xBAspruceFract*(1-xSMI)/0.2093014
  #INTENSITY <- 1/(1+exp(3.9725-2.9673*SHI))
  #INTENSITY[xBAspruceFract<0.05] <- 0
  #pX <- calculatePerCols(outX = data.table(segID=sampleX$segID,INTENSITY))
  #varNam <-  "BBintensity"
  #assign(varNam,pX)
  #if(toRaster){
  #  save(list=varNam,
  #       file=paste0(path_output,"weatherStation",station_id,"/",
  #                   varNam,
  #                   "_harscen",harvScen,
  #                   "_harInten",harvInten,"_",
  #                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  #}
  #pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  #pX <- c(var = varNam, pX)
  #output <- rbind(output, pX)
  #colnames(output) <- names(pX)
  
  ## BB expected damage area
  #SBBprob <- region$multiOut[,,45,1,2]
  #SBBdamArea <- SBBprob*INTENSITY#*sampleX$area
  #pX <- calculatePerCols(outX = data.table(segID=sampleX$segID,SBBdamArea))
  #varNam <- "ExpectedBBdamAreaFraction"
  #pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)*100
  #pX <- c(var = varNam, pX)
  #output <- rbind(output, pX)
  #colnames(output) <- names(pX)
  
  ## BB simulated damage area
  SBBReactionBA <-  apply(region$multiOut[,,"grossGrowth/bb BA disturbed",,2],1:2,sum)
  BA <- apply(region$multiOut[,,"grossGrowth/bb BA disturbed",,1],1:2,sum)
  Vrw <- apply(region$multiOut[,,"VroundWood",,1],1:2,sum)[,-1]
  Vrw <- cbind(Vrw,Vrw[,ncol(Vrw)])
  areaSample <- array(areas,c(dim(SBBReactionBA))) # Segment areas where damage happened
  areaSample[SBBReactionBA==0] <- 0
  #if(clcut==-1){ # if no clearcut, calculate only the 
  areaSample[SBBReactionBA>0 & Vrw==0] <- areaSample[SBBReactionBA>0 & Vrw==0]*
    SBBReactionBA[SBBReactionBA>0 & Vrw==0]/BA[SBBReactionBA>0 & Vrw==0]
  areaSample[BA==0] <- 0
  #}
  #pX <- calculatePerCols(outX = data.table(segID=sampleX$segID, SBBReactionBA))
  pX <- calculatePerCols(outX = data.table(segID=sampleX$segID, areaSample))
  varNam <- "simBBdamArea%"
  pX <- colSums(pX[,-1])/sum(sampleX$area)*100
  #pX <- colSums(pX[,-1])
  #pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)*100
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  
  #### pFire
  outX <- data.table(segID=sampleX$segID,region$multiOut[,,47,1,2])
  pX <- calculatePerCols(outX = outX)
  varNam <-  "pFire"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
       file=paste0(path_output,"weatherStation",station_id,"/",
                   varNam,
                   "_harscen",harvScen,
                   "_harInten",harvInten,"_",
                   rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  
  # pWind
  outX <- data.table(segID=sampleX$segID,region$outDist[,,"wrisk"])
  pX <- calculatePerCols(outX = outX)
  varNam <-  "pWind"
  assign(varNam,pX)
  if(toRaster){
    save(list=varNam,
         file=paste0(path_output,"weatherStation",station_id,"/",
                     varNam,
                     "_harscen",harvScen,
                     "_harInten",harvInten,"_",
                     rcpfile,"_Nswitch",restrictionSwitch,".rdata"))
  }
  pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  
  ## BB simulated damage area
  WindReactionV <-  region$outDist[,,"damvol"]
  V <- apply(region$multiOut[,,"V",,1],1:2,sum)
  WindReactionSalvLog <- region$outDist[,,"salvlog"]
  areaSample <- array(areas,c(dim(WindReactionV))) # Segment areas where damage happened
#  areaSample[WindReactionSalvLog==0] <- 0 # if no wind damage, set to zero
  areaSample[WindReactionV==0 & WindReactionSalvLog==0] <- 0 # if no wind damage, set to zero
  # if no salvage logging but damage, Vdamage/V*area
  areaSample[WindReactionV>0 & WindReactionSalvLog==0] <- areaSample[WindReactionV>0 & WindReactionSalvLog==0]*
    WindReactionV[WindReactionV>0 & WindReactionSalvLog==0]/V[WindReactionV>0 & WindReactionSalvLog==0]
  areaSample[V==0] <-0
  #pX <- calculatePerCols(outX = data.table(segID=sampleX$segID, SBBReactionBA))
  pX <- calculatePerCols(outX = data.table(segID=sampleX$segID, areaSample))
  varNam <- "simWinddamArea%"
  pX <- colSums(pX[,-1])/sum(sampleX$area)*100
  #pX <- colSums(pX[,-1])
  #pX <- colSums(pX[,-1]*matrix(sampleX$area,nrow(pX),ncol(pX)-1))/sum(sampleX$area)*100
  pX <- c(var = varNam, pX)
  output <- rbind(output, pX)
  colnames(output) <- names(pX)
  
  
    ####SBBbp
#  outX <- data.table(segID=sampleX$segID,SBBbp)
#  pX <- calculatePerCols(outX = outX)
#  pX <- colMeans(pX)
#  pX[1] <- "SBBbp"
#  output <- rbind(output, pX)
#  colnames(output) <- names(pX)
  
  ####PI
#  outX <- data.table(segID=sampleX$segID,PI)
#  pX <- calculatePerCols(outX = outX)
#  pX <- colMeans(pX)
#  pX[1] <- "sbbPI"
#  output <- rbind(output, pX)
#  colnames(output) <- names(pX)

  ####pSBB damage
#  outX <- data.table(segID=sampleX$segID,pSBB)
#  pX <- calculatePerCols(outX = outX)
#  pX <- colMeans(pX)
#  pX[1] <- "pSBBdamage"
#  output <- rbind(output, pX)
#  colnames(output) <- names(pX)

  gc()
  
  return(output)
} 


BASpFun <- function(modOut,SpID){
  segID <- modOut$siteInfo[,1]
  oo <- data.table(which(modOut$multiOut[,,4,,1]==SpID,arr.ind=T))
  setnames(oo,c("site","year","layer"))
  vx <-modOut$multiOut[,,13,,1][as.matrix(oo)]
  oo$VSp <- vx
  setkey(oo,site,year)
  ff <- oo[,sum(VSp),by=.(site,year)]
  VspMat <- matrix(0,modOut$nSites,modOut$maxYears)
  VspMat[as.matrix(ff[,1:2])] <- unlist(ff[,3])
  outX <- data.table(segID=segID,VspMat)
  return(outX)
}



####test plot
testPlot <- function(outX,titleX,areas){
  cc <- data.table(rbind(cbind(1:nYears,apply(outX[,2:(nYears+1)],2,min,na.rm=T),"min"),
                         cbind(1:nYears,apply(outX[,2:(nYears+1)],2,max,na.rm=T),"max"),
                         cbind(1:nYears,apply(outX[,2:(nYears+1)],2,median,na.rm=T),"median"),
                         cbind(1:nYears,apply(outX[,2:(nYears+1)],2,mean,na.rm=T),"aritMean"),
                         cbind(1:nYears,apply((outX[,2:(nYears+1)]*areas/sum(areas)),2,sum,na.rm=T),"regionMean")))
  setnames(cc,c("simYear","value","metric"))
  # cc$metric=as.factor(cc$metric)
  cc$metric=factor(cc$metric)
  cc$value=as.double(cc$value)
  cc$simYear <- as.double(cc$simYear)
  cc <- cc[order(simYear)]
  testP <- ggplot(data=cc, aes(x=simYear, y=value, col=metric,group=metric)) +
    geom_line()+
    geom_point() + ggtitle(titleX)
  print(testP)
}


####Function to process NEP for drained peatlands (used in 2.1_procNep.r)
processPeat <- function(peatXf, fertf, npp_lit, nepf, peatval, fertval) {
  # peatXf = raster with peat soils
  # fertf = raster with soilType
  # npp_lit = raster of npp - litterfall (NEP= NPP - coeffSoil - lit)
  # nepf= raster with nep
  # peatval = ID to identify the drained peatlands -> tells which peat soil you want to treat
  # fertval = soilType ID -> tells which siteType you want to treat
  
  # rasters may be off by a couple pixels, resize:
  if (any(dim(fertf) < dim(peatXf))) {peatXf <- crop(peatXf,fertf)} 
  if (any(dim(peatXf) < dim(fertf))) {fertf <- crop(fertf,peatXf)}
  if (any(dim(fertf) < dim(npp_lit))) {npp_lit <- crop(npp_lit,fertf)} 
  if (any(dim(peatXf) < dim(npp_lit))) {npp_lit <- crop(npp_lit,peatXf)}
  if (any(dim(fertf) < dim(nepf))) {nepf <- crop(nepf,fertf)} 
  if (any(dim(peatXf) < dim(nepf))) {nepf <- crop(nepf,peatXf)}
  # mask out pixels where peatXf == peatval and fertx == fertval
  drPeatNeg <- peatXf == peatval & fertf == fertval  ###selecting the pixels that match the conditions of peat and siteType
  drPeatNeg[drPeatNeg==0] <- NA  ### assign NA to the remaining pixels
  drPeat <- mask(npp_lit, drPeatNeg)  ###raster with only the pixel of interest
  
  ###calculate the new NEP according to the siteType (fertval)
  if (fertval < 3) {         
    drPeat <- drPeat - 240  
  } else if (fertval >= 3) {
    drPeat <- drPeat + 70
  }
  return(merge(drPeat,nepf))
}



#####functions to calculate Mortality related metrics as in 
###15Silva Fennica vol. 54 no. 5 article id 10414 ?? Siipilehto et al. ?? Stand-level mortality models for Nordic boreal ...
pMort <- function(modOut,ageClass, rangeYear=5){
  endX <- rangeYear:dim(modOut)[2]
  startX <- endX-(rangeYear-1)
  pMortX <- rep(0.,length(endX))
  
  for(i in 1:length(startX)){
    ageX <-rowMeans(modOut[,startX[i]:endX[i],7,1,1])
    cX <- which(ageX %in% ageClass)
    # outX <- modOut[cX,,,,]
    mortX <- data.table(which(modOut[cX,startX[i]:endX[i],42,,1]>0,arr.ind=T))
    nMort <- length(unique(mortX$site))
    pMortX[i] <- nMort/length(cX)
  }
  return(pMortX)
}

###Function to calculate the probability of a mortality (pM) event occuring
# Arguments: 
# modOut = output array from a PREBAS multisite runs: $multiOut
# rangeYear = number of years  for which to calculate pM
# sp = species/layer for which to calculate pM it can be a vector for combinations of species
# pureFor = proportion of Basal area to consider as pure stands
# mixFor = it works only for mixed forests, it is the minimum proportion of basal area for the species of interest

pMort2 <- function(modOut,ageClass, rangeYear=5,sp,pureFor,mixFor){
  endX <- rangeYear:dim(modOut)[2]
  startX <- endX-(rangeYear-1)
  pMortX <- nSites <- rep(0.,length(endX))
  
  for(i in 1:length(startX)){
    ageX <-rowMeans(modOut[,startX[i]:endX[i],7,1,1])
    pBA <- apply(modOut[,startX[i]:endX[i],13,,1],c(1,3),mean)
    pBA <- pBA/rowSums(pBA)
    if(length(sp)==1){
      selX <- which(ageX %in% ageClass & pBA[,sp]>pureFor)
    }else{
      selX <- which(ageX %in% ageClass & rowSums(pBA[,sp])>mixFor &
                      pBA[,1]<pureFor & pBA[,2]<pureFor)  
    }
    
    # outX <- modOut[cX,,,,]
    mortX <- data.table(which(modOut[selX,startX[i]:endX[i],42,,1]>0,arr.ind=T))
    nMort <- length(unique(mortX$site))
    pMortX[i] <- nMort/length(selX)
    nSites[i] <- length(selX)
  }
  return(list(pMort=pMortX,nSites=nSites))
}


###function to calculate the mortality probability along some variable classes
# Arguments: 
# modOut = output array from a PREBAS multisite runs
# rangeYear = number of years  for which to calculate pM
# minX = minimum value for the variable class
# maxX = maximum value for the variable class
# stepX = class step
# varX = variable ID of PREBAS output (see varNames)
# funX = function to use to aggregate the data (mean or sum) mean for age and DBH, sum for BA, stemNumber
pMortVarX <- function(modOut,minX,maxX,stepX,varX,funX,rangeYear=5){
  endX <- rangeYear:dim(modOut)[2]
  startX <- endX-(rangeYear-1)
  seqX <- seq(minX,maxX,by=stepX)
  nClass <- length(seqX)+1
  pMortX <- nData <- matrix(0.,length(endX),nClass)
  for(i in 1:length(startX)){
    varXs<-apply(modOut[,startX[i]:endX[i],varX,,1],1:2,funX)
    varXs <- rowMeans(varXs)
    for(ij in 1:nClass){
      if(ij==1) cX <- which(varXs <= seqX[ij])
      if(ij>1 & ij<nClass) cX <- which(varXs <= seqX[ij] & varXs > seqX[ij-1])
      if(ij==nClass) cX <- which(varXs > seqX[ij-1])
      # outX <- modOut[cX,,,,]
      if(length(cX)>0.){
        mortX <- data.table(which(modOut[cX,startX[i]:endX[i],42,,1]>0,arr.ind=T))
        nMort <- length(unique(mortX$site))
        nData[i,ij] <- length(cX)
        pMortX[i,ij] <- nMort/length(cX)
      }
    }
  }
  return(list(pMort=pMortX,nData=nData,classes=seqX))
}


###function to calculate basal area of dead trees along some variable classes
# Arguments: 
# modOut = output array from a PREBAS multisite runs: $multiOut
# rangeYear = number of years  for which to calculate pM
# minX = minimum value for the variable class
# maxX = maximum value for the variable class
# stepX = class step
# varX = variable ID of PREBAS output (see varNames)
# funX = function to use to aggregate the data (mean or sum) mean for age and DBH, sum for BA, stemNumber
baMortVarX <- function(modOut,minX,maxX,stepX,varX,funX,rangeYear=5){
  nYears <- dim(modOut)[2]
  nMort <- modOut[,2:nYears,42,,1]/modOut[,1:(nYears-1),30,,1]*modOut[,1:(nYears-1),17,,1]
  nMort[which(is.na(nMort))] <- 0.
  baMort <- nMort * modOut[,1:(nYears-1),35,,1]
  baTot <- apply(modOut[,1:(nYears-1),13,,1],1:2,sum)
  
  endX <- rangeYear:(nYears-1)
  startX <- endX-(rangeYear-1)
  seqX <- seq(minX,maxX,by=stepX)
  nClass <- length(seqX)+1
  baTotX <- baMortX <- nData <- matrix(0.,length(endX),nClass)
  # oo <- modOut
  modOut <- modOut[,2:nYears,,,]
  for(i in 1:length(startX)){
    varXs<-apply(modOut[,startX[i]:endX[i],varX,,1],1:2,funX)
    varXs <- rowMeans(varXs)
    for(ij in 1:nClass){
      if(ij==1) cX <- which(varXs <= seqX[ij])
      if(ij>1 & ij<nClass) cX <- which(varXs <= seqX[ij] & varXs > seqX[ij-1])
      if(ij==nClass) cX <- which(varXs > seqX[ij-1])
      # outX <- modOut[cX,,,,]
      if(length(cX)>0.){
        baX <- sum(baMort[cX,startX[i]:endX[i],])/length(cX)
        baTx <- sum(baTot[cX,startX[i]:endX[i]])/rangeYear/length(cX)
        nData[i,ij] <- length(cX)
        baMortX[i,ij] <- baX
        baTotX[i,ij] <- baTx
      }
    }
  }
  return(list(baMort=baMortX,nData=nData,classes=seqX,
              baTot=baTotX))
}


###function to calculate the mortality probability for species proportion
# Arguments: 
# modOut = output array from a PREBAS multisite runs
# rangeYear = number of years  for which to calculate pM
# minX = minimum species cover
# maxX = maximum species cover
# stepX = class step
pMortSpecies <- function(modOut,minX=0.1,maxX=0.9,stepX=0.1,rangeYear=5){
  endX <- rangeYear:dim(modOut)[2]
  startX <- endX-(rangeYear-1)
  seqX <- seq(minX,maxX,by=stepX)
  nClass <- length(seqX)+1
  pMortXpine <- nDataPine <- 
    pMortXspruce <- nDataSpruce <- 
    pMortXbirch <- nDataBirch <- matrix(0.,length(endX),nClass)
  totBA <- apply(modOut[,,13,,1],1:2,sum)
  pBApine <- modOut[,,13,1,1]/totBA
  pBAspruce <- modOut[,,13,2,1]/totBA
  pBAbirch <- modOut[,,13,3,1]/totBA
  for(i in 1:length(startX)){
    subPine <-rowMeans(pBApine[,startX[i]:endX[i]],na.rm=T)
    subSpruce <-rowMeans(pBAspruce[,startX[i]:endX[i]],na.rm=T)
    subBirch <-rowMeans(pBAbirch[,startX[i]:endX[i]],na.rm=T)
    for(ij in 1:nClass){
      if(ij==1){
        cXpine <- which(subPine <= seqX[ij])
        cXspruce <- which(subSpruce <= seqX[ij])
        cXbirch <- which(subBirch <= seqX[ij])
      } 
      if(ij>1 & ij<nClass){
        cXpine <- which(subPine <= seqX[ij] & subPine > seqX[ij-1])
        cXspruce <- which(subSpruce <= seqX[ij] & subSpruce > seqX[ij-1])
        cXbirch <- which(subBirch <= seqX[ij] & subBirch > seqX[ij-1])
      } 
      if(ij==nClass){
        cXpine <- which(subPine > seqX[ij])
        cXspruce <- which(subSpruce > seqX[ij])
        cXbirch <- which(subBirch > seqX[ij])
      } 
      # outX <- modOut[cX,,,,]
      if(length(cXpine)>0.){
        mortX <- data.table(which(modOut[cXpine,startX[i]:endX[i],42,,1]>0,arr.ind=T))
        nMort <- length(unique(mortX$site))
        nDataPine[i,ij] <- length(cXpine)
        pMortXpine[i,ij] <- nMort/length(cXpine)
      }
      if(length(cXspruce)>0.){
        mortX <- data.table(which(modOut[cXspruce,startX[i]:endX[i],42,,1]>0,arr.ind=T))
        nMort <- length(unique(mortX$site))
        nDataSpruce[i,ij] <- length(cXspruce)
        pMortXspruce[i,ij] <- nMort/length(cXspruce)
      }
      if(length(cXbirch)>0.){
        mortX <- data.table(which(modOut[cXbirch,startX[i]:endX[i],42,,1]>0,arr.ind=T))
        nMort <- length(unique(mortX$site))
        nDataBirch[i,ij] <- length(cXbirch)
        pMortXbirch[i,ij] <- nMort/length(cXbirch)
      }
    }
  }
  return(list(pMortPine=pMortXpine,nDataPine=nDataPine,
              pMortSpruce=pMortXspruce,nDataSpruce=nDataSpruce,
              pMortBirch=pMortXbirch,nDataBirch=nDataBirch))
}


#### function to calculate the new parameters of the ClearCuts
#### increasing the rotation length
#### out=multi prebas run output
calNewDclcut <- function(out,
                         ClCut_pine,
                         ClCut_spruce,
                         ClCut_birch,
                         fact=0.25){
  newClCut_pine <- ClCut_pine
  newClCut_spruce <- ClCut_spruce 
  newClCut_birch <- ClCut_birch
  
  nSites <- dim(out$multiOut)[1]
  ETSmean = round(mean(out$multiOut[,,5,1,1]))
  domSp <- rep(NA,nSites)
  domX <- apply(out$multiOut[,,13,,1],1:2, which.max)
  domPos <- apply(domX,1,FUN=function(x) which.max(table(x)))
  for(i in 1:nSites) domSp[i] <- out$multiOut[i,1,4,domPos[i],1]
  siteType <- out$siteInfo[,3]
  
  sitesP3 <- which(siteType<=3 & domSp==1)
  sitesP4 <- which(siteType==4 & domSp==1)
  sitesP5 <- which(siteType>=5 & domSp==1)
  sitesSP2 <- which(siteType<=2 & domSp==2)
  sitesSP3 <- which(siteType>=3 & domSp==2)
  sitesB2 <- which(siteType<=2 & domSp==3)
  sitesB3 <- which(siteType>=3 & domSp==3)
  
  pdf(paste0(pathX,"ClCutplots_maak",r_no,".pdf"))
  for(j in 1:7){
    sites <- get(spSite[j])
    spX <- spXs[j]
    dClcut <- get(tabX[j])[indX[j],c(1,3)]
    aClcut <- get(tabX[j])[indX[j],c(2,4)]
    
    dataX <- data.table(age=as.vector(out$multiOut[sites,,7,spX,1]),
                        d=as.vector(out$multiOut[sites,,12,spX,1]))
    
    dataX <- dataX[age>0. & d>0]
    
    fitMod = nlsLM(d ~ a*(1-exp(b*age))^c,
                   start = list(a = 60,b=-0.07,c=1.185),
                   data = dataX)
    
    modD <- data.table(age=seq(0,300,0.5))    
    modD$d=predict(fitMod, list(age = modD$age))
    
    px <- coef(fitMod)
    a=px[1];b=px[2];c=px[3]
    # 
    dd=dClcut
    aa= log(1-(dd/a)^(1/c))/b
    if(any(is.na(aa))) aa[is.na(aa)] <- aClcut[is.na(aa)]
    predict(fitMod, list(age = aa))
    
    if(fact<2 & fact>0.){
      age2 <- (1+fact) * aa
    }else{
      age2 <- aa + fact 
    } 
    d2 <- predict(fitMod, list(age = age2))
    d2
    
    dataX[,plot(age,d,pch='.',ylim=c(0,45),xlim=c(0,300))]
    lines(modD$age,modD$d,col=4)
    points(aa,dd,col=3,pch=c(1,20))
    points(age2,d2,col=2,pch=c(1,20))
    abline(v=aClcut,col=3,lty=1:2)
    abline(v=aClcut*1.25,col=2,lty=1:2)
    legend("bottomright",cex=0.8,
           pch=c(1,1,1,20,1,NA),
           legend=c("standard","+25%",
                    "ETs<1000","ETS>1000",
                    "D","age"),
           col= c(3,2,1,1,1,1),
           lty=c(NA,NA,NA,NA,NA,1)
    )
    legend("topleft",cex=0.5,
           c(paste0("ETSmean = ",ETSmean))
    )
    
    if(tabX[j]=="ClCut_pine") newClCut_pine[indX[j],c(1,3)] <- d2
    if(tabX[j]=="ClCut_spruce") newClCut_spruce[indX[j],c(1,3)] <- d2
    if(tabX[j]=="ClCut_birch") newClCut_birch[indX[j],c(1,3)] <- d2
    
    print(j)
  }
  dev.off()
  if(fact<2 & fact>0.){
    newClCut_pine[,c(2,4)] <- ClCut_pine[,c(2,4)]*(1+fact)
    newClCut_spruce[,c(2,4)] <- ClCut_spruce[,c(2,4)]*(1+fact)
    newClCut_birch[,c(2,4)] <- ClCut_birch[,c(2,4)]*(1+fact)
  }else{
    newClCut_pine[,c(2,4)] <- ClCut_pine[,c(2,4)]+fact
    newClCut_spruce[,c(2,4)] <- ClCut_spruce[,c(2,4)]+fact
    newClCut_birch[,c(2,4)] <- ClCut_birch[,c(2,4)]+fact
  } 
  return(list(ClCut_pine=newClCut_pine,
              ClCut_spruce=newClCut_spruce,
              ClCut_birch=newClCut_birch))
}

updatePclcut <- function(initPrebas,pClCut){
  nSites <- initPrebas$nSites
  ClCut <- initPrebas$ClCut
  inDclct <- initPrebas$inDclct
  ETSmean <- rowMeans(initPrebas$ETSy)
  ETSthres <- 1000
  climIDs <- initPrebas$siteInfo[,2]
  siteType <- initPrebas$siteInfo[,3]
  inDclct <- initPrebas$inDclct
  inAclct <- initPrebas$inAclct
  for(i in 1: nSites){
    if(ClCut[i]==1) inDclct[i,] <-
        c(ClCutD_Pine(ETSmean[climIDs[i]],ETSthres,siteType[i],pClcut= pClCut$ClCut_pine),
          ClCutD_Spruce(ETSmean[climIDs[i]],ETSthres,siteType[i],pClcut= pClCut$ClCut_spruce),
          ClCutD_Birch(ETSmean[climIDs[i]],ETSthres,siteType[i],pClcut= pClCut$ClCut_birch),
          0,0,0,0,0,0,0)  ###"fasy","pipi","eugl","rops","popu",'eugrur','piab(DE)')
    if(ClCut[i]==1) inAclct[i,] <-
        c(ClCutA_Pine(ETSmean[climIDs[i]],ETSthres,siteType[i],pClcut= pClCut$ClCut_pine),
          ClCutA_Spruce(ETSmean[climIDs[i]],ETSthres,siteType[i],pClcut= pClCut$ClCut_spruce),
          ClCutA_Birch(ETSmean[climIDs[i]],ETSthres,siteType[i],pClcut= pClCut$ClCut_birch),
          80,50,13,30,50,13,120)  ###"fasy","pipi","eugl","rops","popu",'eugrur','piab(DE)')
  }
  return(list(inDclct=inDclct,inAclct=inAclct))
}

#returns a the dominant species or the age of dominant species for each site at each year
###varX="species" -> returns the dominant species
###varX="age" -> returns the age of dominant layer
domFun <- function(modOut,varX="species"){
  nSites <- modOut$nSites
  nYears <- modOut$maxYears
  segID <- modOut$siteInfo[,1]
  
  oo <- as.vector(apply(modOut$multiOut[,,30,1:3,1],1:2,which.max))  
  oo <- cbind(rep(1:nSites,nYears),
              rep(1:nYears,each=nSites),
              oo)
  if(varX=="species") domX <- matrix(modOut$multiOut[,,4,1:3,1][oo],
                                     nSites,nYears)
  if(varX=="age") domX <- matrix(modOut$multiOut[,,7,1:3,1][oo],
                                 nSites,nYears)
  outX <- data.table(segID=segID,domX)
}


###retunrs the Volume of deciduous
##modOut -> multiPREBAS output
vDecFun <- function(modOut){
  segID <- modOut$siteInfo[,1]
  oo <- data.table(which(modOut$multiOut[,,4,,1]==3,arr.ind=T))
  setnames(oo,c("site","year","layer"))
  vx <-modOut$multiOut[,,30,,1][as.matrix(oo)]
  oo$Vdec <- vx
  setkey(oo,site,year)
  ff <- oo[,sum(Vdec),by=.(site,year)]
  VdecMat <- matrix(0,modOut$nSites,modOut$maxYears)
  VdecMat[as.matrix(ff[,1:2])] <- unlist(ff[,3])
  outX <- data.table(segID=segID,VdecMat)
}


#####extract model output as baweighted mean or sum according to funX
##modOut -> multiPREBAS output
##varSel -> variable ID 
##funX -> function to summarize the output accross layers (sum or BA weighted mean (baWmean))
outProcFun <- function(modOut,varSel,funX="baWmean"){
  segID <- modOut$siteInfo[,1]
  marginX <- 1:2
  if(funX=="baWmean"){
    outX <- data.table(segID=segID,baWmean(modOut,varSel))
  }
  if(funX=="sum"){
    outX <- data.table(segID=segID,apply(modOut$multiOut[,,varSel,,1],marginX,sum))
  }
  setnames(outX,c("segID",1:modOut$maxYears))
  return(outX)
}

SBBbivoltinePotential <- function(initPrebas=initPrebas,nYears){
  # climid, year, date, vars c(PAR, TAir, VPD, Precip, CO2)
  #initPrebas <- sampleXs0$initPrebas$weather
  SBB_funct <- "lange" # "lange" or "seidl"

  nT <- nrow(initPrebas$weather)
  SBBgen <- matrix(0,nT,nYears)
  
  for(yi in 1:nYears){
    wi <- 2 # Tair
    Ti <- initPrebas$weather[,yi,,wi]
    stageNames <- c("flight","egg","larvae","pupae","adult")
    Talpha <- c(5, 10.6, 8.2, 9.9, 3.2) # Temp.threshold values for phases
    nStages <- length(Talpha) # Number of development stages
    pulpaeStages <- 1:3*nStages-2 # stage from which pulpae development begins: If the next stage is reached, adult generation is generated
    Falpha <- c(110, 51.8, 204.4, 57.7, 238.5) # Degree sum limits
    Talpha <- rep(Talpha,3)
    Falpha <- rep(Falpha,3)
    
    stageID <- matrix(0,nT,ncol(Ti))
    #Fid <- matrix(0,nT,nYears)
    Titmp <-Ti
    for(ij in 1:nrow(Titmp)){
      Titmp[ij,1:which(Ti[1,]>(19.5-5))[1]]<-0
    }
    
    if(SBB_funct=="lange"){
      for(alpha in 1:length(Talpha)){
        # print(alpha)
        Falphai <- Falpha[alpha]
        Di <- t(apply((Titmp-Talpha[alpha])*(Titmp>=Talpha[alpha]),1,cumsum))
        for(ij in 1:nrow(Di)){
          stageFinish <- which(Di[ij,]/Falphai>=1)
          stageID[ij,Di[ij,]>0] <- alpha -1 + min(1,Di[ij,ncol(Di)]/Falphai)
          if(length(stageFinish)>0){  
            stageFinish <- stageFinish[1]
            #print(stageFinish)
            #  Fid[ij,(which(0<Di[ij,])[1]):stageFinish] <- alpha
            Titmp[ij,1:stageFinish] <- 0
          } else { 
            Titmp[ij,] <- 0
          }
        }
        #  if(alpha==11) break()
        #  print(rbind(#Tii<-Ti[ij,],Titmpi=Titmp[ij,],
        #    Dii=Di[ij,],stagei=stageID[ij,],Fi=Fid[ij,])[,70:300])
      }
      for(alpha in pulpaeStages){
        ss<-matrix(1,nrow(Di),2)
        ss[,2]<-stageID[,ncol(stageID)]-alpha
        SBBgen[stageID[,ncol(stageID)]>alpha,yi]<-SBBgen[stageID[,ncol(stageID)]>alpha,yi]+apply(ss,1,min)[stageID[,ncol(stageID)]>alpha]
      }
    } else if(SBB_funct=="seidl"){
      SBBgen[,yi]<-apply((Titmp-5)*(Titmp>8.3),1,sum)/557
      SBBgen[(SBBgen%%1<0.6)]<-floor(SBBgen[(SBBgen%%1<0.6)])
    }
  }
  #if(SBB_funct=="lange") SBBgen<-SBBgen-1

  #print(SBBgen)
  return(SBBgen)
  
}

SBB_predisposition <- function(modOut){
  
  ba <- data.table(apply(modOut$multiOut[,,13,,1],1:2,sum))
  age <-data.table(apply(modOut$multiOut[,,7,,1]*(modOut$multiOut[,,4,,1]==2),1:2,mean)) # should it be the age of spruce or the whole stand?
  ba_spruce <-data.table(apply(modOut$multiOut[,,13,,1]*(modOut$multiOut[,,4,,1]==2),1:2,sum))
  share_spruce <- ba_spruce/ba
  share_spruce[share_spruce=="NaN"]<-0
  ba <- ba_spruce
  aSW<-apply(modOut$multiOut[,,40,,1],1:2,sum)
  aSW[aSW>1]<-1
  aSW <- 1 - aSW
  
  aSW_lims <- c(0,0.05, 0.1, 0.15, 0.2, 0.3, 1.01)
  aSW_facts <- c(0.0, 0.05,0.4,0.6,0.9,0.95,1)
  PIdrought <- matrix(aSW_facts[findInterval(as.matrix(aSW),aSW_lims)],nrow = dim(modOut$multiOut)[1],ncol=nYears)
  
  share_lims <- c(0,0.1, 0.25, 0.50, 0.7, 1.01)
  share_facts <- c(0.08, 0.17, 0.5, 0.83, 1.00)
  PIspruce <- matrix(share_facts[findInterval(as.matrix(share_spruce),share_lims)],nrow = dim(modOut$multiOut)[1],ncol=nYears)
  
  age_lims <- c(0,60,80,100,1e4)
  age_facts <- c(0.2, 0.6, 0.9, 1.0)
  PIage <- matrix(age_facts[findInterval(as.matrix(age),age_lims)],nrow = dim(modOut$multiOut)[1],ncol=nYears)
  
  ba_lims <- c(0, 10, 20, 30, 40, 60, 1e3)
  ba_facts <- c(0.9, 0.5, 0.3, 0.2, 0.3, 0.4)
  PIba <- matrix(ba_facts[findInterval(as.matrix(ba),ba_lims)],nrow = dim(modOut$multiOut)[1],ncol=nYears)
  
  #PI = 0.3*PIspruce + 0.25*PIage + 0.15*PIba + 0.3*PIdrought
  PI = (0.3*PIspruce + 0.25*PIage + 0.15*PIba)/0.7*PIdrought
  #PI = PIspruce*PIage*PIba*PIdrought
  return(PI)
}

SBB_damage_prob <- function(PI,SBBbp,clim_ids){
  # calculate probability for SBB damage from PI and SBB bivoltine potential

  SBBgen <- SBBbp
  # Seidl et al. 2007:
  GEN <- 0*SBBgen
  GEN[SBBgen==0]<-0
  GEN[SBBgen==1]<-0.1
  GEN[SBBgen>=0.6 & SBBgen<=1]<-0.2
  GEN[SBBgen==2]<-0.6
  GEN[SBBgen>2]<-1
  x1 <- -1.51
  x2 <- 1.65
  pBByr <- 1 - exp(x1*PI^x2)^GEN[clim_ids,]
  return(pBByr)
}



validationPeriodEstimates <- function(sampleXs0){ # Calculate yearly statistics for validation period
  valPeriod <- 2016:2023
  if(dev.interactive()) dev.off()
  figdim <- c(floor(sqrt(length(valPeriod)))+1,ceiling(sqrt(length(valPeriod))))
  
  # Vspruce
  colsi <- c(1:length(valPeriod)+1)
  outXi <- vSpFun(sampleXs0$region,SpID = 2)[,..colsi]
  outXi <- array(unlist(outXi), c(dim(outXi)))
  par(mfrow = figdim)
  drawRaster(outXi=outXi, varNams = paste0("Vspruce_",valPeriod),r_nos = r_nos_stations[[station_id]])
  
  # pSBB
  par(mfrow = figdim)
  outXi <- sampleXs0$region$multiOut[,valPeriod-2015,"Rh/SBBpob[layer_1]",1,2]
  Ep_SBBarea <- colSums(ops[[1]]$area*outXi)
  drawRaster(outXi=outXi, varNams = paste0("pSBB_",valPeriod),r_nos = r_nos_stations[[station_id]])
  
  # SBB intensity
  par(mfrow = figdim)
  outXi <- sampleXs0$region$multiOut[,valPeriod-2015,48,1,2]
  E_SBBIntensity <- colSums(ops[[1]]$area*outXi)
  drawRaster(outXi=outXi, varNams = paste0("SBBintensity_",valPeriod),r_nos = r_nos_stations[[station_id]])
  
  # BA disturbed
  par(mfrow = figdim)
  outXi <- sampleXs0$region$multiOut[,valPeriod-2015,"grossGrowth/bb BA disturbed",1,2]
  E_BAdisturbed <- colSums(ops[[1]]$area*outXi)
  drawRaster(outXi=outXi, varNams = paste0("BAdisturbed_",valPeriod),r_nos = r_nos_stations[[station_id]])
  
  # 1-SMI
  par(mfrow = figdim)
  outXi <- 1-sampleXs0$region$multiOut[,valPeriod-2015,"NEP/SMI[layer_1]",1,2]
  Ep_SBBBAdisturbed <- colSums(ops[[1]]$area*outXi)
  drawRaster(outXi=outXi, varNams = paste0("(1-SMI)_",valPeriod),r_nos = r_nos_stations[[station_id]])
  
  lopeta <- tahan
}


drawRaster <- function(outXi, varNams, r_nos){
  data.IDs_rnos <- data.frame()
  for(r_noi in r_nos){  
    load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_noi,"_IDsTab.rdata"))
    #data.IDs$segID <- data.IDs$maakuntaID
    data.IDs_rnos <- rbind(data.IDs_rnos, data.IDs)
  }
  data.IDs<-data.IDs_rnos
  rm("data.IDs_rnos")
  gc()
  
  data.IDs <- data.IDs[segID!=0]
  setkey(data.IDs,segID)
  
  outX <- cbind(outXi,ops[[1]])
  colnames(outX)[1:ncol(outXi)] <- varNams
  setkey(outX,segID)
  setkey(data.IDs,segID)
  
  tabX <- merge(data.IDs,outX)
  colnames(tabX)[colnames(tabX)=="x.x"]<-"x"
  colnames(tabX)[colnames(tabX)=="y.x"]<-"y"
  #dev.off()
  ndat <- sample(1:nrow(tabX),1000)
  #plot(tabX$x[ndat],tabX$y[ndat],col="blue")
  
  tabX$x <- tabX$x - stations[station_id,"x_UTM"]
  tabX$y <- tabX$y - stations[station_id,"y_UTM"]
  #tabX <- merge(outX,data.IDs)
  #colnames(tabX)[colnames(tabX)=="x.y"]<-"x"
  #colnames(tabX)[colnames(tabX)=="y.y"]<-"y"
  rm(outX);gc()
  
  #filee <- paste0("../rasters/",stations$name[station_id],"_",varX,"_harscen",harvScen,"_harInten",harvInten,"_",rcps,"_Nswitch",restrictionSwitch)
  
  xi <- which(colnames(tabX) == "x")
  yi <- which(colnames(tabX) == "y")
  
  h <- hist(outXi, 5, plot=F)
  library(RColorBrewer)
  colorsi <- cols <- brewer.pal(length(h$breaks),"YlOrRd")
  #colorsi <- c(terrain.colors(length(h$breaks)))
  if(min(outXi)<1e-5) colorsi <- c("gray",colorsi[-1])
  
  legendi <- F
  for(ij in 1:length(varNams)){
    j <- which(colnames(tabX) == varNams[ij])
    colsi <- c(xi,yi,j)
    rastX <- rasterFromXYZ(tabX[,..colsi])
    crs(rastX) <- crsX
    if(ij == length(varNams)) legendi <- T
    plot(rastX, breaks = h$breaks, 
         col = colorsi, 
         #xlab="m", ylab="m",
         main = varNams[ij], yaxt="n", xaxt="n",
         legend = legendi, axes = FALSE, box = F)
    #box(col = "white")
    #writeRaster(rastX,filename = paste0(filee,"_per1.tiff"),overwrite=T)
  }  
  rm(tabX);
  rm(rastX);gc()
}

