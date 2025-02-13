rm(list=ls())
gc()
if(dev.interactive()) dev.off()
outType <- "testRun"
neighborIDs <- T
toFile <- T
asParallel <- F
if(neighborIDs) asParallel <- T
if(!exists("toFile")) toFile <- T
nSegs <- 200
if(toFile) nSegs <- 10000
ttAll <- T # T = Add data outside forest declarations to the simulations
outputPlot <- T # T = if use existing data for plotting and function generation, F if sample construction
#if(neighborIDs) outputPlot <- F

set.seed(1)
toMem <- ls()
#nYears <- 2024-2015

workdir <-"/scratch/project_2000994/PREBASruns/PREBAStesting/"
setwd(workdir)
library(sf)
library(data.table)
library("readxl")
library("raster")
library("terra")

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex =  cex.cor * (1 + r) / 2)
}
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
}
panel.lm <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = "black", ...) {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  abline(stats::lm(y ~ x),  col = col.smooth, ...)
}

rnames <- c("Uusimaa","Ahvenanmaa","Keski-Pohjanmaa","Pirkanmaa","Etel%C3%A4-Karjala","Keski-Suomi",
            "Pohjois-Savo","Lappi_E","Lappi_P","Kanta-H%C3%A4me","Pohjanmaa","Varsinais-Suomi",
            "Etel%C3%A4-Pohjanmaa","P%C3%A4ij%C3%A4t-H%C3%A4me","Satakunta","Kymenlaakso",
            "Kainuu","Etel%C3%A4-Savo","Pohjois-Karjala","Pohjois-Pohjanmaa")
regnames <- c("Uusimaa","Ahvenanmaa","Keski-Pohjanmaa","Pirkanmaa","Etela-Karjala","Keski-Suomi",
              "Pohjois-Savo","Lappi_E","Lappi_P","Kanta-Hame","Pohjanmaa","Varsinais-Suomi",
              "Etela-Pohjanmaa","Paijat-Hame","Satakunta","Kymenlaakso",
              "Kainuu","Etela-Savo","Pohjois-Karjala","Pohjois-Pohjanmaa")
CSCrun <- T
#vPREBAS <- "newVersion" #
vPREBAS <- "master"   

savepath = "/scratch/project_2000994/PREBASruns/PREBAStesting/MKIdata/"

ContinueIterations <- F

#r_noi <- r_no
rnos <- c(1:8,8:19)
rids <- rids0 <- c(1,3:length(rnos))
#rids0 <- c(20,19,8,17,7)
#rids0 <- c(6,18,4,9,13)
#rids0 <- c(12,15,5,10,11)
#rids0 <- c(1,14,16,3)
if(!outputPlot){
  rids <- rids0
}
#regnames[rids]

if(!toFile) rids <- c(1,5)
#if(toFile) pdf(paste0(savepath,"regionalCharacteristics.pdf"))
#rids <- c(1,5)
r_noi <- rids[1]

nYears <- 2023-2015
outputs <- data.frame()
rids0 <- 1
if(ContinueIterations){
  load("MKIdata/SBB_sample_output.rdata")
  rids0 <- max(outputs$reg)
}
r_noi <- rids[rids0]
dam_names <- c("alldeclarations","all","SBB","wind","fire","moose")
dam_indexs <- c("all","0","1602","1504","1503","1650")

inds <- 1

if(outputPlot){
  ################################################################################
  TestaaSBBkoodi <- T
  fmi_from_allas <- T
  Continue_outputs <- T
  
  NEW <- T
  if(!NEW){# Load region specific data
    load(file="../PREBAStesting/MKIdata/samples.rdata")
    samples <- ops
    rm(list="ops")
    gc()
  }
  if(Continue_outputs){
    load(paste0(savepath,"SBB_sample_output.rdata"))
    ij0 <- which(max(outputs$reg)==rids)+1
  } else {
    outputs <- list()
    ij0 <- 1
  }
  ij <- ij0
  r_noi <- rids[ij]
  
  for(ij in ij0:length(rids)){#rids[rids0:length(rids)]){ #######################################
    setwd(workdir)
    time0 <- Sys.time()
    toMem <- ls()
    r_no <- rnos[rids[ij]]
    r_noi <- rids[ij]
    if(NEW){
      load(paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/MKIdata/samples_",r_noi,".rdata"))
      ops <- list(dataS)
      rm("dataS")
      gc()  
    } else {
      ops <- list(samples[[which(as.numeric(names(samples))==r_noi)]])
      #ops <- list(samples[[which(as.numeric(names(samples))==r_noi)]][1:1000,])
    }
    #ops[[1]] <- ops[[1]][1:10,]
    
    print(paste("Region",r_no))
    print(paste("SBB area",sum(ops[[1]]$area[which(ops[[1]]$forestdamagequalifier=="1602")])))
    #source("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Settings_IBCCarbon.R", local=T)  
    
    print(fmi_from_allas)
    # fmi data from allas
    if(fmi_from_allas){
      toMemFmi <- ls()
      source("0.5_get_fmi_from_allas.R")
      repo <- "ForModLabUHel/fmi.weather.finland"
      file_path <- "r/init_setup.R"
      branch <- "main"
      # Get init functions from github
      init_funs <- fetch_file_from_github(repo, file_path, branch)
      eval(parse(text = init_funs))
      rm(init_funs, file_path, repo)
      # SET PARAMETERS
      resolution <- 1 # Resolution in km (1, 5 or 9)
      years <- c(2015:2024) # For which years to extract (1961:2023 are full years)
      save_path <- paste0(getwd()) # Where to save the extracted data.table as .rdata
      repo_url <- "https://github.com/ForModLabUHel/fmi.weather.finland.git" # Project repository to use
      format_to_prebas <- T # TRUE for Prebas format, FALSE for raw data. Default is TRUE.
      
      req_coords_dt <- data.table(
        id = 1:nrow(ops[[1]]),
        E = ops[[1]]$x,
        N = ops[[1]]$y
      )
      req_coords <- as.matrix(req_coords_dt[, c("E", "N")]) # The coords are passed as a matrix
      
      # Set parameters
      params <- list(req_coords = req_coords, resolution = resolution, years = years)
      
      # Combine arguments
      setup_and_run_args <- c(params, list(save_path = save_path, repo_url = repo_url, format_to_prebas = format_to_prebas))
      
      # RUN
      result <- do.call(setup_and_run, setup_and_run_args)
      
      # Change file name
      file.rename(list.files(path=workdir, pattern="fmi_vars_", all.files=FALSE,full.names=FALSE)[1],
                  "fmi_vars_PREBAS.rdata")
      file.rename(list.files(path=workdir, pattern="climID_lookup_", all.files=FALSE,full.names=FALSE)[1],
                  "climID_lookup.rdata")
      rm(list = setdiff(ls(), toMemFmi))
      gc()
    }
    
    # PREBAS runs
    source("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Settings_IBCCarbon.R", local=T)  
    area_tot <<- sum(data.all$area)
    #rm(list="data.all")
    #gc()
    climatepath_orig = "/scratch/project_2000994/RCP/"
    station_id <- "tmp"
    climScen <- 0
    nSegs <- nrow(ops[[1]])
    print(paste("Spruce = 0 segments", length(which(ops[[1]]$spruce==0 & ops[[1]]$ba==0))))
    print(paste("SBB area",sum(ops[[1]]$area[which(ops[[1]]$forestdamagequalifier=="1602")])))
    #print(head(ops[[1]]))
    
    #    source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
    #source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
    source("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/functions_IBSCarbon.R", local=T)
    source("~/adaptFirst_runs/functions.R", local=T)
    #source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
    mortMod <<- 1
    nYears <<- 2050-2015
    endingYear <<- nYears + startingYear
    byManual <- T
    if(byManual){
      sampleID<-1 
      disturbanceON = c("bb","fire","wind")
      uncRCP=0
      rcps = "CurrClim_fmi"
      #rcps = "CurrClim"
      harvScen<-"Base"
      harvInten<-"Base"
      easyInit=FALSE
      forceSaveInitSoil=F 
      cons10run = F
      procDrPeat=F
      ingrowth <- F
      coeffPeat1=-240
      coeffPeat2=70
      coefCH4 = 0.34#g m-2 y-1
      coefN20_1 = 0.23
      coefN20_2 = 0.077#g m-2 y-1
      landClassUnman=NULL
      compHarvX = 0
      funPreb = regionPrebas
      initSoilCreStart=NULL
      outModReStart=NULL
      reStartYear=1
      sampleX=NULL
      deltaID <- 1
      P0currclim=NA 
      fT0=NA
      TminTmax <- NA
      #clcuts <<- 1
    }
    disturbanceON <- c("fire","wind","bb")
    mortMod <<- 1
    #ops <<- list(ops[[1]][1:(nSegs/2),],ops[[1]][((nSegs/2)+1):nSegs,])
    nSegs2 <- dim(ops[[1]])[1]
    
    #out <- runModelAdapt(1,sampleID=1, outType = outType, rcps = "CurrClim", harvScen="Base", harvInten="Base", forceSaveInitSoil=T)
    #out_fmi <- runModelAdapt(1,sampleID=1, outType = outType, rcps = "CurrClim_fmi", harvScen="Base", harvInten="Base",forceSaveInitSoil=T)
    
    sampleIDs <- 1:length(ops)
    lapply(sampleIDs, 
           function(jx) { 
             runModelAdapt(1,sampleID=jx, outType = outType, rcps = "CurrClim_fmi",
                           harvScen="Base", harvInten="Base",
                           forceSaveInitSoil=T)
             gc()
           })
    print(paste("SBB area",sum(ops[[1]]$area[which(ops[[1]]$forestdamagequalifier=="1602")])))
    #for(sampleID in sampleIDs){
    #  sampleXs0 <- runModelAdapt(1,sampleID=sampleID, outType = outType, rcps = "CurrClim_fmi",
    #                             harvScen="Base", harvInten="Base",#ingrowth = F, clcut = clcuts,
    #                             forceSaveInitSoil=T)#, disturbanceON = disturbanceON)
    #}
    byManual <- T
    if(byManual){
      sampleID<-1 
      #disturbanceON = "bb"
      uncRCP=0
      rcps = "CurrClim_fmi"
      harvScen<-"NoHarv"
      harvInten<-"NoHarv"
      easyInit=FALSE
      forceSaveInitSoil=F 
      cons10run = F
      procDrPeat=F
      coeffPeat1=-240
      ingrowth <- F
      coeffPeat2=70
      coefCH4 = 0.34#g m-2 y-1
      coefN20_1 = 0.23
      coefN20_2 = 0.077#g m-2 y-1
      landClassUnman=NULL
      compHarvX = 0
      funPreb = regionPrebas
      initSoilCreStart=NULL
      outModReStart=NULL
      reStartYear=1
      sampleX=NULL
      deltaID <- 1
      P0currclim=NA 
      fT0=NA
    }
    
    nYears <<- 2024-2015
    endingYear <<- nYears + startingYear
    clcuts <<- 1
    sampleXs0 <- lapply(sampleIDs, 
                        function(jx) { 
                          runModelAdapt(1,sampleID = jx, outType=outType,  
                                        harvScen="NoHarv",rcps = "CurrClim_fmi", # clcut = clcuts,
                                        harvInten="NoHarv")
                          #ingrowth = T, 
                          #clcut = -1, disturbanceON = disturbanceON)
                        })
    
    print(paste("SBB area",sum(ops[[1]]$area[which(ops[[1]]$forestdamagequalifier=="1602")])))
    print("SMIs:")
    print(sampleXs0[[1]]$region$multiOut[1,,"NEP/SMI[layer_1]",1,2])
    print("Range SMIs:")
    print(range(sampleXs0[[1]]$region$multiOut[,,"NEP/SMI[layer_1]",1,2]))
    print("")
    
    multiout <- array(0,
                      c(nSegs,dim(sampleXs0[[1]]$region$multiOut)[-1]),
                      dimnames = dimnames(sampleXs0[[1]]$region$multiOut))
    for(iii in 1:dim(multiout)[3]){
      for(iji in 1:dim(multiout)[4]){
        for(kki in 1:dim(multiout)[5]){
          for(zzi in 1:length(sampleIDs)){
            multiout[(((zzi-1)*nSegs2)+1):(nSegs2*zzi),,iii,iji,kki] <- 
              sampleXs0[[zzi]]$region$multiOut[,,iii,iji,kki]
          }
        }
      }
    }
    #plaah <- otapwst
    print(multiout[6,,"BA",,1])
    any(multiout[,,"W_wsap/fireRisk[layer_1]",1,2]>0)
    file.remove(paste0(workdir,"fmi_vars_PREBAS.rdata"))
    file.remove(paste0(workdir,"climID_lookup.rdata"))
    
    #################################################
    #bb risk is in 45 (first layer,status=2)(1,2)
    pSBB <- multiout[,,45,1,2]
    #bb potential impact is in 48 (2) (share of basal area damaged) (first layer,status=2)(1,2)
    damageSBB <- multiout[,,48,1,2]
    #bb disturbed basal area is in 43
    damBASBB <- multiout[,,43,1,2]
    ################################################
    damagedYearlyAreaSim <- colSums(damageSBB*pSBB*ops[[1]]$area)
    names(damagedYearlyAreaSim) <- (startingYear+1):endingYear 
    #ops <- list(samples[[which(as.numeric(names(samples))==r_noi)]])
    
    damagedYearlyAreaObs <- array(0,c(1,length(unique(ops[[1]]$dam_year))))
    damYears <- sort(unique(ops[[1]]$dam_year))
    for(tii in 1:length(damYears)){
      damagedYearlyAreaObs[tii] <- sum(ops[[1]]$area[
        which(ops[[1]]$forestdamagequalifier=="1602" & 
                ops[[1]]$dam_year==damYears[tii])])
    }
    dimnames(damagedYearlyAreaObs)[[2]] <- damYears
    print("Observed bb area")
    print(0.1*damagedYearlyAreaObs)
    print("simulated bb area")
    print(damagedYearlyAreaSim[which((startingYear+1):endingYear%in%unique(ops[[1]]$dam_year))])
    gc()
    source("/scratch/project_2000994/PREBASruns/PREBAStesting/dayl_calc.R", local = T)
    lat <- array(0,c(nSegs,1))
    for(zzi in 1:length(sampleIDs)){
      lat[(((zzi-1)*nSegs2)+1):(nSegs2*zzi)] <- 
        sampleXs0[[zzi]]$region$latitude
    }
    #lat <- sampleXs0$region$latitude#location$y
    
    daylights <- dayl(lat) # climIDs x 365 matrix of daylight length
    daylights <- daylights[[2]]
    daylights[is.na(daylights)] <- 0
    #################################################
    weatherIDs <- 0
    siteInfo <- array(0,c(nSegs,dim(sampleXs0[[1]]$region$siteInfo)[2]), 
                      dimnames = dimnames(sampleXs0[[1]]$region$siteInfo))
    for(zzi in 1:length(sampleIDs)){
      weatherIDs <- weatherIDs + dim(sampleXs0[[zzi]]$region$weather)[1]
      siteInfo[((zzi-1)*nSegs2+1):(zzi*nSegs2),] <- sampleXs0[[zzi]]$region$siteInfo
    }
    weather <- array(0,
                     c(weatherIDs,dim(sampleXs0[[zzi]]$region$weather)[-1]),
                     dimnames = dimnames(sampleXs0[[1]]$region$weather))
    for(iii in 1:dim(weather)[3]){
      for(iji in 1:dim(weather)[4]){
        weatherID <- 0
        for(zzi in 1:length(sampleIDs)){
          weather[(weatherID+1):(weatherID + dim(sampleXs0[[zzi]]$region$weather)[1]),,iii,iji] <- 
            sampleXs0[[zzi]]$region$weather[,,iii,iji]
          weatherID <- weatherID + dim(sampleXs0[[zzi]]$region$weather)[1]
        }
      }
    }
    
    
    tmpsum <- weather[,,,2] # daily temperature for each year
    tmpsum <- tmpsum[siteInfo[,"climID"],,]
    # From publication:
    Tpar <- data.table(alpha=0.02876507, beta=3.5922336,
                       gamma=1.24657367, Tmax=40.9958913, T0=30.4, DTL=8.3)
    Teff <- function(temp) {
      as.integer(temp > Tpar$DTL & temp < 38.9)*(Tpar$T0-Tpar$DTL)*
        (exp(Tpar$alpha*temp)- exp(Tpar$alpha*Tpar$Tmax-(Tpar$Tmax-temp)/Tpar$beta)-Tpar$gamma)
    }
    effBTS <- apply(tmpsum, 1:2, Teff)
    #effBTS <- effBTS[,sampleXs0$region$siteInfo[,"climID"],]
    daylights <- array(t(daylights),dim(effBTS))
    daylights[daylights<14.5] <- 0    
    daylights[daylights>=14.5] <- 1
    #daylights <- matrix(as.integer(daylights >= 14.5),365,dim(effBTS)[2],byrow = F)
    #daylights <- array(daylights, dim(effBTS))
    effBTS <- effBTS*daylights
    effBTS <- apply(effBTS,2:3,sum)
    TsumSBB <- effBTS/557 #[sampleXs0$region$siteInfo[,"climID"],]/557
    rm(list = c("effBTS","daylights"))
    gc()
    
    if(TestaaSBBkoodi){
      ###test daylight & TsumSBB function
      si <-190
      ti <- 8
      TsumSBBtest <- rep(0,3)
      daylight <- rep(0,365)
      temp <- weather[siteInfo[,"climID"][si],ti,,2]
      #temp <- tmpsum[si,ti,]
      xx <- .Fortran("TsumSBBfun", lat=as.double(lat[si]),
                     temp=as.double(temp),TsumSBB=as.double(TsumSBBtest[1]))
      
      print(xx$TsumSBB)
      print(TsumSBB[si,ti])
      
    }
    #lopetatahan <- lopeta
    
    #############################################################
    
    sampleXs <- multiout
    
    timeCol <- as.numeric(ops[[1]]$dam_year)-2015
    apick <- function(a,timeCol){
      nI <- nrow(a)
      y <- array(0,c(nI,1))
      for(nn in 1:nI){
        y[nn] <- a[nn,timeCol[nn]]
      }
      return(y)
    }
    
    nYears <- sampleXs0[[1]]$initPrebas$maxYears
    varsSBB <- c("age","D","BA","V","N")
    output <- data.frame(matrix(0,nSegs,length(varsSBB)))
    # spruce specific vaiables
    ind <- 1
    SpID <- 2
    varsi <- varsSBB[1]
    for(varsi in varsSBB){
      oo <- data.table(which(sampleXs[,,"species",,1]==SpID,arr.ind=T))
      setnames(oo,c("site","year","layer"))
      a <- sampleXs[,,varsi,,1][as.matrix(oo)]
      oo$VSp <- a
      setkey(oo,site,year)
      if(varsi%in%(c("BA","V","N"))){    
        ff <- oo[,sum(VSp),by=.(site,year)]
      } else {
        ff <- oo[,max(VSp),by=.(site,year)]
      }
      a <- matrix(0,nSegs,nYears)
      a[as.matrix(ff[,1:2])] <- unlist(ff[,3])
      assign(paste0(varsi,"s_spruce"), a)
      output[,ind] <- data.table(apick(a,timeCol))#a[col(a)==timeCol])
      #    if(varsi=="aSW") output[which(output[,ind]>1),ind] <- NA 
      ind <- ind+1
    }
    varsSBB <- paste0(varsSBB,"_spruce")
    colnames(output) <- paste0(varsSBB)
    # Total N
    a<-apply(sampleXs[,,"N",,1],1:2,sum)
    assign("Nsums", a)
    Nsum <- apick(a,timeCol)
    output <- cbind(output, Nsum=Nsum)#data.table("BAsum" = a[col(a)==timeCol]))
    # Total V
    a<-apply(sampleXs[,,"V",,1],1:2,sum)
    assign("Vsum", a)
    Vsum <- apick(a,timeCol)
    output <- cbind(output, Vsum=Vsum)#data.table("BAsum" = a[col(a)==timeCol]))
    # Total BA
    a<-apply(sampleXs[,,"BA",,1],1:2,sum)
    assign("BAsums", a)
    BAsum <- apick(a,timeCol)
    output <- cbind(output, BAsum=BAsum)#data.table("BAsum" = a[col(a)==timeCol]))
    # spruce BA fraction
    a <- BAs_spruce/BAsums
    a[is.na(a)] <- 0
    assign("BAspruceFracts", a)
    BAspruceFract <- apick(a, timeCol)
    output <- cbind(output, BAspruceFract = BAspruceFract)#output[,"BA_spruce"]/output[,"BAsum"]))
    print(summary(output))
    #    # Add layers for temperature sum
    # Add layer for Tsum for SBB at Tprevprev3
    TsumSBBs <- TsumSBB
    x = data.table(apick(TsumSBB, timeCol-3))
    colnames(x) <- "TsumSBBprev3"
    output <- cbind(output,x)
    # Add layer for Tsum for SBB at Tprevprev3
    x = data.table(apick(TsumSBB, timeCol-2))
    colnames(x) <- "TsumSBBprev2"
    output <- cbind(output,x)
    # Add layer for Tsum for SBB at Tprevprev
    x = data.table(apick(TsumSBB, timeCol-1))
    colnames(x) <- "TsumSBBprev"
    output <- cbind(output,x)
    # Add layer for Tsum for SBB at Tpred
    x = data.table(apick(TsumSBB, timeCol))
    colnames(x) <- "TsumSBB"
    output <- cbind(output,x)
    # sitetype
    output <- cbind(output,data.table(sitetype = sampleXs[,1,"sitetype",1,1]))
    # colnames to memory
    varsSBB <- colnames(output)
    
    # also previous year in weather based data
    varsSBBw <- c("SMI","aSW","ETS","ET_preles")
    ind <- ncol(output)
    for(varsi in varsSBBw){
      # SMI
      if(varsi=="SMI"){ 
        a <- sampleXs[,,46,1,2]
      } else {
        a <- sampleXs[,,varsi,1,1]
      }  
      assign(paste0(varsi,"s"), a)
      
      a3 = data.table(apick(a, timeCol-3))
      colnames(a3) <- paste0(varsi,"_Tprev3")
      a2 = data.table(apick(a, timeCol-2))
      colnames(a2) <- paste0(varsi,"_Tprev2")
      a1 = data.table(apick(a, timeCol-1))
      colnames(a1) <- paste0(varsi,"_Tprev")
      a0 = data.table(apick(a, timeCol))
      colnames(a0) <- paste0(varsi,"_T")
      
      if(varsi=="aSW"){ 
        a0[which(a0>1)] <- NA 
        a1[which(a1>1)] <- NA 
        a2[which(a2>1)] <- NA 
        a3[which(a3>1)] <- NA 
      }
      output <- cbind(output,a3,a2,a1,a0)
      ind <- ncol(output)
    }
    
    nvarsSBB <- ncol(output)
    namesSBB <- colnames(output)
    
    #########################################################
    
    if(TestaaSBBkoodi){
      standInfo <- multiout[si,ti,c(4,7,13),,1]
      nLayers <- ncol(standInfo)
      spruceIDs <- c(2,10)
      nSpIDs <- length(spruceIDs)
      spruceStandVars <- rep(0,3)
      rBAspruce <- 0
      
      spruceVarsTest <- .Fortran("spruceVars",
                                 standInfo = as.matrix(standInfo),
                                 nLayers = as.integer(nLayers),
                                 spruceIDs = as.integer(spruceIDs),
                                 nSpIDs=as.integer(nSpIDs),
                                 spruceStandVars = as.double(spruceStandVars),
                                 rBAspruce = as.double(rBAspruce))
      
      print(spruceVarsTest$spruceStandVars)
      print(c(BAs_spruce[si,ti],ages_spruce[si,ti],BAsums[si,ti]))
    }
    
    
    
    # PI for BA spruceFract
    #x0 = 0.6
    #k = -10.
    #PI_spruceFract = 1./(1.+exp(k* (BAspruceFract - x0)))
    x <- BAspruceFracts 
    x0 <- 0.4
    k <- -10
    fromFortran <- T
    if(fromFortran) {   
      x0 <- 0.6 # if from fortran
      k <- -10
    }
    fPI <- function(x){
      y <- 0
      xi <-which(!is.na(x))
      y[xi] <- 1/(1+exp(k* (x[xi] - x0)))
    }
    PI_spruceFract <- apply(x,1:2,fPI)
    
    #x0 = 85.
    #k = -0.05
    #PI_agespruce = 1.0/(1.+exp(k* (age_spruce - x0)))
    
    x <- ages_spruce 
    x0 <- 60
    k <- -0.1
    fPI <- function(x){
      y <- 0
      xi <-which(!is.na(x))
      y[xi] <- 0.2 + 0.8/(1+exp(k* (x[xi] - x0)))
    }
    if(fromFortran) {   
      x0 = 85.
      k = -0.05
      fPI <- function(x){
        y <- 0
        xi <-which(!is.na(x))
        y[xi] <- 1/(1+exp(k* (x[xi] - x0)))
      }
    }
    PI_agespruce <- apply(x,1:2,fPI)
    
    x <- BAs_spruce
    x0 <- 1
    k <- -1
    #    x0 <- 20
    #    k <- -.1
    fPI <- function(x){
      y <- 0
      xi <-which(!is.na(x))
      y[xi] <- 1/(1+exp(k* (x[xi] - x0)))
    }
    if(fromFortran) {   
      x0 = 17.
      k = -0.250
      fPI <- function(x){
        y <- 0
        xi <-which(!is.na(x))
        y[xi] <- 1/(1+exp(k* (x[xi] - x0)))
      }
    }
    
    PI_BAspruce <- apply(x,1:2,fPI)
    
    
    x <- SMIs
    x0 <- 0.89
    k <- 100
    #    x0 <- 0.92
    #    k <- 100
    fPI <- function(x){
      y <- 0
      xi <-which(!is.na(x))
      y[xi] <- 1/(1+exp(k* (x[xi] - x0)))
    }
    
    if(fromFortran) {   
      x0 = 0.955
      k = 400.
      fPI <- function(x){
        y <- 0
        xi <-which(!is.na(x))
        y[xi] <- 1/(1+exp(k* (x[xi] - x0)))
      }
    }
    PI_SMI <- apply(x,1:2,fPI)
    
    TsumsMean <- c(1.418036, 1.632891, 1.408584)
    C <- matrix(c(10.16699, -22.40190, -13.48438,
                  0.00000,  20.23612, -17.89729,
                  0.00000,   0.00000,  31.11091),3,3,byrow = T)
    TsumsMean <- c(1.418036, 1.632891)
    C <- matrix(c(10.16699, -22.40190,
                  0.00000,  20.23612),2,2,byrow = T)
    
    TsumsT <- array(0,c(nrow(TsumSBBs),nYears))    
    #ti <- 1
    TsumSBBss <- cbind(TsumSBBs[,1],TsumSBBs[,1],TsumSBBs)
    for(tii in 1:(nYears)){
      Tsums <- cbind(TsumSBBss[,tii], TsumSBBss[,tii+1])#, TsumSBBss[,ti+2])
      dX <- Tsums-matrix(TsumsMean,nrow=nrow(Tsums),
                         ncol=ncol(Tsums),byrow = T)
      Tsums <- exp(-rowSums((dX%*%C)^2)/4)
      Tsums[Tsums<0.0001]<-0
      Tsums[Tsums>=0.0001]<-1
      #      Tsums[which(rowSums(matrix(as.integer(dX > 0),nrow(dX),3,byrow = F))==3)] <- 1
      Tsums[which(rowSums(matrix(as.integer(dX > 0),nrow(dX),2,byrow = F))==2)] <- 1
      
      TsumsT[,tii] <- Tsums
    }
    
    aspruceshare <- 0.3
    aage <- 0.25
    aBA <- 0.15
    adrought <- max(1-aspruceshare-aage-aBA,0)
    
    #    PI_SMITprev <- cbind(PI_SMI[,1],PI_SMI[,-nYears])
    PI_SMITx <- PI_SMI[,c(1,1:(ncol(PI_SMI)-1))]
    PI <- (aspruceshare*PI_spruceFract+
             aage*PI_agespruce+
             aBA*PI_BAspruce+
             adrought*PI_SMITx)
    
    # GEN is the bark beetle generation index, which depends on temperature
    Tsum <- TsumSBBss[,c(-1,-nYears)] # The previous year TsumSBB is used for bark beetle generations
    if(fromFortran) Tsum <- TsumSBBs[,c(1,1:(ncol(TsumSBBs)-1))]
    gen <- array(0,dim(Tsum))  
    gen[Tsum<1] <-0
    gen[Tsum>=1 & Tsum<1.5] <- 0.1
    gen[Tsum>=1.5 & Tsum<2] <- 0.2
    gen[Tsum>=2 & Tsum<2.5] <- 0.6
    gen[Tsum>=2.5] <- 1
    f0 <- 1
    if(fromFortran) {
      TsumSBBss <- array(0,dim(TsumSBBs))
      for(tii in 1:ncol(TsumSBBs)){
        TsumSBBss[,tii] <- rowSums(TsumSBBs[,c(max(1,tii-3),max(1,tii-2),max(1,tii-1),tii)])/4
      }
      x0 = 1.49
      k = -10. 
      fPI <- function(x){
        y <- 0
        xi <-which(!is.na(x))
        y[xi] <- 1/(1+exp(k* (x[xi] - x0)))
      }      
      f0 <- apply(TsumSBBss,1:2,fPI)
      
    }
    
    # probability function coefficients from Seild et al. 2007
    x1 <- -1.51
    x2 <- 1.65
    
    # SBB probability
    #    prob <-(1-exp(TsumsT*x1*PI^x2)^gen)
    prob <-(1-exp(x1*PI^x2)^gen)*f0
    # (1.0d0-exp(x1*PI**x2)**gen) * f0
    
    
    # indices for SBB and no SBB
    nSBB <- which(ops[[1]]$forestdamagequalifier=="1602")
    noSBB <- setdiff(1:nSegs,nSBB)
    
    ############################# 
    if(TestaaSBBkoodi){
      pBBx=rep(0.,5)
      TsumSBBsi <- c(TsumSBB[si,ti-3],TsumSBB[si,ti-2],TsumSBB[si,ti-1],TsumSBB[si,ti])
      BA_spruceX <- BAs_spruce[si,ti]
      BAtotX <- BAsums[si,ti]
      age_spruceX <- ages_spruce[si,ti]
      SMIX <- SMIs[si,ti-1]
      # (1.0d0-exp(x1*PI**x2)**gen) * f0,  PI_spruceFract, PI_agespruce, PI_BAspruce, PI_SMITprev
      ff <- .Fortran("riskBB",pBB=as.double(pBBx),
                     TsumSBBs=as.double(TsumSBBsi),
                     BA_spruce=as.double(BA_spruceX),
                     BAtot=as.double(BAtotX),
                     age_spruce=as.double(age_spruceX),
                     SMI=as.double(SMIX)) 
      
      print(paste("FF:",ff$pBB))
      print(paste("VJ:",c(prob[si,ti],PI_spruceFract[si,ti], PI_agespruce[si,ti], PI_BAspruce[si,ti], PI_SMI[si,ti-1])))
      
      print(paste("Prebas",pSBB[si,ti]))
      
      nSBB <- which(ops[[1]]$forestdamagequalifier=="1602")
      noSBB <- setdiff(1:nSegs,nSBB)
      if(length(nSBB)>0){
        print("mean prob, SBB segments:")
        print(colMeans(prob[nSBB,]))
      }
      print("mean prob, noSBB segments:")
      print(colMeans(prob[noSBB,]))
      if(length(nSBB)>0){
        print("FM mean prob, SBB segments:")
        print(colMeans(pSBB[nSBB,]))
      }
      print("FM mean prob, noSBB segments:")
      print(colMeans(pSBB[noSBB,]))
    }
    ns <- sample(1:nSegs,min(nSegs,2000))
    par(mfrow=c(1,1))
    
    #########################################################
    
    # Add ops[[1]]-data and declaration data to the outputs
    #if(ttAll){
    output <- cbind(ops[[1]][,c("segID","N","area","fert","minpeat","landclass",
                                "cons","forestdamagequalifier","dam_id","dam_year",
                                "damSegId","cuttingpurpose",
                                "cuttingrealizationpractice")],
                    output)
    #} else {
    #  output <- cbind(ops[[1]][,c("segID","N","area","fert","minpeat","landclass","cons")],
    ##                  rbind(tt[ni,])[,c("forestdamagequalifier","dam_id","damSegID",
    #                                    "cuttingpurpose","cuttingrealizationpractice")],
    #                  output)
    
    # neighboring data
    if(neighborIDs){
      # add info about neighboring SBB damages  
      colsNeigh <- c(which(names(ops[[1]])=="minDist"):which(names(ops[[1]])=="clearcutYear_Wind"))
      output <- cbind(output,ops[[1]][,..colsNeigh])
      nvarsSBB <- nvarsSBB+length(colsNeigh)
    }
    
    output_mem <- output#[,-2]
    
    output2 <-output[output[,.I[which.max(BA_spruce)], by=dam_id]$V1,] # one row for each damage polygon
    #rrows <- which(output$damSegID %in% output2$damSegID)
    #output3 <- output2[rrows[match(output$damSegID,output2$damSegID)],]
    #  output$damNeighbor[is.na(damNeighbor)]=1
    varsSBB <- colnames(output)[(ncol(output)-nvarsSBB+1):ncol(output)]
    if(length(nSBB)>0){
      # as pixels
      output <- output2
      output <- output[rep(1:nrow(output), as.vector(output$N)),]  
      plotPairs <- F
      if(plotPairs){
        ns <- sample(1:nrow(output),5000)
        pairs(output[ns,..varsSBB],na.rm=T,
              upper.panel = panel.cor,
              diag.panel  = panel.hist,
              lower.panel = panel.smooth,
              pch = ".",
              main="all pixels")
        pairs(output[output$forestdamagequalifier=="1602",..varsSBB],na.rm=T,
              upper.panel = panel.cor,
              diag.panel  = panel.hist,
              lower.panel = panel.smooth,
              pch = ".",
              main="SBB pixels")
      }
      ncols <- ncol(output)
      for(ij in (ncols-length(varsSBB)+1):ncols){
        if(length(which(!is.na(output[,..ij])))>0){
          par(mfrow = c(3,1))
          a <- hist(as.numeric(as.matrix(output[,..ij])),plot=F)
          barplot(a$counts,names.arg = a$mids,xlab=colnames(output)[ij],
                  main = "all sim data",
                  ylab="ha sim")
          b <- hist(as.numeric(as.matrix(output[output$forestdamagequalifier==1602 & 
                                                  output$BA_spruce>0,..ij])),
                    breaks=a$breaks,plot=F)
          barplot(b$counts,names.arg = a$mids,xlab=colnames(output)[ij],
                  main = "SBB sim",
                  breaks=a$breaks,
                  ylab="ha sim")
          barplot(b$counts/a$counts*100,names.arg = a$mids,
                  main=paste(regnames[r_noi],"/","SBB %"),
                  xlab=colnames(output)[ij],ylab="% of declarations")
        }
      }
    }
    #    outputs <- rbind(outputs, cbind(data.table(reg = r_noi),output_mem))
    output_mem <- cbind(data.table(reg = r_noi),output_mem)
    #return(output_mem)
    outputs <- rbind(outputs,output_mem)
    if(toFile){ 
      save(outputs,file=paste0(savepath,"SBB_sample_output.rdata"))
      print(paste0("outputs saved for region ",r_no,"/",regnames[r_noi]))
    }  
    print(Sys.time()-time0)
    rm(list=setdiff(ls(),c("outputs",toMem)))
    gc()
  }  
}
if(FALSE){
  rnos <- c(1:8,8:19)
  rids <- c(1,3:length(rnos))
  load("damInfo_5rdata")
  regnames <- c("Uusimaa","Ahvenanmaa","Keski-Pohjanmaa","Pirkanmaa","Etela-Karjala","Keski-Suomi",
                "Pohjois-Savo","Lappi_E","Lappi_P","Kanta-Hame","Pohjanmaa","Varsinais-Suomi",
                "Etela-Pohjanmaa","Paijat-Hame","Satakunta","Kymenlaakso",
                "Kainuu","Etela-Savo","Pohjois-Karjala","Pohjois-Pohjanmaa")
  
  declInformation <- array(0,c(length(rids),dim(damInfo)[2]+1,dim(damInfo)[1]),
                           dimnames = list(regnames[rids],c("rno",colnames(damInfo)),rownames(damInfo)))
  ii <- 1
  for(ii in 1:length(rids)){
    r_noi <- rids[ii]
    load(paste0("damInfo_",r_noi,"rdata"))
    for(ij in 1:nrow(damInfo)){
      declInformation[ii,,ij] <- c(r_noi,damInfo[ij,])
    }
  }
  dimnames(declInformation)[[3]][7:8] <- dimnames(declInformation)[[3]][c(8,7)]
  
  #declInformation <- declInformation[,,c(1:6,8,7)]
  decltypes <- dimnames(declInformation)[[3]]
  
  if(dev.interactive()==T) dev.off()
  
  pdf(file="declInformation.pdf")
  par(mar=c(11,4,8,4))
  par(mfrow=c(1,1))
  for(ij in 1:dim(declInformation)[3]){
    #for(ti in 1:3){
    title <- paste(decltypes[ij])
    #if(ti==3){    
    if(ij%%2==0){
      tmp <- declInformation[,-1,ij]/declInformation[,-1,ij-1]*100
      tmp[which(is.na(tmp))] <- 0
      barplot(height = t(tmp), beside=T, las=2,
              names.arg = row.names(declInformation), 
              main=title, ylab="%")
    } else {
      barplot(height = t(declInformation[,-1,ij]), beside=T, las=2,
              names.arg = row.names(declInformation), ylab="ha", 
              main=title)
    }
  }
  dev.off()
}

