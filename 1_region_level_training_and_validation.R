#rm(list=ls())
#gc()
if(dev.interactive()) dev.off()
if(!exists("fmi_from_allas")) fmi_from_allas <- F
if(!exists("weighted")) weighted <- F
if(!exists("sbatches")) sbatches <- F
if(!exists("onlyValidationset")) onlyValidationset <- T
if(!exists("newSamples")) newSamples <- T # T if generate new samples, F if only run PREBAS with previous samples
outType <- "testRun"
if(!exists("neighborIDs")) neighborIDs <- F
toFile <- T
asParallel <- F
if(neighborIDs) asParallel <- T
if(!exists("toFile")) toFile <- T
if(!exists("nSegs")) nSegs <- 30000
#if(toFile) nSegs <- 30000
ttAll <- T # T = Add data outside forest declarations to the simulations

set.seed(1)
toMem <- ls()

workdir <-"/scratch/project_2000994/PREBASruns/PREBAStesting/"
setwd(workdir)
library(sf)
library(data.table)
library("readxl")
library("raster")
library("terra")

rnames <- c("Uusimaa","Ahvenanmaa","Keski-Pohjanmaa","Pirkanmaa","Etel%C3%A4-Karjala","Keski-Suomi",
            "Pohjois-Savo","Lappi_E","Lappi_P","Kanta-H%C3%A4me","Pohjanmaa","Varsinais-Suomi",
            "Etel%C3%A4-Pohjanmaa","P%C3%A4ij%C3%A4t-H%C3%A4me","Satakunta","Kymenlaakso",
            "Kainuu","Etel%C3%A4-Savo","Pohjois-Karjala","Pohjois-Pohjanmaa")
regnames <- c("Uusimaa","Ahvenanmaa","Keski-Pohjanmaa","Pirkanmaa","Etela-Karjala","Keski-Suomi",
              "Pohjois-Savo","Lappi_E","Lappi_P","Kanta-Hame","Pohjanmaa","Varsinais-Suomi",
              "Etela-Pohjanmaa","Paijat-Hame","Satakunta","Kymenlaakso",
              "Kainuu","Etela-Savo","Pohjois-Karjala","Pohjois-Pohjanmaa")
CSCrun <- T

if(!exists("vPREBAS")) vPREBAS <- "newVersion" #
#vPREBAS <- "master"   

savepath = "/scratch/project_2000994/PREBASruns/PREBAStesting/MKIdata/"

ContinueIterations <- F

#r_noi <- r_no
rnos <- c(1:8,8:19)
rids <- rids0 <- c(1,3:length(rnos))
#rids <- rids0 <- rids[10:length(rnos)]
if(!exists("setX")) setX<-1
if(!exists("ridsi")) ridsi <- 1:length(rids)
rids <- rids0 <- rids[ridsi]
#if(setX==1){
#  rids <- rids0 <- rids[5:7]
#}  else if(setX==2){
#  rids <- rids0 <- rids[16:length(rnos)]
#}
#rids0 <- c(20,19,8,17,7)
#rids0 <- c(6,18,4,9,13)
#rids0 <- c(12,15,5,10,11)
#rids0 <- c(1,14,16,3)
#regnames[rids]

#if(!toFile) rids <- c(1,5)
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
calculateOPSdata  <-  function(r_noi, nSegs=1000, neighborIDs=T, weighted = T, climScen=0){
  toMem <- ls()
  r_no <- rnos[r_noi]
  
  print(paste("region",r_no))
  print(paste("Neighbor information =",neighborIDs))
  
  # PREBAS run for a sample
  landClassX <- 1:2
  #    devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/settings.r")
  source("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Settings_IBCCarbon.R", local=TRUE)  
  vars_to_prebas <- colnames(data.all)
  areatot <- sum(data.all$area)
  
  if(climScen==0 | r_no %in% c(8,9)){
    
    fname <- paste0("DeclaredDamages_",dam_names[inds],"_rno",r_no,"_",regnames[r_noi],".rdata")
    load(file=paste0(savepath,"/",fname))
    print(paste("File",fname,"opened"))
    
    XYdamages <- XYdamages[dam_year>2018 & dam_year<2024,]
    XYdamages <- XYdamages[which(!is.na(pine)&!is.na(spruce)&!is.na(birch)),] # remove data where there was NAs
    gc()
    
    IDsUniq <- unique(XYdamages[,c("dam_id","segID")])
    IDsUniq[,damSegID:=1:nrow(IDsUniq)]
    
    # Give segments that are split with declaration polygons, a unique name
    XYdamages[,damSegID := match(paste(XYdamages$dam_id,XYdamages$segID),
                                 paste(IDsUniq$dam_id,IDsUniq$segID))]
    areas <- XYdamages[, .N, by = list(damSegID)]
  }
  if(climScen==0){
    XYdam_uniqueSegm <- XYdamages[XYdamages[,.I[which.max(spruce)], by=damSegID]$V1] # one row for each segment
    XYdam_uniqueSegm[,N := areas[match(areas$damSegID, XYdam_uniqueSegm$damSegID),"N"]]
    XYdam_uniqueSegm[,area := N*16^2/100^2]
    
    #if(weighted) tt <- tt[tt$dam_year>=2019,]
    XYdam_uniqueSegm <- XYdam_uniqueSegm[segID%in%data.all$segID,] # some damage polygons are outside the landclasses 1:2
    gc()
    
    if(TRUE){
      ntmp <- match(XYdam_uniqueSegm$segID,data.all$segID)
      ntmp2 <- which(!(colnames(XYdam_uniqueSegm)%in%colnames(data.all)))
      ntmp2 <- cbind(data.all[ntmp,], XYdam_uniqueSegm[,..ntmp2])
      ntmp2[,c("N","x","y","area")] <- XYdam_uniqueSegm[,c("N","x","y","area")]
      XYdam_uniqueSegm <- ntmp2
      rm(list=c("ntmp","ntmp2"))
      gc()
    } else  {   # columns that are not in XYdam tables
      not_in_XYdam_uniqueSegm <- colnames(data.all)[which(!(colnames(data.all)%in%
                                                              colnames(XYdam_uniqueSegm)))]
      # find data for these columns from data.all
      XYdam_uniqueSegm <- cbind(XYdam_uniqueSegm,
                                data.all[match(XYdam_uniqueSegm$segID,data.all$segID),..not_in_XYdam_uniqueSegm])
    }
    ######################### clearcuts?
    # types of cuttings classified as clearcut
    cuttinginpractise <- c(5, 8, 16, 17, 19, 21, 22, 24)
    yrs <- sort(as.numeric(unique(XYdam_uniqueSegm$dam_year)))
    
    # function for yearly damaged areas
    areasDam <- function(x,id) sum(XYdam_uniqueSegm$area[
      which(XYdam_uniqueSegm$forestdamagequalifier==id & XYdam_uniqueSegm$dam_year==x)])
    # function for yearly damaged areas that are clearcut
    areasDamClct <- function(x,id) sum(XYdam_uniqueSegm$area[
      which(XYdam_uniqueSegm$forestdamagequalifier==id & XYdam_uniqueSegm$dam_year==x  & 
              XYdam_uniqueSegm$cuttingrealizationpractice%in%cuttinginpractise)])
    
    # Table for saving the damage infor
    damInfo <- array(0,c(8,length(yrs)),
                     dimnames = list(c("alldeclarea","alldeclclctarea","SBBarea","SBBclctarea","windarea",
                                       "windclctarea","firearea","fireclctarea"),paste0("Year",yrs)))
    for(y in 1:length(yrs)) damInfo[1,y] <- sum(XYdam_uniqueSegm$area[which(XYdam_uniqueSegm$dam_year==yrs[y])])
    for(y in 1:length(yrs)) damInfo[2,y] <- sum(XYdam_uniqueSegm$area[which(XYdam_uniqueSegm$dam_year==yrs[y] & XYdam_uniqueSegm$cuttingrealizationpractice%in%cuttinginpractise)])
    for(y in 1:length(yrs)) damInfo[3,y] <- areasDam(yrs[y],dam_indexs[which(dam_names=="SBB")])
    for(y in 1:length(yrs)) damInfo[4,y] <- areasDamClct(yrs[y],dam_indexs[which(dam_names=="SBB")])
    for(y in 1:length(yrs)) damInfo[5,y] <- areasDam(yrs[y],dam_indexs[which(dam_names=="wind")])
    for(y in 1:length(yrs)) damInfo[6,y] <- areasDamClct(yrs[y],dam_indexs[which(dam_names=="wind")])
    for(y in 1:length(yrs)) damInfo[7,y] <- areasDam(yrs[y],dam_indexs[which(dam_names=="fire")])
    for(y in 1:length(yrs)) damInfo[8,y] <- areasDamClct(yrs[y],dam_indexs[which(dam_names=="fire")])
    #if(weighted) save(damInfo,file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/damInfo_",r_noi,".rdata"))
  }
  ########################
  #nt <- 1:nrow(XYdam_uniqueSegm)
  load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
  data.IDs$segID <- data.IDs$maakuntaID
  data.IDs <- data.IDs[segID!=0]
  setkey(data.IDs,segID)
  setkey(data.all,segID)
  tabX <- merge(data.IDs,data.all) # coords of the segments in sample outside declarations
  x <- tabX$x[match(data.all$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"x"]
  y <- tabX$y[match(data.all$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"y"]
  rm("tabX"); gc()
  data.all[,x:=x]
  data.all[,y:=y]

  gc()
  KUVA <- F
  if(KUVA){
    setkey(data.all,segID)
    tabX <- merge(data.IDs,data.all) # coords of the segments in sample outside declarations
    x <- tabX$x[match(data.all$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"x"]
    y <- tabX$y[match(data.all$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"y"]
    ni <- sample(1:length(x),1000, replace = F)
    plot(x[ni],y[ni])
  }
  if(r_no %in% c(8,9)){ # Lappi E and P: divide data.all according to y-coordinate
    setkey(data.all,segID)
    tabX <- merge(data.IDs,data.all) # coords of the segments in sample outside declarations
    x <- tabX$x[match(data.all$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"x"]
    y <- tabX$y[match(data.all$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"y"]
    rm("tabX"); gc()
    xybb <- array(0,c(nrow(XYdamages),2))  
    xybb[,1] <- XYdamages$x
    xybb[,2] <- XYdamages$y
    xx <- bbox(xybb)
    rm("xybb"); gc()
    if(r_no == 8){ 
      yLappiEP <- xx[2,2] 
      data.all <- data.all[y<=yLappiEP,]
      gc()
    }
    if(r_no == 9) { 
      yLappiEP <- xx[2,1] 
      data.all <- data.all[y>=yLappiEP,]
      gc()
    }
  }
  
  ###### combine data.all and XYdam_uniqueSegm
  if(climScen==0){
    # Which columns are in 
    XYdam_uniqueSegm_names_in <- which(colnames(XYdam_uniqueSegm)%in%vars_to_prebas)
    dataAll_names_in <- which(vars_to_prebas%in%colnames(XYdam_uniqueSegm))
    
    # rows in data.all, which are in XYdam_uniqueSegm data
    ni <- which((data.all$segID %in% XYdam_uniqueSegm$segID))
    
    # for the XYdam_uniqueSegm data, add data.all information -> All variables already included!
    #tmp <- data.all[ni[match(XYdam_uniqueSegm$segID, data.all$segID[ni])],]
    # Sort columns of XYdam_uniqueSegm -> first data.all-columns, then decl information
    nicols <- which(!colnames(XYdam_uniqueSegm)%in%vars_to_prebas)
    XYcols <- colnames(XYdam_uniqueSegm)[nicols]
    XYdam_uniqueSegm <- XYdam_uniqueSegm[,c(..vars_to_prebas,..XYcols)]
    # get the sorter column ids  
    nicols <- which(!colnames(XYdam_uniqueSegm)%in%vars_to_prebas)
    XYcols <- colnames(XYdam_uniqueSegm)[nicols]
    #XYdam_uniqueSegm <- cbind(tmp,XYdam_uniqueSegm[,..nicols])
    
    # rows in data.all, which are not in XYdam_uniqueSegm data
    ni <- which(!(data.all$segID %in% XYdam_uniqueSegm$segID))
    data.all <- cbind(data.all[ni,],array(0,c(length(ni), length(XYcols))))
    colnames(data.all)[nicols] <- colnames(XYdam_uniqueSegm)[nicols]
    data.all$dam_year <- sample(yrs,nrow(data.all),replace = T) # give random year for monitoring
    data.all$dam_id <- data.all$damSegID <- paste0("0",1:nrow(data.all)) # damageid starting with 0
    
    # Combine the two sets
    data.all <- rbind(XYdam_uniqueSegm, data.all)
  }
  
  # validation set as a completelly random set
  ni <- sample(1:nrow(data.all), nSegs, replace=F)
  sampleValidation <- data.all[ni,]
  
  if(climScen==0){
    # training set from the rest of segments: 
    data.all <- data.all[setdiff(1:nrow(data.all),ni),]
    gc()
    
    weightedTraining <- weighted
    print(paste("weighted Training =",weighted))
    if(weightedTraining){
      # damaged segments, max nSegs/3 rows
      damID <- dam_indexs[dam_names=="SBB"]
      ni <- which(data.all$forestdamagequalifier==damID)
      nsDam <- round(nSegs/3)
      if(length(ni) > nsDam) ni <- ni[sample(1:length(ni), nsDam, replace = F)]
      sampleTraining <- data.all[ni,]
      # Sample from outside damaged segments
      ns <- nSegs - length(ni)
      ni <- which(!(data.all$forestdamagequalifier%in%c(damID)))
      ni <- ni[sample(1:length(ni),ns, replace = F)]
      sampleTraining <- rbind(sampleTraining, data.all[ni,])
      rm("data.all")
      gc()
    } else {
      print("No weighting on training set.")
      ni <- sample(1:nrow(data.all), nSegs, replace=F)
      sampleTraining <- data.all[ni,]
    }
    # X and y coordinates for the samples
    if(FALSE){
      setkey(data.IDs,segID)
      xycols <- which(colnames(sampleTraining)%in%c("x","y"))
      #nicols <- setdiff(1:ncol(sampleTraining),which(colnames(sampleTraining)%in%c("x","y")))
      #sampleTraining <- sampleTraining[,..nicols]
      ni <- which(sampleTraining[,c("x","y")]==0, arr.ind=T)[,1]
      tmp <- sampleTraining[ni,]
      setkey(tmp,segID)
      tabX <- merge(data.IDs,tmp) # coords of the segments in sample outside declarations
      x <- tabX$x[match(tmp$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"x"]
      y <- tabX$y[match(tmp$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"y"]
      #setkey(sampleTraining,segID)
      #tabX <- merge(data.IDs,sampleTraining) # coords of the segments in sample outside declarations
      #x <- tabX$x[match(sampleTraining$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"x"]
      #y <- tabX$y[match(sampleTraining$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"y"]
      sampleTraining[ni,x:=x]
      sampleTraining[ni,y:=y]
      rm(list=c("x","y","tabX")); gc()
      if(KUVA){
        ni <- sample(1:nSegs,1000,replace = F)
        points(sampleTraining$x[ni],sampleTraining$y[ni], col="blue")
      }
      if(length(which(is.na(sampleTraining$x)))>0){
        nNas <- which(is.na(sampleTraining$x))
        print(sampleTraining[nNas,])
      }
    }
    if(FALSE){
      setkey(data.IDs,segID)
      nicols <- setdiff(1:ncol(sampleValidation),which(colnames(sampleValidation)%in%c("x","y")))
      sampleValidation <- sampleValidation[,..nicols]
      setkey(sampleValidation,segID)
      tabX <- merge(data.IDs,sampleValidation) # coords of the segments in sample outside declarations
      x <- tabX$x[match(sampleValidation$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"x"]
      y <- tabX$y[match(sampleValidation$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"y"]
      sampleValidation[,x:=x]
      sampleValidation[,y:=y]
      rm(list=c("x","y","tabX")); gc()
      if(KUVA){
        points(sampleValidation$x[ni],sampleValidation$y[ni],col="red")
      }
    }  
  }
  #####
  # For the sample, indexes about neighbors
  print(paste("Calculate neighbor information =",neighborIDs))
  if(neighborIDs){# & weighted){
    iks <- 1:2
    if(climScen>0) iks <- 2
    ik <- iks[1]
    for(ik in iks){
      if(ik==1) {dataS <- sampleTraining
      print("neighbors for training set")}
      if(ik==2){ dataS <- sampleValidation
      print("neighbors for validation set")}
      # all declarations, all coordinates and years of interest
      ntmp <- which(XYdamages$cuttingrealizationpractice%in%cuttinginpractise | 
                      XYdamages$forestdamagequalifier==dam_indexs[dam_names=="SBB"] |
                      XYdamages$forestdamagequalifier==dam_indexs[dam_names=="wind"])
      ntmp <- ntmp[order(XYdamages$y[ntmp],XYdamages$x[ntmp],decreasing = T)]
      yyd<-XYdamages$y[ntmp]
      xxd<-XYdamages$x[ntmp]
      dam_idd <- XYdamages$dam_id[ntmp]
      dam_yeard <- as.numeric(XYdamages$dam_year[ntmp])
      dam_indd <- XYdamages$forestdamagequalifier[ntmp]
      dam_crpd <- XYdamages$cuttingrealizationpractice[ntmp] # cutting realisation practise
      
      # closest distance and clearcut year for each segment in sample 
      # closest distance and clearcut year for each segment in sample 
      # closest distance and clearcut year for each segment in sample 
      # closest south distance and clearcut yearfor each segment in sample 
      clearcutNeighbor_SBB <- clearcutNeighbor_wind <- clearcutNeighbor <- clearcutNeighbor_south <- matrix(NA,nSegs,2) 
      
      id_ij <- 1
      ij <- 1
      print("Start finding neighbor information")
      timeT <- Sys.time()
      neighborSets <- sampleSets <- 0
      
      #setkey(dataS,segID)
      #tabX <- merge(data.IDs,dataS) 
      dimNams <- c("minDist","clearcutYear", "minDist_SBB", "clearcutYear_SBB", "minDist_Wind", "clearcutYear_Wind", "minDist_s", "clearcutYear_s")
      #  outputNeighbor <- array(0,c(8,nSegs),dimnames = list(dimNams,c(1:nSegs)))
      timeT <- Sys.time()
      damCCutInt <- as.integer(dam_crpd%in%cuttinginpractise)
      damCCutInt[damCCutInt==0] <- 1e12
      print(paste("CCut obs. n =",length(which(damCCutInt==1))))
      damBBInt <- as.integer(dam_indd==dam_indexs[dam_names=="SBB"])
      damBBInt[is.na(damBBInt)] <- 1e12
      damBBInt[damBBInt==0] <-1e12
      print(paste("SBB obs. n =",length(which(damBBInt==1))))
      damWInt <- as.integer(dam_indd==dam_indexs[dam_names=="wind"])
      damWInt[is.na(damWInt)] <- 1e12
      damWInt[damWInt==0] <-1e12
      print(paste("Wind obs. n =",length(which(damWInt==1))))
      #source("../PREBAStesting/0.5_functions.R", local = T)
      source("~/disturbance_modules/0.5_functions_updated.R", local=T)
      outputNeighbor <- apply(data.table(c(1:nSegs)),1,neighborsAll,dataS=dataS,
                              damCCutInt=damCCutInt, damBBInt=damBBInt, damWInt=damWInt)
      outputNeighbor <- data.table(t(outputNeighbor))
      colnames(outputNeighbor) <- dimNams
      outputNeighbor[outputNeighbor==1e12] <- NA
      print(outputNeighbor)
      #rm(list=c("xx","yy","dx"))
      gc()
      if(ik==1){
        sampleTraining <- cbind(dataS,outputNeighbor)  
        save(sampleTraining,file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/MKIdata/sample_train_",r_noi,".rdata"))
      } 
      if(ik==2){
        sampleValidation <- cbind(dataS,outputNeighbor)
        save(sampleValidation,file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/MKIdata/sample_valid_",r_noi,".rdata"))
      } 
    }
  }
  if(!toFile & climScen==0) head(sampleTraining)
  if(!toFile) head(sampleValidation)
  #samples <- list(sampleTraining, sampleValidation)
  #names(samples) <- c("sampleTraining", "sampleValidation")
  if(climScen==0) {
    samples <- list(sampleTraining, sampleValidation, damInfo,areatot, vars_to_prebas)
    names(samples) <- c("sampleTraining", "sampleValidation","damInfo","areatot","vars_to_prebas")
  }
  if(climScen>0){
    samples <- list(sampleValidation, areatot, vars_to_prebas)
    names(samples) <- c("sampleValidation","areatot","vars_to_prebas")
  } 
  if(weighted){
    save(samples, damInfo, 
         file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/MKIdata/samples_train_valid_",r_noi,".rdata"))
    print("sample data saved.")
  }
  return(samples)
  rm(list=setdiff(ls(),toMem))
  gc()
  
}

###############################################################################
TestaaSBBkoodi=F
trainingSetCreation <- function(r_noi, sampleXs, dataS, startingYear=2015, endingYear=2050,
                                TestaaSBBkoodi=F, neighborIDs=F,setname="training"){ # Inputs to modelfitting
  r_no <- rnos[r_noi]
  print(paste("Calculate training set for region",r_no))
  multiout <- sampleXs$region$multiOut
  #print(multiout[6,,"BA",,1])
  any(multiout[,,"W_wsap/fireRisk[layer_1]",1,2]>0)
  outputs <- data.frame

  #################################################
  #bb risk is in 45 (first layer,status=2)(1,2)
  pSBB <- multiout[,,45,1,2]
  #bb potential impact is in 48 (2) (share of basal area damaged) (first layer,status=2)(1,2)
  damageSBB <- multiout[,,48,1,2]
  #bb disturbed basal area is in 43
  damBASBB <- multiout[,,43,1,2]
  ################################################
  damagedYearlyAreaSim <- colSums(damageSBB*pSBB*dataS$area)
  names(damagedYearlyAreaSim) <- 2015+(1:length(damagedYearlyAreaSim))#(startingYear+1):endingYear 
  #ops <- list(samples[[which(as.numeric(names(samples))==r_noi)]])
  
  damagedYearlyAreaObs <- array(0,c(1,length(unique(dataS$dam_year))))
  damYears <- sort(unique(dataS$dam_year))
  for(tii in 1:length(damYears)){
    damagedYearlyAreaObs[tii] <- sum(dataS$area[
      which(dataS$forestdamagequalifier=="1602" & 
              dataS$dam_year==damYears[tii])])
  }
  dimnames(damagedYearlyAreaObs)[[2]] <- damYears
  print("Observed bb area")
  print(0.1*damagedYearlyAreaObs)
  print("simulated bb area")
  print(damagedYearlyAreaSim[which((startingYear+1):endingYear%in%unique(dataS$dam_year))])
  gc()
  source("/scratch/project_2000994/PREBASruns/PREBAStesting/dayl_calc.R", local = T)
  lat <- sampleXs$region$latitude
  print("latitudes")
  print(lat[1:10])
  daylights <- dayl(lat) # climIDs x 365 matrix of daylight length
  daylights <- daylights[[2]]
  daylights[is.na(daylights)] <- 0
  #################################################
  weatherIDs <- 0
  siteInfo <- sampleXs$region$siteInfo
  weather <- sampleXs$region$weather

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
  daylights <- array(t(daylights),dim(effBTS))
  daylights[daylights<14.5] <- 0    
  daylights[daylights>=14.5] <- 1
  effBTS <- effBTS*daylights
  effBTS <- apply(effBTS,2:3,sum)
  TsumSBB <- effBTS/557
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

  #############################################################
  
  #sampleXs <- multiout
  
  timeCol <- as.numeric(dataS$dam_year)-2015
  apick <- function(a,timeCol){
    nI <- nrow(a)
    y <- array(0,c(nI,1))
    for(nn in 1:nI){
      y[nn] <- a[nn,timeCol[nn]]
    }
    return(y)
  }
  
  nYears <- sampleXs$initPrebas$maxYears
  varsSBB <- c("age","D","BA","V","N")
  output <- data.frame(matrix(0,nSegs,length(varsSBB)))
  # spruce specific vaiables
  ind <- 1
  SpID <- 2
  varsi <- varsSBB[1]
  for(varsi in varsSBB){
    oo <- data.table(which(multiout[,,"species",,1]==SpID,arr.ind=T))
    setnames(oo,c("site","year","layer"))
    a <- multiout[,,varsi,,1][as.matrix(oo)]
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
  a<-apply(multiout[,,"N",,1],1:2,sum)
  assign("Nsums", a)
  Nsum <- apick(a,timeCol)
  output <- cbind(output, Nsum=Nsum)#data.table("BAsum" = a[col(a)==timeCol]))
  # Total V
  a<-apply(multiout[,,"V",,1],1:2,sum)
  assign("Vsum", a)
  Vsum <- apick(a,timeCol)
  output <- cbind(output, Vsum=Vsum)#data.table("BAsum" = a[col(a)==timeCol]))
  # Total BA
  a<-apply(multiout[,,"BA",,1],1:2,sum)
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
  output <- cbind(output,data.table(sitetype = multiout[,1,"sitetype",1,1]))
  # colnames to memory
  varsSBB <- colnames(output)
  
  # also previous year in weather based data
  varsSBBw <- c("SMI","aSW","ETS","ET_preles")
  ind <- ncol(output)
  for(varsi in varsSBBw){
    # SMI
    if(varsi=="SMI"){ 
      a <- multiout[,,46,1,2]
    } else {
      a <- multiout[,,varsi,1,1]
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
  

  TsumsT <- array(0,c(nrow(TsumSBBs),nYears))    
  #ti <- 1
  TsumSBBss <- cbind(TsumSBBs[,1],TsumSBBs[,1],TsumSBBs)

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
  nSBB <- which(dataS$forestdamagequalifier=="1602")
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
    
    nSBB <- which(dataS$forestdamagequalifier=="1602")
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
  
  # Add dataS-data and declaration data to the outputs
  #if(ttAll){
  output <- cbind(dataS[,c("segID","N","area","fert","minpeat","landclass",
                              "cons","forestdamagequalifier","dam_id","dam_year",
                              "damSegID","cuttingpurpose",
                              "cuttingrealizationpractice")],
                  output)
  
  # neighboring data
  print(paste("NeighborIDs",neighborIDs))
  if(neighborIDs){
    # add info about neighboring SBB damages  
    colsNeigh <- c((which(names(dataS)=="y")+1):ncol(dataS))
    output <- cbind(output,dataS[,..colsNeigh])
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
  outputs <- cbind(data.table(reg = r_noi),output_mem)
  #return(output_mem)
  #outputs <- rbind(outputs,output_mem)
  if(toFile){ 
    save(outputs,file=paste0(savepath,"SBB_sample_",setname,"_rno",r_noi,"_",regnames[r_noi],".rdata"))
    print(paste0("outputs saved for set",setname," region ",r_noi,"/",regnames[r_noi]))
  }  
  #print(Sys.time()-time0)
  rm(list=setdiff(ls(),toMem))
  gc()
  #return(outputs)
  
} 
    
#############################################################################
ij <- 1
#fmi_from_allas <- F
#if(!exists("fmi_from_allas")) fmi_from_allas <- T
if(!exists("weighted")) weighted <- F

calculateStatistics <- function(ij, fmi_from_allas=F, weighted = F, neighborIDs=F,
                                outputs = outputs, climScen=0, disturbanceON=NA,
                                newSamples=T){
  print(paste("Run climScen",climScen))
  climScen0 <- climScen
  if(weighted) set.seed(10) # for training and testing, keep the samples same
  toMem <- ls()
  r_noi <- rids[ij]
  r_no <- rnos[r_noi]
  print(paste("Region", r_no,"sample"))
  if(newSamples){
    print("Generate new sample")
    sample <- calculateOPSdata(r_noi,nSegs = nSegs, neighborIDs = neighborIDs, weighted = weighted,climScen = climScen)
  } else {
    print("Load sample")
    load(paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/MKIdata/samples_train_valid_",r_noi,".rdata"))
    assign("sample",samples)
  }
  
  if(!newSamples | climScen>0){
    gc()
    
    setid0 <- 1  
    if(onlyValidationset) setid0 <- 2 
    setid <- setid0
    setid1 <- 2
    if(climScen>0) setid1 <- 1 
    
    for(setid in setid0:setid1){ # Go through training and validation sets
      toMem3 <- ls()
      dataS <- sample[[setid]] # setid=1 for training, setid=2 for validation  
      cols_in_prebas <- which(colnames(dataS)%in%c(sample$vars_to_prebas,"x","y"))
      area_tot <- totArea <- sample$areatot
      # PREBAS runs
      source("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Settings_IBCCarbon.R", local=T)  
      setwd("/scratch/project_2000994/PREBASruns/PREBAStesting/")
      rm("data.all")
      gc()
      
      print(paste("Region",r_no,"/",names(sample)[setid]))
      print(paste("SBB area",sum(dataS$area[which(dataS$forestdamagequalifier=="1602")])))
      
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
          id = 1:nrow(dataS),
          E = dataS$x,
          N = dataS$y
        )
        req_coords <- as.matrix(req_coords_dt[, c("E", "N")]) # The coords are passed as a matrix
        
        # Set parameters
        params <- list(req_coords = req_coords, resolution = resolution, years = years)
        
        # Combine arguments
        setup_and_run_args <- c(params, list(save_path = save_path, repo_url = repo_url, format_to_prebas = format_to_prebas))
        
        # RUN
        result <- do.call(setup_and_run, setup_and_run_args)
        
        rm(list = setdiff(ls(), toMemFmi))
        gc()
        
        # Change file name
        setXi <- paste0(setX,"_",r_noi)
        fmi_vars_PREBAS_file <<- paste0("fmi_varsPREBAS",setXi,".rdata")
        climID_lookup_file <<- paste0("climID_lookupPREBAS",setXi,".rdata")
        file.rename(list.files(path=workdir, pattern="fmi_vars_", all.files=FALSE,full.names=FALSE)[1],
                    fmi_vars_PREBAS_file)
        file.rename(list.files(path=workdir, pattern="climID_lookup_", all.files=FALSE,full.names=FALSE)[1],
                    paste0(workdir,climID_lookup_file))
      }
      setwd(workdir)
      
      climatepath_orig = "/scratch/project_2000994/RCP/"
      station_id <- "tmp"
      if(setX==2) station_id <- paste0("tmp2","_",r_noi)
      #climScen <- 0
      nSegs <- nrow(dataS)#dataS)
      print(paste("Spruce = 0 segments", length(which(dataS$spruce==0 & dataS$ba==0))))
      print(paste("SBB area",sum(dataS$area[which(dataS$forestdamagequalifier=="1602")])))
      
      if(!exists("rcps")) rcps <- "CurrClim"
      if(fmi_from_allas) rcps <- "CurrClim_fmi"
      source("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/functions_IBSCarbon.R", local=T)
      source("~/finruns_to_update/functions.R", local=T)
      mortMod <<- 1
      nYears <<- 2050-2015
      endingYear <<- nYears + startingYear
      byManual <- T
      mortMod <<- 1
      
      deltaID=1; sampleID=1; sampleX <-dataS;RCP=0; easyInit=FALSE; 
      initSoilCreStart=NULL; CO2fixed=0;forceSaveInitSoil=F; cons10run = F; procDrPeat=F;coeffPeat1=-240;coeffPeat2=70;coefCH4 = 0.34; coefN20_1 = 0.23;coefN20_2 = 0.077; landClassUnman=NULL;compHarvX = 0;P0currclim=NA;fT0=NA;TminTmax = NA;toRaster=F; ingrowth = F; clcut = 1
      
      #if(!fmi_from_allas){
      rcps0 = "CurrClim"; harvScen="Base"; harvInten="Base"; forceSaveInitSoil=T
      toMem2 <- ls()
      out <-   runModel(1,sampleID=1, outType = outType, 
                        rcps = rcps0,  
                        harvScen="Base",sampleX = dataS, 
                        harvInten="Base", forceSaveInitSoil=T)
      print("grossgrowth")
      print(round(apply(out$region$multiOut[1,1:10,"grossGrowth",,1],1,sum),1))
      print("V")
      print(round(apply(out$region$multiOut[1,1:10,"V",,1],1,sum),1))
      rm(list=setdiff(ls(), toMem2))
      gc()
      
      climScen <- climScen0
      print(paste("SBB area",sum(dataS$area[which(dataS$forestdamagequalifier=="1602")])))
      if(fmi_from_allas) rcps <- "CurrClim_fmi"
      
      if(climScen==0){
        print("Run validation period.")
        nYears <<- 2024-2015
        endingYear <<- nYears + startingYear
        clcuts <<- 1
        disturbanceON <- NA
        if(setid==2) disturbanceON <- disturbanceON0 # "bb" # c("fire","wind","bb")
        source("~/finruns_to_update/functions.R", local=T)
        toMem2 <- ls()
        sampleXs <-   runModel(1,sampleID=1, outType = outType, 
                               rcps = rcps, climScen=climScen,
                               harvScen="NoHarv",sampleX = dataS, 
                               harvInten="NoHarv", 
                               disturbanceON = disturbanceON)
        rm(list=setdiff(ls(),c(toMem2,"sampleXs")))
        gc()
        print("grossgrowth")
        print(round(apply(sampleXs$region$multiOut[1,1:6,"grossGrowth",,1],1,sum),1))
        print("V")
        print(round(apply(sampleXs$region$multiOut[1,1:6,"V",,1],1,sum),1))
        print("bb prob")
        print(sampleXs$region$multiOut[1,1:6,"Rh/SBBpob[layer_1]",1,2])
        print("fire prob")
        print(sampleXs$region$multiOut[1,1:6,"W_wsap/fireRisk[layer_1]",1,2])
        print("wind prob")
        print(sampleXs$region$outDist[1,1:6,"wrisk"])
        
        print(paste("SBB area",sum(dataS$area[which(dataS$forestdamagequalifier=="1602")])))
        #which(sampleXs$region$outDist[,,8]==1,arr.ind=T)
        #print("SMIs:")
        #print(sampleXs$region$multiOut[1,,"NEP/SMI[layer_1]",1,2])
        #print("BBprob:")
        #print(sampleXs$region$multiOut[1,,"Rh/SBBpob[layer_1]",1,2])
        #print("fireprob:")
        #print(sampleXs$region$multiOut[1,,"W_wsap/fireRisk[layer_1]",1,2])
        #print("windprob:")
        #print(sampleXs$region$outDist[1,,"wrisk"])
        
        # SBB simulated damage segments
        #SBBReactionBA <-  apply(sampleXs$region$multiOut[,,"grossGrowth/bb BA disturbed",,2],1:2,sum)
        SBBReactionBA <-  apply(sampleXs$region$multiOut[,,43,,2],1:2,sum)
        any(SBBReactionBA>0)
        Intensitybb <- sampleXs$region$multiOut[,,48,1,2]
        BA <- apply(sampleXs$region$multiOut[,,"BA",,1],1:2,sum)
        Vspruce <- array(unlist(vSpFun(sampleXs$region,2)[,-1]),dim(BA))
        BAspruce <- array(unlist(BASpFun(sampleXs$region,2)[,-1]),dim(BA))
        V <- apply(sampleXs$region$multiOut[,,"V",,1],1:2,sum)
        Vrw <- apply(sampleXs$region$multiOut[,,"VroundWood",,1],1:2,sum)[,-1]
        Vrw <- cbind(Vrw,Vrw[,ncol(Vrw)]) # harvests done next year
        Vrw <- pmax(Vrw,apply(sampleXs$region$multiOut[,,"VroundWood",,1],1:2,sum))
        Ven <- apply(sampleXs$region$multiEnergyWood[,,,1],1:2,sum)[,-1]
        Ven <- cbind(Ven,Ven[,ncol(Ven)])
        Ven <- pmax(Ven,apply(sampleXs$region$multiEnergyWood[,,,1],1:2,sum))
        Vrw <- Vrw+Ven
        
        ## wind
        # all segment areas as initial values for the damages
        areaSamplew <- array(dataS$area,c(dim(SBBReactionBA))) # Segment areas where damage happened
        areaSamplew[sampleXs$region$outDist[,,"damvol"]==0] <- 0 # if no damage happened, damaged area = 0
        #areaSamplebb[SBBReactionBA>0 & Vrw==0] <- areaSamplebb[SBBReactionBA>0 & Vrw==0]*
        #  SBBReactionBA[SBBReactionBA>0 & Vrw==0]/BA[SBBReactionBA>0 & Vrw==0] # if no clearcut, only part of area damaged
        areaSamplew[BA==0] <- 0 # if no basal area, no damaged area
        areaSamplewHarv <- areaSamplew
        areaSamplewHarv[sampleXs$region$outDist[,,"salvlog"]==0] <- 0
        
        ## bb  
        # all segment areas as initial values for the damages
        areaSamplebb <- array(dataS$area,c(dim(SBBReactionBA))) # Segment areas where damage happened
        areaSamplebb[SBBReactionBA==0] <- 0 # if no damage happened, damaged area = 0
        #areaSamplebb[SBBReactionBA>0 & Vrw==0] <- areaSamplebb[SBBReactionBA>0 & Vrw==0]*
        #  SBBReactionBA[SBBReactionBA>0 & Vrw==0]/BA[SBBReactionBA>0 & Vrw==0] # if no clearcut, only part of area damaged
        areaSamplebb[BA==0] <- 0 # if no basal area, no damaged area
        areaSamplebbHarv <- areaSamplebb
        areaSamplebbHarv[Vrw==0] <- 0
        ##
        
        years <- (startingYear+1):endingYear
        nbb <- which(SBBReactionBA>0)
        allDamagesbb <- NA
        if(FALSE & length(nbb>0)){
          areabb <- array(dataS$area,c(dim(SBBReactionBA)))[SBBReactionBA>0] # Segment areas where damage happened
          areabbHarv <- areaSamplebbHarv[SBBReactionBA>0]
          yearsSamplebb <- t(array(years, c(dim(SBBReactionBA)[c(2,1)])))[SBBReactionBA>0] # years for damage
          VSamplebb <- V[SBBReactionBA>0]
          BASamplebb <- BA[SBBReactionBA>0]
          VspruceSamplebb <- c(Vspruce[SBBReactionBA>0])
          BAspruceSamplebb <- c(BAspruce[SBBReactionBA>0])
          damageIntenSamplebb <- Intensitybb[SBBReactionBA>0]
          damageBASamplebb <- SBBReactionBA[SBBReactionBA>0]
          BASamplebb <- BA[SBBReactionBA>0]
          allDamagesbb <- data.table(BAspruceSamplebb,VspruceSamplebb,BASamplebb,VSamplebb,
                                     damageIntenSamplebb,damageBASamplebb,BASamplebb,
                                     areabb,areabbHarv,yearsSamplebb)
        }
        years <- years[years>2018 & years<2024]
        bb_dam_area <- array(0,c(length(years),3),dimnames = list(years,c("sampledata","sim_harvested","sim_all")))
        w_dam_area <- array(0,c(length(years),3),dimnames = list(years,c("sampledata","sim_harvested","sim_all")))
        probs_segm <- array(0,c(length(years),6),dimnames = list(years,c("min_pw_decl_segm","pw_median_decl","pw_median_decl/median_all",
                                                                         "min_pbb_decl_segm","pbb_median_decl","pbb_median_decl/median_all")))
        
        #ppprob <- sampleXs$region$multiOut[,years-2015,"Rh/SBBpob[layer_1]",1,2]
        #indx <- max.col(ppprob, ties.method='first')
        #ppprobmax <- ppprob[cbind(1:nrow(ppprob), indx)]
        
        ti <- 1
        for(ti in 1:length(years)){
          yeari <- years[ti]
          nibb <- which(dataS$forestdamagequalifier=="1602" & as.numeric(dataS$dam_year)==yeari)
          niw <- which(dataS$forestdamagequalifier=="1504" & as.numeric(dataS$dam_year)==yeari)
          pw <- sampleXs$region$outDist[,yeari-2015,"wrisk"]
          pbb <- sampleXs$region$multiOut[,yeari-2015,"Rh/SBBpob[layer_1]",1,2]
          #hist(pbb,50)
          if(length(niw)>0){
            probs_segm[ti,1:3] <- c(min(pw[niw]),median(pw[niw]),
                                    median(pw[niw])/median(pw[setdiff(1:nSegs,niw)]))}
          if(length(nibb)>0){
            probs_segm[ti,4:6] <- c(min(pbb[nibb]),median(pbb[nibb]),
                                    median(pbb[nibb])/median(pbb[setdiff(1:nSegs,nibb)]))
            }
          bb_dam_area[ti,] <- c(sum(dataS$area[which(dataS$forestdamagequalifier=="1602" &
                                                       as.numeric(dataS$dam_year)==yeari)]),
                                #sum(sampleXs$region$multiOut[,yeari-2015,"Rh/SBBpob[layer_1]",1,2]*dataS$area),
                                sum(areaSamplebbHarv[,yeari-2015]),
                                sum(areaSamplebb[,yeari-2015]))
          w_dam_area[ti,] <- c(sum(dataS$area[which(dataS$forestdamagequalifier=="1504" &
                                                      as.numeric(dataS$dam_year)==yeari)]),
                               #sum(sampleXs$region$outDist[1,yeari-2015,"wrisk"]*dataS$area),
                               sum(areaSamplewHarv[,yeari-2015]),
                               sum(areaSamplew[,yeari-2015]))
        }
        print(paste("Region",r_no,"/",regnames[r_noi],": bb damaged segment area"))
        print(bb_dam_area)
        print(colSums(bb_dam_area))
        print(paste("Region",r_no,"/",regnames[r_noi],": wind damaged segment area"))
        print(w_dam_area)
        print(colSums(w_dam_area))
        print(paste("Region",r_no,"/",regnames[r_noi],": decl segment probabilities"))
        print(probs_segm)
        sampleArea <- sum(dataS$area)
        out <- list(bb_dam_area,w_dam_area,probs_segm,regnames[r_noi],allDamagesbb, sampleArea)
        #save(out, file = paste0("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Results/validation_stats_rno",r_no,"_",names(sample)[setid],".rdata"))
        save(out, file = paste0(savepath,"validation_stats_rno",r_noi,"_",names(sample)[setid],".rdata"))
        print(paste0("saved file validation_stats_rno",r_noi,"_",names(sample)[setid],".rdata"))
        
        
        #if(setid==1){ 
        nams <- c("training","validation")
        toMem2 <- ls()
        outputs <- trainingSetCreation(r_noi, sampleXs, dataS, neighborIDs = neighborIDs,
                                       startingYear = startingYear, endingYear= endingYear,
                                       TestaaSBBkoodi=F,nams[setid])
        rm(list=setdiff(ls(),c(toMem2)))
        #}
        
      } else if(climScen>0){
        print("Run scenarios.")
        nYears <<- 2100-2015
        endingYear <<- nYears + startingYear
        disturbanceON <- disturbanceON0 #c("fire","wind","bb")
        HarvScens <- c("NoHarv","baseTapio","Base")
        hi <- 3; climi <- 1
        outputnames <- paste0(rep(HarvScens,each=3),1:3)
        #out <- array(0,c(1+length(outputnames),nYears),dimnames = list(c("totarea",outputnames),(startingYear+1):endingYear))
        
        out <- data.frame()
        for(hi in 1:3){
          for(climi in 1:3){
            toMemiter <- ls()
            harvScen <- HarvScens[hi]
            harvInten <- "Base"
            climScen <- climi
            clcuts <<- 1
            print(paste("climScen changed to",climScen))
            rcps <<- rcpsFile <-paste0(climMod[ClimModid],rcpx[climi])
            rcpsName <- rcps
            source("~/finruns_to_update/functions.R")#, local=T)
            toMem2 <- ls()
            sampleXs <-   runModel(1,sampleID=1, outType = outType, #RCP=climScen,
                                   rcps = rcps, climScen = climi, 
                                   harvScen = harvScen, harvInten = harvInten, #clcut = clcut,
                                   sampleX = dataS, ingrowth=T,
                                   disturbanceON = disturbanceON)
            rm(list=setdiff(ls(),c(toMem2,"sampleXs")))
            gc()
            print("grossgrowth")
            print(round(apply(sampleXs$region$multiOut[1,1:6,"grossGrowth",,1],1,sum),1))
            print("bb prob")
            print(sampleXs$region$multiOut[1,1:6,"Rh/SBBpob[layer_1]",1,2])
            print("fire prob")
            print(sampleXs$region$multiOut[1,1:6,"W_wsap/fireRisk[layer_1]",1,2])
            print("wind prob")
            print(sampleXs$region$outDist[1,1:6,"wrisk"])
            
            print(paste("SBB area",sum(dataS$area[which(dataS$forestdamagequalifier=="1602")])))
            
            # SBB simulated damage segments
            #SBBReactionBA <-  apply(sampleXs$region$multiOut[,,"grossGrowth/bb BA disturbed",,2],1:2,sum)
            SBBReactionBA <-  apply(sampleXs$region$multiOut[,,43,,2],1:2,sum)
            any(SBBReactionBA>0)
            Intensitybb <- sampleXs$region$multiOut[,,48,1,2]
            BA <- apply(sampleXs$region$multiOut[,,"BA",,1],1:2,sum)
            # BA_1 to calculate how big area actually died by sbb
            BA_1 <- apply(sampleXs$region$multiOut[,,"BA",,1],1:2,sum)[,-ncol(BA)]
            BA_1 <- cbind(BA_1[,ncol(BA_1)],BA_1)
            Vmort <- apply(sampleXs$region$multiOut[,,"Vmort",,1],1:2,sum)
            Vspruce <- array(unlist(vSpFun(sampleXs$region,2)[,-1]),dim(BA))
            BAspruce <- array(unlist(BASpFun(sampleXs$region,2)[,-1]),dim(BA))
            V <- apply(sampleXs$region$multiOut[,,"V",,1],1:2,sum)
            deadwoodVolume <- apply(sampleXs$region$multiOut[,,"DeadWoodVolume",,1],1:2,sum)
            grossgrowth <- apply(sampleXs$region$multiOut[,,"grossGrowth",,1],1:2,sum)
            harvNextYear <- T
            if(harvNextYear){
              Vrw <- apply(sampleXs$region$multiOut[,,"VroundWood",,1],1:2,sum)[,-1]
              Vrw <- cbind(Vrw,Vrw[,ncol(Vrw)]) # harvests done next year
              Vrw <- pmax(Vrw,apply(sampleXs$region$multiOut[,,"VroundWood",,1],1:2,sum))
              Ven <- apply(sampleXs$region$multiEnergyWood[,,,1],1:2,sum)[,-1]
              Ven <- cbind(Ven,Ven[,ncol(Ven)])
              Ven <- pmax(Ven,apply(sampleXs$region$multiEnergyWood[,,,1],1:2,sum))
            } else {
              Vrw <- apply(sampleXs$region$multiOut[,,"VroundWood",,1],1:2,sum)
              Ven <- apply(sampleXs$region$multiEnergyWood[,,,1],1:2,sum)
            }
            Vrw <- Vrw+Ven
            
            print("Any nonzero harvests?")
            print(any(Vrw>0))
            #id <- 11
            #V[which(SBBReactionBA[,id]>0),(id-2):(id+3)]
            #Vrw[which(SBBReactionBA[,id]>0),(id-2):(id+3)]
            #SBBReactionBA[which(SBBReactionBA[,id]>0),(id-2):(id+3)]
            
            ## bb  
            # all segment areas as initial values for the damages
            areaSamplebb <- array(dataS$area,c(dim(SBBReactionBA))) # Segment areas where damage happened
            areaSamplebb[SBBReactionBA==0] <- 0 # if no damage happened, damaged area = 0
            #areaSamplebb[SBBReactionBA>0 & Vrw==0] <- areaSamplebb[SBBReactionBA>0 & Vrw==0]*
            #  SBBReactionBA[SBBReactionBA>0 & Vrw==0]/BA[SBBReactionBA>0 & Vrw==0] # if no clearcut, only part of area damaged
            areaSamplebb[BA==0] <- 0 # if no basal area, no damaged area
            
            areaSamplebbHarv <- areaSamplebb
            #areaSamplebbHarv[sampleXs$region$outDist[,,"mgmtreact"]==0]<-0
            areaSamplebbHarv[Vrw==0] <- 0
            
            areaSampleActualbbdamage <- areaSamplebb*SBBReactionBA/BA_1
            areaSampleActualbbdamage[BA_1==0] <- 0
            ##
            
            nbb <- which(SBBReactionBA>0)
            if(FALSE){# length(nbb)>0){
              years <- (startingYear+1):endingYear
              allDamagesbb <- NA
              areabb <- array(dataS$area,c(dim(SBBReactionBA)))[SBBReactionBA>0] # Segment areas where damage happened
              areabbHarv <- areaSamplebbHarv[SBBReactionBA>0]
              yearsSamplebb <- t(array(years, c(dim(SBBReactionBA)[c(2,1)])))[SBBReactionBA>0] # years for damage
              VSamplebb <- V[SBBReactionBA>0]
              BASamplebb <- BA[SBBReactionBA>0]
              VspruceSamplebb <- c(Vspruce[SBBReactionBA>0])
              BAspruceSamplebb <- c(BAspruce[SBBReactionBA>0])
              damageIntenSamplebb <- Intensitybb[SBBReactionBA>0]
              damageBASamplebb <- SBBReactionBA[SBBReactionBA>0]
              BASamplebb <- BA[SBBReactionBA>0]
              allDamagesbb <- data.table(BAspruceSamplebb,VspruceSamplebb,BASamplebb,VSamplebb,
                                         damageIntenSamplebb,damageBASamplebb,BASamplebb,
                                         areabb,areabbHarv,yearsSamplebb)
            } 
            if(hi==1 & climi==1) out <- data.table(c(scen = "all",var="area", totArea = rep(sum(dataS$area),nYears)))
            out <- cbind(out, data.table(c(scen = paste0(harvScen,"_clim",climScen),
                                           var ="bbdamage", bbdamage = colSums(areaSampleActualbbdamage))))
            out <- cbind(out, data.table(c(scen = paste0(harvScen,"_clim",climScen),
                                           var ="bbHarvsegm", bbHarvsegm = colSums(areaSamplebbHarv))))
            out <- cbind(out, data.table(c(scen = paste0(harvScen,"_clim",climScen),
                                           var ="bbsegm", bbsegm = colSums(areaSamplebb))))
            out <- cbind(out, data.table(c(scen = paste0(harvScen,"_clim",climScen),
                                           var ="V", V = colSums(V*dataS$area)/sum(dataS$area))))
            out <- cbind(out, data.table(c(scen = paste0(harvScen,"_clim",climScen),
                                           var ="Vspruce", Vspruce = colSums(Vspruce*dataS$area)/sum(dataS$area))))
            out <- cbind(out, data.table(c(scen = paste0(harvScen,"_clim",climScen),
                                           var ="grossgrowth", grossgrowth = colSums(grossgrowth*dataS$area)/sum(dataS$area))))
            Vrw <- apply(sampleXs$region$multiOut[,,"VroundWood",,1],1:2,sum)
            Ven <- apply(sampleXs$region$multiEnergyWood[,,,1],1:2,sum)
            Vrw <- Vrw+Ven
            out <- cbind(out, data.table(c(scen = paste0(harvScen,"_clim",climScen),
                                           var ="Vharv", Vharv = colSums(Vrw*dataS$area)/sum(dataS$area))))
            out <- cbind(out, data.table(c(scen = paste0(harvScen,"_clim",climScen),
                                           var ="Vmort", Vmort = colSums(Vmort*dataS$area)/sum(dataS$area))))
            out <- cbind(out, data.table(c(scen = paste0(harvScen,"_clim",climScen),
                                           var ="deadwoodVolume", deadwoodVolume = colSums(deadwoodVolume*dataS$area)/sum(dataS$area))))
            print(head(out))
            print(tail(out))
            # print(out[,7])
            rm(list=setdiff(ls(),toMemiter))
            gc()
          }
        }
        filee <- 
          if(is.na(disturbanceON[1])){
            save(out,file = paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/bbScenarios/output_",r_noi,"_",regnames[r_noi],"_nodist.rdata"))
          } else {
            save(out,file = paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/bbScenarios/output_",r_noi,"_",regnames[r_noi],".rdata"))
          }
      }
      
      rm(list=setdiff(ls(),c(toMem3)))
      gc()
    }
    #return(outputs)
    rm(list=setdiff(ls(),c(toMem)))#,"out","outputs")))
    gc()
    if(fmi_from_allas){
      file.remove(paste0(workdir,fmi_vars_PREBAS_file))
      file.remove(paste0(workdir,climID_lookup_file))
    }
  }
}
#calculateStatistics(1, fmi_from_allas = fmi_from_allas, weighted = F)
if(!sbatches){
  output_stats <- lapply(1:length(rids), function(jx) {
    print(paste("ClimScen value", climScen))
    #     print(paste0("region list: ",which(rids==20),"/",length(rids)))
    calculateStatistics(jx, fmi_from_allas = fmi_from_allas, neighborIDs=neighborIDs,
                        weighted = weighted, climScen = climScen, 
                        disturbanceON = disturbanceON0, newSamples = newSamples)
  })      
} else {
  library(parallel)
  output_stats <- mclapply(1:length(rids), function(jx) {
    print(paste("ClimScen value", climScen))
    #     print(paste0("region list: ",which(rids==20),"/",length(rids)))
    calculateStatistics(jx, fmi_from_allas = fmi_from_allas, neighborIDs=neighborIDs,
                        weighted = weighted, climScen = climScen, 
                        disturbanceON = disturbanceON0, newSamples = newSamples)
  }, mc.cores = 5,mc.silent=FALSE)      
}

save(output_stats, file="testitiedosto.rdata")
#if(fmi_from_allas){
#  file.remove(paste0(workdir,"fmi_vars_PREBAS.rdata"))
#  file.remove(paste0(workdir,"climID_lookup.rdata"))
#}
#save(output_stats, file = "/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Results/validation_stats.rdata")
#any(sampleXs$regio$multiOut[,,"grossGrowth/bb BA disturbed",,2]>0)

