#rm(list=ls())
#gc()
if(dev.interactive()) dev.off()
outType <- "testRun"
neighborIDs <- T
toFile <- F
asParallel <- F
if(neighborIDs) asParallel <- T
if(!exists("toFile")) toFile <- T
nSegs <- 10000
if(toFile) nSegs <- 20000
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

vPREBAS <- "newVersion" #
#vPREBAS <- "master"   

savepath = "/scratch/project_2000994/PREBASruns/PREBAStesting/MKIdata/"

ContinueIterations <- F

#r_noi <- r_no
rnos <- c(1:8,8:19)
rids <- rids0 <- c(1,3:length(rnos))
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
calculateOPSdata  <-  function(r_noi, neighborIDs=T, weighted = T){
  toMem <- ls()
  r_no <- rnos[r_noi]
  
  print(paste("region",r_no))
  print(paste("Neighbor information =",neighborIDs))
  
  fname <- paste0("DeclaredDamages_",dam_names[inds],"_rno",r_no,"_",regnames[r_noi],".rdata")
  load(file=paste0(savepath,"/",fname))
  print(paste("File",fname,"opened"))

  tt <- XYdamages[dam_year<2024,]   
  if(!weighted) tt <- tt[tt$dam_year>2014,]
  tt$pine[is.na(tt$pine)] <- 0
  tt$spruce[is.na(tt$spruce)] <- 0
  tt$birch[is.na(tt$birch)] <- 0
  ttsum = tt$pine+tt$spruce+tt$birch
  tt[,pine := pine/ttsum]
  tt[,spruce := spruce/ttsum]
  tt[,birch:= birch/ttsum]
  tt[ttsum==0,c("pine","spruce","birch")]<-0
  XYdamages <- tt
  
  # PREBAS run for a sample
  landClassX <- 1:2
  #print(r_no)
  
  #    devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/settings.r")
  source("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Settings_IBCCarbon.R", local=TRUE)  
  vars_to_prebas <- colnames(data.all)
  
  IDsUniq <- unique(XYdamages[,c("dam_id","segID")])
  IDsUniq[,damSegID:=1:nrow(IDsUniq)]
  
  # Give segments that are split with declaration polygons, a unique name
  XYdamages[,damSegID := match(paste(XYdamages$dam_id,XYdamages$segID),
                               paste(IDsUniq$dam_id,IDsUniq$segID))]
  areas <- XYdamages[, .N, by = list(damSegID)]
  
  tt <-XYdamages[XYdamages[,.I[which.max(spruce)], by=damSegID]$V1] # one row for each segment
  tt[,area := areas[match(areas$damSegID, tt$damSegID),N]*16^2/100^2]
  
  if(weighted) tt <- tt[tt$dam_year>=2019,]
  tt <- tt[segID%in%data.all$segID,] # some damage polygons are outside the landclasses 1:2
  #tt[,climID:=data.all$climID[match(tt$segID,data.all$segID)]] # cliID from data.all
  #tt[,cons:=data.all$cons[match(tt$segID,data.all$segID)]] # cons from data.all
  
  not_in_tt <- colnames(data.all)[which(!(colnames(data.all)%in%colnames(tt)))]
  tt <- cbind(tt,data.all[match(tt$segID,data.all$segID),..not_in_tt])
  #tt[,cons:=data.all$cons[match(tt$segID,data.all$segID)]] # cliID from data.all
  
  #########################
  cuttinginpractise <- c(5, 8, 16, 17, 19, 21, 22, 24)
  data.allA <- cbind(data.all[c(match(tt$segID,data.all$segID)),],
                     forestdamagequalifier = tt$forestdamagequalifier,              
                     cuttingpurpose = tt$cuttingpurpose,              
                     cuttingrealizationpractice = tt$cuttingrealizationpractice,
                     dam_year = tt$dam_year,
                     damSegId = tt$damSegID,
                     dam_id = tt$dam_id)
  data.allA <- tt
  yrs <- 2019:2023
  #c(5,#avohakkuu
  #                       22) #hyonteistuhoalue, uudistamishakkuu
  
  areasDamClct <- function(x,id) sum(data.allA$area[which(data.allA$forestdamagequalifier==id & data.allA$dam_year==x  & data.allA$cuttingrealizationpractice%in%cuttinginpractise)])
  areasDam <- function(x,id) sum(data.allA$area[which(data.allA$forestdamagequalifier==id & data.allA$dam_year==x)])
  
  damInfo <- array(0,c(8,length(yrs)),
                   dimnames = list(c("alldeclarea","alldeclclctarea","SBBarea","SBBclctarea","windarea",
                                     "windclctarea","fireclctarea","firearea"),paste0("Year",yrs)))
  for(y in 1:length(yrs)) damInfo[1,y] <- sum(data.allA$area[which(data.allA$dam_year==yrs[y])])
  for(y in 1:length(yrs)) damInfo[2,y] <- sum(data.allA$area[which(data.allA$dam_year==yrs[y] & data.allA$cuttingrealizationpractice%in%cuttinginpractise)])
  for(y in 1:length(yrs)) damInfo[3,y] <- areasDam(yrs[y],dam_indexs[which(dam_names=="SBB")])
  for(y in 1:length(yrs)) damInfo[4,y] <- areasDamClct(yrs[y],dam_indexs[which(dam_names=="SBB")])
  for(y in 1:length(yrs)) damInfo[5,y] <- areasDam(yrs[y],dam_indexs[which(dam_names=="wind")])
  for(y in 1:length(yrs)) damInfo[6,y] <- areasDamClct(yrs[y],dam_indexs[which(dam_names=="wind")])
  for(y in 1:length(yrs)) damInfo[7,y] <- areasDam(yrs[y],dam_indexs[which(dam_names=="fire")])
  for(y in 1:length(yrs)) damInfo[8,y] <- areasDamClct(yrs[y],dam_indexs[which(dam_names=="fire")])
  if(weighted) save(damInfo,file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/damInfo_",r_noi,".rdata"))

  ########################
  nt <- 1:nrow(tt)
  # Pick samples of disturbances in the sample
  ni <- which(tt$forestdamagequalifier[nt]%in%c("1602","1504","1503"))
  print(paste("all disturbed (wind, fire, bb) segments n =",length(ni)))
  nSegs2 <- nSegsd <- nSegs/2
  if(ttAll)  nSegsd <- 3*nSegs/8
  if(length(ni)>nSegsd) ni <- sample(ni,nSegsd,replace = F)
  
  ## bb samples
  #nitmp <- which(tt$forestdamagequalifier[nt]%in%c("1602"))
  #print(paste("dist. bb n =",length(nitmp)))
  # # if(length(nitmp)>nSegs2) nitmp <- sample(nitmp,nSegs2,replace = F)#[1:nSegs2]
  #ni <- nitmp
  ## wind  
  #nitmp <-which(tt$forestdamagequalifier[nt]%in%c("1504"))
  #print(paste("dist. wind n =",length(nitmp)))
  ##if(length(nitmp)>nSegs2) nitmp <- sample(nitmp,nSegs2,replace = F)#[1:nSegs2]
  #ni <- c(ni,nitmp)
  ## fire
  #nitmp <-which(tt$forestdamagequalifier[nt]%in%c("1503"))
  #print(paste("dist. fire n =",length(nitmp)))
  ##if(length(nitmp)>nSegs2) nitmp <- sample(nitmp,nSegs2,replace = F)#[1:nSegs2]
  #ni <- c(ni,nitmp)
  
  if(weighted){
    nt <- nt[ni]
    ni <- c(nt,sample(setdiff(1:nrow(tt),ni),nSegs2-min(nSegs2,length(ni))))
  }
  # Add data outside the declarations as sample set
  load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
  data.IDs$segID <- data.IDs$maakuntaID
  data.IDs <- data.IDs[segID!=0]
  setkey(data.IDs,segID)
  gc()
  setkey(data.all,segID)
  
  if(r_no %in% c(8,9)){ # Lappi E and P: divide data.all according to y-coordinate
    tabX <- merge(data.IDs,data.all) # coords of the segments in sample outside declarations
    x <- tabX$x[match(data.all$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"x"]
    y <- tabX$y[match(data.all$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"y"]
    
    xybb <- array(0,c(nrow(XYdamages),2))  
    xybb[,1] <- XYdamages$x
    xybb[,2] <- XYdamages$y
    xx <- bbox(xybb)
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
  
  if(weighted){
    rrows <- which(!(data.all$segID %in% XYdamages$segID))# niAll
    niAll <- sample(rrows,nSegs2)
    tt2 <- data.table(matrix(NA,nrow = nSegs2,ncol = ncol(tt)))
    colnames(tt2) <- colnames(tt)
    tt2$forestdamagequalifier = 0
    tt2$cuttingrealizationpractice <- 0
    tt2$cuttingpurpose = 0
    tt2$dam_year <- sample(unique(tt$dam_year),nSegs2,replace = T)
    tt2$dam_id <- tt2$damSegID <- paste0("0",1:nSegs2)
    
    dataS <- cbind(data.all[c(match(tt$segID[ni],data.all$segID),niAll),],
                   forestdamagequalifier = c(tt$forestdamagequalifier[ni],tt2$forestdamagequalifier),              
                   cuttingpurpose = c(tt$cuttingpurpose[ni],tt2$cuttingpurpose),              
                   cuttingrealizationpractice = c(tt$cuttingrealizationpractice[ni],tt2$cuttingrealizationpractice),
                   dam_year = c(tt$dam_year[ni],tt2$dam_year),
                   damSegId = c(tt$damSegID[ni],tt2$damSegID),
                   dam_id = c(tt$dam_id[ni],tt2$dam_id))
    dataS[1:nSegs2,"area"] <- tt$area[ni] # areas of the MKI based segments in the declaration data
    
    setkey(dataS,segID)
    
    tabX <- merge(data.IDs,dataS) # coords of the segments in sample outside declarations
    rm(list = "data.IDs")
    gc()
    x <- tabX$x[match(dataS$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"x"]
    y <- tabX$y[match(dataS$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"y"]
    ni_sample <- which(dataS$forestdamagequalifier==0)
    ni_decl <- setdiff(1:nrow(dataS),ni_sample)
    dataS <- cbind(dataS,data.table(x,y))
    dataS[ni_decl,c("x","y")] <- XYdamages[match(dataS$dam_id[ni_decl],XYdamages$dam_id),c("x","y")]
    
    set_thin_PROJ6_warnings(TRUE)
    xy <- dataS[,c("segID","x","y")]
    coordinates(xy) <- c("x","y")
    proj4string(xy) <- crsX
    #cord = SpatialPoints(xy, proj4string=CRS("+init=EPSG:3067"))
    location<-as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
    dataS$lat <- location$y
    
    Ntot <- sum(data.all$area)/(0.16^2)
    gc()
    
  } else {
    tt_names_in <- which(colnames(tt)%in%colnames(data.all))
    dataAll_names_in <- which(colnames(data.all)%in%colnames(tt))
    
    tmp_inds <- 1:(nrow(tt)+nrow(data.all))
    ni <- sample(1:length(tmp_inds), nSegs, replace=F)
    tmp <- rbind(tt[,..tt_names_in], data.all[,..dataAll_names_in])[ni,]
    
    rrows_decl <- which(ni<=nrow(tt))
    #rrows_decl <- which((tmp$segID %in% tt$segID))# niAll
    rrows_sample <- setdiff(1:length(ni), rrows_decl)
    #rrows_sample <- which(!(tmp$segID %in% tt$segID))# niAll

    setkey(tmp,segID)
    tabX <- merge(data.IDs,tmp) # coords of the segments in sample outside declarations
    x <- tabX$x[match(tmp$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"x"]
    y <- tabX$y[match(tmp$segID,tabX$segID)]#tabX[tabX[,I(which.max(y)),by=damSegId]$V1,"y"]
    tmp[,x:=x]
    tmp[,y:=y]
    
    #tt2 <- data.table(matrix(NA,nrow = length(rrows_sample),ncol = ncol(tt)))
    #colnames(tt2) <- colnames(tt)
    #tt2$forestdamagequalifier = 0
    #tt2$cuttingrealizationpractice <- 0
    #tt2$cuttingpurpose = 0
    #tt2$dam_year <- sample(unique(tt$dam_year),nrow(tt2),replace = T)
    #tt2$dam_id <- tt2$damSegID <- paste0("0",1:nrow(tt2))
    
    
    tmp <- cbind(tmp, forestdamagequalifier = 0, cuttingrealizationpractice = 0,
                 cuttingpurpose = 0, 
                 dam_year = sample(unique(tt$dam_year),nrow(tmp),replace = T),
                 dam_id = paste0("0",1:nrow(tmp)),
                 damSegID = paste0("0",1:nrow(tmp)))
    #tt_names_in <- which(colnames(tt)%in%colnames(tmp))
    #tmp_names_in <- which(colnames(tmp)%in%colnames(tt))
#tt[ni[rrows_decl],]
  #  rrows <- which(tt$segID %in% tmp$segID[rrows_decl])
    #ttmp <- tt[match(tt$segID[rrows],tmp$segID[rrows_decl]),]
    tmp[rrows_decl, "forestdamagequalifier"] <- tt$forestdamagequalifier[ni[rrows_decl]]
    tmp[rrows_decl, "cuttingrealizationpractice"] <- tt$cuttingrealizationpractice[ni[rrows_decl]]
    tmp[rrows_decl, "cuttingpurpose"] <- tt$cuttingpurpose[ni[rrows_decl]]
    tmp[rrows_decl, "dam_year"] <- tt$dam_year[ni[rrows_decl]]
    tmp[rrows_decl, "damSegID"] <- tt$damSegID[ni[rrows_decl]]
    tmp[rrows_decl, "dam_id"] <- tt$dam_id[ni[rrows_decl]]
    tmp[rrows_decl, "climID"] <- tt$climID[ni[rrows_decl]]
    
    #colnames(tt)[match(colnames(tt)[tt_names_in], colnames(tmp))]
    #colnames(tmp)[tmp_names_in]
    
    #tmp[rrows_decl,] <- tt[match(tt$segID[rrows],tmp$segID[rrows_decl]),]
    #dataS <- tt[match(tt$segID[rrows],tmp$segID[rrows_decl]),]
    
    dataS <- tmp
  }
  #  load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
  #  data.IDs$segID <- data.IDs$maakuntaID
  #  data.IDs <- data.IDs[segID!=0]
  #  setkey(data.IDs,segID)
  #  gc()
  #print(range(data.all$segID))
  #print(range(dataS$segID))
  #if(!toFile) print(dataS$x)
  if(length(which(is.na(dataS$x)))>0){
    nNas <- which(is.na(dataS$x))
    print(dataS[nNas,])
  }
  #####
  # For the sample, indexes about neighbors
  print(paste("Calculate neighbor information =",neighborIDs))
  if(neighborIDs & weighted){
    # all declarations, all coordinates and years of interest
    ntmp <- which(XYdamages$cuttingrealizationpractice%in%cuttinginpractise | 
                    XYdamages$forestdamagequalifier==dam_indexs[dam_names=="SBB"] |
                    XYdamages$forestdamagequalifier==dam_indexs[dam_names=="wind"])
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
    TEST<-T
    if(TEST){ # use the function formulation of the neighbor validation
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
      source("../PREBAStesting/0.5_functions.R", local = T)
      outputNeighbor <- apply(data.table(c(1:nSegs)),1,neighborsAll,dataS=dataS,
                              damCCutInt=damCCutInt, damBBInt=damBBInt, damWInt=damWInt)
      outputNeighbor <- data.table(t(outputNeighbor))
      colnames(outputNeighbor) <- dimNams
      outputNeighbor[outputNeighbor==1e12] <- NA
      #rm(list=c("xx","yy","dx"))
      gc()
      #print(dataS$x)
      dataS <- cbind(dataS,outputNeighbor)
    } else {
      for(ij in 1:nSegs){
        if((id_ij+1)%%min(1000,round(nSegs/4))==1){
          timeT2 <- Sys.time()
          tperiter<- timeT2-timeT
          print(paste0("neighborinfo rno",r_noi,": ",id_ij,"/",nSegs,
                       ", time ",round(tperiter,2),
                       ", segment sizes: ",round(sampleSets/(id_ij-1),2),
                       ", neightborsets: ",round(neighborSets/(id_ij-1),2),
                       ", neighbors: ",
                       round(sum(as.integer(clearcutNeighbor[1:(id_ij-1),1]<1e10))/(id_ij-1)*100,2),
                       "%"))
          neighborSets <- sampleSets <- 0
          timeT <- timeT2
        }
        id_ij <- id_ij+1
        
        damSegIDij <- dataS$damSegId[ij]
        if(is.na(dataS$forestdamagequalifier[ij]) | dataS$forestdamagequalifier[ij]>0){ # in declarations
          iks <- which(XYdamages$damSegID==damSegIDij)
          dam_yeari <- as.numeric(XYdamages$dam_year[iks[1]]) # damage year of the dam_id
          xi <- XYdamages$x[iks]
          yi <- XYdamages$y[iks]
        } else if(dataS$forestdamagequalifier[ij]==0){          # in sample
          iks <- which(tabX$damSegId==damSegIDij) # all segments in the dam_id
          dam_yeari <- as.numeric(tabX$dam_year[iks[1]]) # damage year of the dam_id
          xi <- tabX$x[iks]
          yi <- tabX$y[iks]
        } else { print("error")}
        sampleSets <- sampleSets+length(iks)
        #if(is.na(dam_yeari)) print(paste(ij,dam_indd[iks],dam_yeari))
        # choose declarations from years dam_yeari-4:dam_yeari-1 and thus also not in declaration ij
        KUVA <- F
        clctDist <- 16*5 # poissible range to be checked
        if(is.na(dam_yeari)){
          ijsd <- NULL  
        } else {
          ijsd <- which(dam_yeard %in% c((dam_yeari-15):(dam_yeari-1)) &
                          xxd<=(max(xi)+clctDist) & xxd>=(min(xi)-clctDist) &
                          yyd<=(max(yi)+clctDist) & yyd>=(min(yi)-clctDist))
        }
        neighborSets <- neighborSets+length(ijsd)
        ijs <- ijsSBB <- ijsWind <- NULL
        if(length(ijsd)>0){
          # clearcuts in neighborhood
          ntmp <- ijsd[which(dam_crpd[ijsd]==5)]
          if(length(ntmp)>0){
            ijs <- ntmp
            if(length(ijs)>0 & KUVA){
              print(paste0("CC neighbors found, ij=",ij))
              par(mfrow=c(1,1))
              plot(c(xxd[ijs],xi)-xi[1],c(yyd[ijs],yi)-yi[1],pch=19,
                   main=paste0("CC, row=",ij,", damSegid=",damSegIDij))
              points(xi-xi[1],yi-yi[1],col="red",pch=19)
              points(xxd[ijsd]-xi[1],yyd[ijsd]-yi[1],col="red",pch=4,cex=1.6)
            }
          }
          # SBB in neighborhood
          ntmp <- ijsd[which(dam_indd[ijsd]==dam_indexs[dam_names=="SBB"])]
          if(length(ntmp)>0){
            ijsSBB <- ntmp
            if(length(ijsSBB)>0 & KUVA){
              print(paste0("SBB neighbors found, ij=",ij))
              par(mfrow=c(1,1))
              plot(c(xxd[ijsSBB],xi)-xi[1],c(yyd[ijsSBB],yi)-yi[1],pch=19,
                   main=paste0("SBB, row=",ij,", damSegid=",damSegIDij))
              points(xi-xi[1],yi-yi[1],col="red",pch=19)
              points(xxd[ijsd]-xi[1],yyd[ijsd]-yi[1],col="red",pch=4,cex=1.6)
            }
          }
          # Wind damage in neighborhood
          ntmp <- ijsd[which(dam_indd[ijsd]==dam_indexs[dam_names=="wind"])]
          if(length(ntmp)>0){
            ijsWind <- ntmp
            if(length(ijsWind)>0 & KUVA){
              print(paste0("Wind neighbors found, ij=",ij))
              par(mfrow=c(1,1))
              plot(c(xxd[ijsWind],xi)-xi[1],c(yyd[ijsWind],yi)-yi[1],pch=19,
                   main=paste0("Wind, row=",ij,", damSegid=",damSegIDij))
              points(xi-xi[1],yi-yi[1],col="red",pch=19)
              points(xxd[ijsd]-xi[1],yyd[ijsd]-yi[1],col="red",pch=4,cex=1.6)
            }
          }
        }
        dclct <- dclct_south <- dclct_SBB <- dclct_Wind <- c(1e12,NA)
        if(length(ijs)>0){
          #if(length(ijs)>10) break()
          #ti <- Sys.time()
          ik<-1
          
          xx <- xxd[ijs]
          yy <- yyd[ijs]
          dam_year <- dam_yeard[ijs]
          
          if(TEST){
            neighbors <-function(ik,clctDist){
              dclct <- dclct_south <- c(1e12,NA)
              dx <- sqrt((xx-xi[ik])^2+(yy-yi[ik])^2)
              if(dx[order(dx)[1]]<dclct[1] & dx[order(dx)[1]]<clctDist & dx[order(dx)[1]]>0){
                dclct <- c(dx[order(dx)[1]],dam_year[order(dx)[1]]-dam_yeari)  
              }
              # close south neighbours with clearcut: distance and time from clearcut
              dy <- yy[which(xx<(xi[ik]+32) & xx>(xi[ik]-32))]
              dy <- yi[ik]-dy[dy<=yi[ik]]
              if(length(dy)>0){
                if(dy[order(dy)[1]]<dclct_south[1] & dy[order(dy)[1]]<clctDist){#dy[order(dy)[1]] 
                  dclct_south <- c(dy[order(dy)[1]],dam_year[order(dy)[1]]-dam_yeari)
                }
              }
              return(c(dclct,dclct_south))
            }
            tmp <- apply(data.table(1:length(iks)),1:2,neighbors,clctDist=clctDist)
            dclct <- tmp[1:2,which.min(tmp[1,,1]),1]
            dclct_south <- tmp[3:4,which.min(tmp[3,,1]),1]
            
          } else {
            for(ik in 1:length(iks)){
              # close neighbours with clearcut: distance and time from clearcut
              dx <- sqrt((xx-xi[ik])^2+(yy-yi[ik])^2)
              if(dx[order(dx)[1]]<dclct[1] & dx[order(dx)[1]]<clctDist & dx[order(dx)[1]]>0){
                dclct <- c(dx[order(dx)[1]],dam_year[order(dx)[1]]-dam_yeari)  
              }
              # close south neighbours with clearcut: distance and time from clearcut
              dy <- yy[which(xx<(xi[ik]+32) & xx>(xi[ik]-32))]
              dy <- yi[ik]-dy[dy<=yi[ik]]
              if(length(dy)>0){
                if(dy[order(dy)[1]]<dclct_south[1] & dy[order(dy)[1]]<clctDist){#dy[order(dy)[1]] 
                  dclct_south <- c(dy[order(dy)[1]],dam_year[order(dy)[1]]-dam_yeari)
                }
              }
            }
          }
          #print(Sys.time()-ti)
          
          #print(paste(ij,dclct))
        }
        if(length(ijsSBB)>0){
          ik<-1
          # close neighbours with SBB: distance and time from clearcut
          xx <- xxd[ijsSBB]
          yy <- yyd[ijsSBB]
          dam_year <- dam_yeard[ijsSBB]
          dx <- sqrt((xx-xi[ik])^2+(yy-yi[ik])^2)
          
          if(TEST){
            neighbors <-function(ik,clctDist){
              dclct_SBB <- c(1e12,NA)
              dx <- sqrt((xx-xi[ik])^2+(yy-yi[ik])^2)
              if(dx[order(dx)[1]]<dclct_SBB[1] & dx[order(dx)[1]]<clctDist & 
                 dx[order(dx)[1]]>0){
                dclct_SBB <- c(dx[order(dx)[1]],dam_year[order(dx)[1]]-dam_yeari)  
              }
              return(dclct_SBB)
            }
            tmp <- apply(data.table(1:length(iks)),1:2,neighbors,clctDist=clctDist)
            dclct_SBB <- tmp[1:2,which.min(tmp[1,,1]),1]
          } else {
            for(ik in 1:length(iks)){
              #SBB <- SBB_dam_ind[ijs]
              if(dx[order(dx)[1]]<dclct_SBB[1] & dx[order(dx)[1]]<clctDist & 
                 dx[order(dx)[1]]>0){
                dclct_SBB <- c(dx[order(dx)[1]],dam_year[order(dx)[1]]-dam_yeari)  
              }
            }
          }
        }
        if(length(ijsWind)>0){
          ik<-1
          xx <- xxd[ijsWind]
          yy <- yyd[ijsWind]
          dam_year <- dam_yeard[ijsWind]
          
          if(TEST){
            neighbors <-function(ik,clctDist){
              dclct_Wind <- c(1e12,NA)
              dx <- sqrt((xx-xi[ik])^2+(yy-yi[ik])^2)
              if(dx[order(dx)[1]]<dclct_Wind[1] & dx[order(dx)[1]]<clctDist & 
                 dx[order(dx)[1]]>0){
                dclct_Wind <- c(dx[order(dx)[1]],dam_year[order(dx)[1]]-dam_yeari)  
              }
              return(dclct_Wind)
            }
            tmp <- apply(data.table(1:length(iks)),1:2,neighbors,clctDist=clctDist)
            dclct_Wind <- tmp[1:2,which.min(tmp[1,,1]),1]
          } else {
            for(ik in 1:length(iks)){
              # close neighbours with clearcut because of SBB: distance and time from clearcut
              dx <- sqrt((xx-xi[ik])^2+(yy-yi[ik])^2)
              if(dx[order(dx)[1]]<dclct_Wind[1] & dx[order(dx)[1]]<clctDist & 
                 dx[order(dx)[1]]>0){
                dclct_Wind <- c(dx[order(dx)[1]],dam_year[order(dx)[1]]-dam_yeari)  
              }
            }
          }
        }    
        # To memory for lines of dam_id == ij
        clearcutNeighbor[ij,1:2]<-dclct#matrix(dclct,length(ijs),2,byrow = T)
        clearcutNeighbor_SBB[ij,1:2]<-dclct_SBB#matrix(dclct_SBB,length(ijs),2,byrow = T)
        clearcutNeighbor_wind[ij,1:2]<-dclct_Wind#matrix(dclct_Wind,length(ijs),2,byrow = T)
        clearcutNeighbor_south[ij,1:2]<-dclct_south#matrix(dclct_south,length(ijs),2,byrow = T)
        #  print(clearcutNeighbor)
      }
      clearcutNeighbor <- data.table(minDist=clearcutNeighbor[,1],clearcutYear = clearcutNeighbor[,2])
      clearcutNeighbor_SBB <- data.table(minDist_SBB=clearcutNeighbor_SBB[,1],clearcutYear_SBB = clearcutNeighbor_SBB[,2])
      clearcutNeighbor_wind <- data.table(minDist_Wind=clearcutNeighbor_wind[,1],clearcutYear_Wind = clearcutNeighbor_wind[,2])
      clearcutNeighbor_south <- data.table(minDist_s=clearcutNeighbor_south[,1],
                                           clearcutYear_s = clearcutNeighbor_south[,2])
      
      clearcutNeighbor[clearcutNeighbor$minDist==1e12,1] <- NA 
      clearcutNeighbor_south[clearcutNeighbor_south$minDist==1e12,1]=NA
      clearcutNeighbor_SBB[clearcutNeighbor_SBB$minDist_SBB==1e12,1]=NA
      clearcutNeighbor_wind[clearcutNeighbor_wind$minDist_Wind==1e12,1]=NA
      rm(list=c("xx","yy","dx"))
      gc()
      #print(dataS$x)
      dataS <- cbind(dataS,data.table(clearcutNeighbor,clearcutNeighbor_SBB,
                                      clearcutNeighbor_south,
                                      clearcutNeighbor_wind))
    }
    
    #print(dataS$x)
    #hist(na.omit(clearcutNeighbor[,1]))
  }
  if(!toFile) head(dataS)
  if(weighted) save(dataS,file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/MKIdata/samples_",r_noi,".rdata"))
  areatot <- sum(data.all$area)
  xx <- list(dataS, damInfo,areatot, vars_to_prebas)
  names(xx) <- c("dataS","damInfo","areatot","vars_to_prebas")
  return(xx)
  rm(list=setdiff(ls(),toMem))
  gc()
  
}

ij <- 1
#fmi_from_allas <- F
if(!exists("fmi_from_allas")) fmi_from_allas <- T
if(!exists("weighted")) weighted <- F
calculateStatistics <- function(ij, fmi_from_allas=F, weighted = F){
  set.seed(1)
  toMem <- ls()
  r_noi <- rids[ij]
  r_no <- rnos[r_noi]
  
  sample <- calculateOPSdata(r_noi,neighborIDs = F, weighted = weighted)
  dataS <- sample$dataS
  cols_in_prebas <- which(colnames(dataS)%in%c(sample$vars_to_prebas,"x","y"))
  ops <<- list(dataS)#list(dataS[,..cols_in_prebas])
  
  setwd("/scratch/project_2000994/PREBASruns/PREBAStesting/")
  
  print(paste("Region",r_no))
  print(paste("SBB area",sum(ops[[1]]$area[which(dataS$forestdamagequalifier=="1602")])))

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
  setwd(workdir)
  
  # PREBAS runs
  source("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Settings_IBCCarbon.R", local=T)  
  area_tot <<- sum(data.all$area)
  #rm(list="data.all")
  #gc()
  climatepath_orig = "/scratch/project_2000994/RCP/"
  station_id <- "tmp"
  climScen <- 0
  nSegs <- nrow(ops[[1]])
  print(paste("Spruce = 0 segments", length(which(dataS$spruce==0 & dataS$ba==0))))
  print(paste("SBB area",sum(dataS$area[which(dataS$forestdamagequalifier=="1602")])))
  
  rcps <- "CurrClim"
  if(fmi_from_allas) rcps <- "CurrClim_fmi"
  #    source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
  #source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
  source("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/functions_IBSCarbon.R", local=T)
  source("~/adaptFirst_runs/functions.R", local=T)
  #source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
  mortMod <<- 1
  nYears <<- 2050-2015
  endingYear <<- nYears + startingYear
  byManual <- T
  mortMod <<- 1
  
  deltaID=1; sampleID=1; climScen=0; easyInit=FALSE; CO2fixed=0;forceSaveInitSoil=F; cons10run = F; procDrPeat=F;coeffPeat1=-240;coeffPeat2=70;coefCH4 = 0.34; coefN20_1 = 0.23;coefN20_2 = 0.077; landClassUnman=NULL;compHarvX = 0;P0currclim=NA;fT0=NA;TminTmax = NA;toRaster=F; disturbanceON = NA; ingrowth = F; clcut = 1

  #if(!fmi_from_allas){
    rcps0 = "CurrClim"; harvScen="Base"; harvInten="Base"; forceSaveInitSoil=T
    out <- runModelAdapt(1,sampleID=1, outType = outType, rcps = rcps0, 
                         harvScen="Base", 
                         harvInten="Base", forceSaveInitSoil=T)
  #} else {
  #  out_fmi <- runModelAdapt(1,sampleID=1, outType = outType, rcps = "CurrClim_fmi", harvScen="Base", harvInten="Base",forceSaveInitSoil=T)
  #}
  #  sampleIDs <- 1:length(ops)
  #  lapply(sampleIDs, 
  #         function(jx) { 
  #           runModelAdapt(1,sampleID=jx, outType = outType, rcps = "CurrClim_fmi",
  #                         harvScen="Base", harvInten="Base",
  #                         forceSaveInitSoil=T)
  #           gc()
  #         })
  print(paste("SBB area",sum(ops[[1]]$area[which(ops[[1]]$forestdamagequalifier=="1602")])))
  if(fmi_from_allas) rcps <- "CurrClim_fmi"
  nYears <<- 2024-2015
  ops <<- list(dataS)#list(dataS[,..cols_in_prebas])
  endingYear <<- nYears + startingYear
  clcuts <<- 1
  disturbanceON <- c("fire","wind","bb")
  source("~/adaptFirst_runs/functions.R", local=T)
  sampleXs <- runModelAdapt(1,sampleID=1, outType = outType, rcps = rcps, 
                            disturbanceON = disturbanceON, 
                            harvScen="NoHarv", harvInten="NoHarv")
  #sampleXs0 <- lapply(sampleIDs, 
  #                    function(jx) { 
  #                      runModelAdapt(1,sampleID = jx, outType=outType,  
  #                                    harvScen="NoHarv",rcps = rcps, # clcut = clcuts,
  #                                    harvInten="NoHarv")
  #                      #ingrowth = T, 
  #                      #clcut = -1, disturbanceON = disturbanceON)
  #                    })
    
  ops <- list(dataS)
  print(paste("SBB area",sum(ops[[1]]$area[which(ops[[1]]$forestdamagequalifier=="1602")])))
  print("SMIs:")
  print(sampleXs$region$multiOut[1,,"NEP/SMI[layer_1]",1,2])
  print("BBprob:")
  print(sampleXs$region$multiOut[1,,"Rh/SBBpob[layer_1]",1,2])
  print("fireprob:")
  print(sampleXs$region$multiOut[1,,"W_wsap/fireRisk[layer_1]",1,2])
  print("windprob:")
  print(sampleXs$region$outDist[1,,"wrisk"])

  # SBB simulated damage segments
  #SBBReactionBA <-  apply(sampleXs$region$multiOut[,,"grossGrowth/bb BA disturbed",,2],1:2,sum)
  SBBReactionBA <-  apply(sampleXs$region$multiOut[,,43,,2],1:2,sum)
  any(SBBReactionBA>0)
  BA <- apply(sampleXs$region$multiOut[,,"BA",,1],1:2,sum)
  Vrw <- apply(sampleXs$region$multiOut[,,"VroundWood",,1],1:2,sum)[,-1]
  Vrw <- cbind(Vrw,Vrw[,ncol(Vrw)])
  
  ## wind
  # all segment areas as initial values for the damages
  areaSamplew <- array(ops[[1]]$area,c(dim(SBBReactionBA))) # Segment areas where damage happened
  areaSamplew[sampleXs$region$outDist[,,"damvol"]==0] <- 0 # if no damage happened, damaged area = 0
  #areaSamplebb[SBBReactionBA>0 & Vrw==0] <- areaSamplebb[SBBReactionBA>0 & Vrw==0]*
  #  SBBReactionBA[SBBReactionBA>0 & Vrw==0]/BA[SBBReactionBA>0 & Vrw==0] # if no clearcut, only part of area damaged
  areaSamplew[BA==0] <- 0 # if no basal area, no damaged area

  ## bb  
  # all segment areas as initial values for the damages
  areaSamplebb <- array(ops[[1]]$area,c(dim(SBBReactionBA))) # Segment areas where damage happened
  areaSamplebb[SBBReactionBA==0] <- 0 # if no damage happened, damaged area = 0
  #areaSamplebb[SBBReactionBA>0 & Vrw==0] <- areaSamplebb[SBBReactionBA>0 & Vrw==0]*
  #  SBBReactionBA[SBBReactionBA>0 & Vrw==0]/BA[SBBReactionBA>0 & Vrw==0] # if no clearcut, only part of area damaged
  areaSamplebb[BA==0] <- 0 # if no basal area, no damaged area
  ##
    
  years <- (startingYear+1):endingYear
  years <- years[years>2018 & years<2024]
  bb_dam_area <- array(0,c(length(years),3),dimnames = list(years,c("sampledata","expected_sim","simulation")))
  w_dam_area <- array(0,c(length(years),3),dimnames = list(years,c("sampledata","expected_sim","simulation")))
  probs_segm <- array(0,c(length(years),6),dimnames = list(years,c("min_pw_decl_segm","pw_median_decl","pw_median_decl/median_all",
                                                                   "min_pbb_decl_segm","pbb_median_decl","pbb_median_decl/median_all")))
  
  ti <- 1
  for(ti in 1:length(years)){
    yeari <- years[ti]
    nibb <- which(ops[[1]]$forestdamagequalifier=="1602" & as.numeric(ops[[1]]$dam_year)==yeari)
    niw <- which(ops[[1]]$forestdamagequalifier=="1504" & as.numeric(ops[[1]]$dam_year)==yeari)
    pw <- sampleXs$region$outDist[,yeari-2015,"wrisk"]
    pbb <- sampleXs$region$multiOut[,yeari-2015,"Rh/SBBpob[layer_1]",1,2]
    #hist(pbb,50)
    if(length(niw)>0){
    probs_segm[ti,1:3] <- c(min(pw[niw]),median(pw[niw]),
                         median(pw[niw])/median(pw[setdiff(1:nSegs,niw)]))}
    if(length(nibb)>0){
    probs_segm[ti,4:6] <- c(min(pbb[nibb]),median(pbb[nibb]),
                         median(pbb[nibb])/median(pbb[setdiff(1:nSegs,nibb)]))}
    bb_dam_area[ti,] <- c(sum(ops[[1]]$area[which(ops[[1]]$forestdamagequalifier=="1602" &
                                                    as.numeric(ops[[1]]$dam_year)==yeari)]),
                          sum(sampleXs$region$multiOut[,yeari-2015,"Rh/SBBpob[layer_1]",1,2]*ops[[1]]$area),
                          sum(areaSamplebb[,yeari-2015]))
    w_dam_area[ti,] <- c(sum(ops[[1]]$area[which(ops[[1]]$forestdamagequalifier=="1504" &
                                                    as.numeric(ops[[1]]$dam_year)==yeari)]),
                          sum(sampleXs$region$outDist[1,yeari-2015,"wrisk"]*ops[[1]]$area),
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
  out <- list(bb_dam_area,w_dam_area,probs_segm,regnames[r_noi])
  return(out)
  rm(list=setdiff(ls(),c(toMem,"out")))
  gc()
}

#calculateStatistics(1, fmi_from_allas = fmi_from_allas, weighted = F)
output_stats <- lapply(1:length(rids), function(jx) {
  #     print(paste0("region list: ",which(rids==20),"/",length(rids)))
  calculateStatistics(jx, fmi_from_allas = fmi_from_allas, weighted = F)
})      

save(output_stats, file = "/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Results/validation_stats.rdata")
#any(sampleXs$regio$multiOut[,,"grossGrowth/bb BA disturbed",,2]>0)

