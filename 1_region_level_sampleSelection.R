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
calculateOPSdata  <-  function(r_noi, neighborIDs=T){
  toMem <- ls()
  r_no <- rnos[r_noi]
  
  print(paste("region",r_no))
  print(paste("Neighbor information =",neighborIDs))
  
  fname <- paste0("DeclaredDamages_",dam_names[inds],"_rno",r_no,"_",regnames[r_noi],".rdata")
  load(file=paste0(savepath,"/",fname))
  print(paste("File",fname,"opened"))
  #  print(head(XYdamages))
  
  FIGU <- T
  if(FIGU){
    brks <- (min(as.numeric(XYdamages$dam_year),na.rm = T)-.5):
      (max(as.numeric(XYdamages$dam_year),na.rm = T)+.5)
    par(mfrow = c(2, 2))
    a<-hist(as.numeric(XYdamages$dam_year),
            breaks=brks,plot=F)
    barplot(a$counts*16*16/100/100,names.arg = a$mids,main=paste(regnames[r_noi],"/",dam_names[inds]),
            xlab="year",ylab="ha")
    
    b<-hist(as.numeric(XYdamages$dam_year[#XYdamages$declartionmaintreespecies=="2" &
      # XYdamages$cuttingpurpose == "6" &
      XYdamages$forestdamagequalifier=="1602"]),
      breaks=brks,plot=F)
    barplot(100*b$counts/a$counts,names.arg = a$mids,main="% of decl / SBB",
            xlab="year",ylab="% of declarations",new=F,col="red")
    
    c<-hist(as.numeric(XYdamages$dam_year[#XYdamages$declartionmaintreespecies=="2" &
      # XYdamages$cuttingpurpose == "6" &
      XYdamages$forestdamagequalifier==dam_indexs[dam_names=="wind"]]),
      breaks=brks,plot=F)
    barplot(100*c$counts/a$counts,names.arg = a$mids,main="% of decl / wind",
            xlab="year",ylab="% of declarations",new=F,col="red")
    
    d<-hist(as.numeric(XYdamages$dam_year[#XYdamages$declartionmaintreespecies=="2" &
      # XYdamages$cuttingpurpose == "6" &
      XYdamages$forestdamagequalifier==dam_indexs[dam_names=="fire"]]),
      breaks=brks,plot=F)
    barplot(100*d$counts/a$counts,names.arg = a$mids,main="% of decl/ fire",
            xlab="year",ylab="% of declarations",new=F,col="red")
  }
  #nSBB <- sum(b$counts[b$mids>=2015])
  tt <- XYdamages[dam_year<2024,]   
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
  
  noPrebas <- F  
  #    devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/settings.r")
  source("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Settings_IBCCarbon.R", local=TRUE)  
  
  IDsUniq <- unique(XYdamages[,c("dam_id","segID")])
  IDsUniq[,damSegID:=1:nrow(IDsUniq)]
  
  # Give segments that are split with declaration polygons, a unique name
  XYdamages[,damSegID := match(paste(XYdamages$dam_id,XYdamages$segID),
                               paste(IDsUniq$dam_id,IDsUniq$segID))]
  areas <- XYdamages[, .N, by = list(damSegID)]
  
  tt <-XYdamages[XYdamages[,.I[which.max(spruce)], by=damSegID]$V1] # one row for each segment
  tt[,area := areas[match(areas$damSegID, tt$damSegID),N]*16^2/100^2]
  
  tt <- tt[tt$dam_year>=2019,]
  tt <- tt[segID%in%data.all$segID,] # some damage polygons are outside the landclasses 1:2
  
  #########################
  cuttinginpractise <- c(5, 8, 16, 17, 19, 21, 22, 24)
  if(FIGU){
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
    save(damInfo,file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/damInfo_",r_noi,".rdata"))
  }
  
  ########################
  nt <- 1:nrow(tt)
  # Pick samples of disturbances in the sample
  ni <- which(tt$forestdamagequalifier[nt]%in%c("1602","1504","1503"))
  print(paste("dist. all n =",length(ni)))
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
  
  nt <- nt[ni]
  
  ni <- c(nt,sample(setdiff(1:nrow(tt),ni),nSegs2-min(nSegs2,length(ni))))
  
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
  
  #  load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
  #  data.IDs$segID <- data.IDs$maakuntaID
  #  data.IDs <- data.IDs[segID!=0]
  #  setkey(data.IDs,segID)
  #  gc()
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
  
  Ntot <- sum(data.all$area)/(0.16^2)
  gc()
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
  if(neighborIDs){
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
  save(dataS,file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/MKIdata/samples_",r_noi,".rdata"))
  return(dataS)
  rm(list=setdiff(ls(),toMem))
  gc()
  
}
if(!outputPlot){
  if(asParallel){
    library("parallel")
    ops <- mclapply(rids, function(jx) {
      #      print(paste0("region list: ",which(rids==20),"/",length(rids)))
      calculateOPSdata(jx)
    }, mc.cores = 5,mc.silent=FALSE)      
  } else {
    ops <- lapply(rids, function(jx) {
      #     print(paste0("region list: ",which(rids==20),"/",length(rids)))
      calculateOPSdata(jx, neighborIDs = neighborIDs)
    })      
  }
  names(ops) <- rids
  if(toFile){ 
    save(ops,file="../PREBAStesting/MKIdata/samples.rdata")
    print("sampleSets saved")
  }
}

if(FALSE){
  load(paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/damInfo_1.rdata"))
  outputStats <- array(0,c(dim(damInfo)[1],dim(damInfo)[2],length(rnos)+1),
                       dimnames = list(rownames(damInfo),
                                       colnames(damInfo),c(regnames,"wholeCountry")))
  for(ij in 1:length(c(1,3:length(rnos)))){
    r_noi <- c(1,3:length(rnos))[ij]
    load(paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/damInfo_",r_noi,".rdata"))
    outputStats[,,r_noi] <- damInfo
  }
  ij <- ij+1
  r_noi <- r_noi+1
  outputStats[,,r_noi] <- apply(outputStats,1:2,sum)
  save(outputStats,file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/outputStats.rdata"))
}

