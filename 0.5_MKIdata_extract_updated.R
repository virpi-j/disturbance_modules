rm(list=ls())
gc()
setwd("/scratch/project_2000994/PREBASruns/PREBAStesting/")
library(sf)
library(data.table)
library("readxl")
library("raster")
library("terra")
toFile <- T


rnames <- c("Uusimaa","Ahvenanmaa","Keski-Pohjanmaa","Pirkanmaa","Etel%C3%A4-Karjala","Keski-Suomi",
            "Pohjois-Savo","Lappi_E","Lappi_P","Kanta-H%C3%A4me","Pohjanmaa","Varsinais-Suomi",
            "Etel%C3%A4-Pohjanmaa","P%C3%A4ij%C3%A4t-H%C3%A4me","Satakunta","Kymenlaakso",
            "Kainuu","Etel%C3%A4-Savo","Pohjois-Karjala","Pohjois-Pohjanmaa")
regnames <- c("Uusimaa","Ahvenanmaa","Keski-Pohjanmaa","Pirkanmaa","Etela-Karjala","Keski-Suomi",
              "Pohjois-Savo","Lappi_E","Lappi_P","Kanta-Hame","Pohjanmaa","Varsinais-Suomi",
              "Etela-Pohjanmaa","Paijat-Hame","Satakunta","Kymenlaakso",
              "Kainuu","Etela-Savo","Pohjois-Karjala","Pohjois-Pohjanmaa")
r_noi <- r_no <- 1
rnos <- c(1:8,8:19)

# path to save the outputs
savepath = "/scratch/project_2000994/PREBASruns/PREBAStesting/MKIdata"

# forest damage type indexes for the disturbaces of interest
dam_names <- c("alldeclarations","all","SBB","wind","fire","moose")
dam_indexs <- c("all","0","1602","1504","1503","1650")

clearcuts <- c(5, 8, 16, 17, 19, 21, 22, 24) # cuttingrealizationpractice codes for clearcuts

regs <- c(1,3:20)
dam_years <- 2011:2024

outputStats <- array(0,c(8,1+max(dam_years)-min(dam_years),length(regs)+1),
                     dimnames = list(paste0(dam_names[c(1,1,3,3,4,4,5,5)],c("_all","_clearcut")),
                                     paste0("year",dam_years),c(regnames[regs],"WholeCountry")))

library(sp)
#library(rgdal)

upm <- cbind(c(61.068067045391004, 28.239737695831963),c(61.161953, 26.820275), c(60.967942, 26.671960),
             c(60.778085, 26.890481),c(61.12872381073421, 28.476726026705737),
             c(61.253919991113364, 28.864573071682752),c(61.43498780382168, 29.343663069179648))
colnames(upm) <- c("UPMKaukas","UPMKymmene","UPMKymi","SEAnjala","SEJoutseno","SEKaukoaa","MBSimpele")
xy <- data.table(ID=1:ncol(upm),x=upm[1,],y=upm[2,])
coordinates(xy) <- c("y","x")
#crsX <- ("+proj=utm +zone=35 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m+no_defs")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84") 
#res <- spTransform(xy, crsX)#CRS(paste("+proj=utm +zone=",35," ellps=WGS84",sep='')))
upms <- coordinates(spTransform(xy, CRS(paste("+proj=utm +zone=",35," ellps=GRS80",sep=''))))
#extent(upms)
#proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
#res <- spTransform(xy, CRS("+proj=utm +zone=51 ellps=WGS84"))
#extent(res)


DeclToRaster <- T # T if update the declaration database
if(DeclToRaster) dam_years <- 2015:2024

r_noi <- 1
#disturbance_extract <- function(r_noi){
for(r_noi in 1:length(regs)){
  r_no <- rnos[regs[r_noi]]
  toMem <- ls()
  print(paste("region",r_no))
  
  download.file(paste0("https://avoin.metsakeskus.fi/aineistot/Metsankayttoilmoitukset/Maakunta/MKI_",
                       rnames[regs[r_noi]],".zip"),destfile = "test.zip")
  unzip("test.zip")
  unlink("test.zip")
  files <- list.files()
  print(files[grepl(".gpkg", files, fixed = TRUE)])
  file.rename(from=files[grepl(".gpkg", files, fixed = TRUE)][1],to="test.gpkg")
  layers = c("forestdamagequalifier","maingroup",#"fertilityclass","soiltype",
             #"declarationmaintreespecies","meanage","meandiameter",
             "cuttingpurpose","cuttingrealizationpractice",
             "declarationarrivaldate","geometry")
  tmp<-st_read("test.gpkg", layer = "forestusedeclaration")[,layers] # read forestdeclaration layer
  files <- list.files()
  unlink(files[grepl(".gpkg", files, fixed = TRUE)]) # remove all gpkg-files
  gc()
  # pick only damage year 2016,2017,... data
  dam_year <- substr(tmp$declarationarrivaldate,1,4)
  dam <- which(as.numeric(dam_year)>=min(dam_years) & as.numeric(dam_year)<=max(dam_years) & tmp$maingroup%in%c("1","2"))
  tmp<- tmp[dam,]
  dam_year <- dam_year[dam]#substr(tmp$declarationarrivaldate[dam],1,4)
  rm(dam)
  gc()

  ti <- 1
  
  for(ti in 1:length(dam_years)){
    #tmp <- XYdamages[which(XYdamages$dam_year==dam_years[ti]),]
    #outputStats[1,ti,r_noi] <- nrow(tmp)*16*16/100/100
    #outputStats[2,ti,r_noi] <- length(which(tmp$cuttingrealizationpractice%in%clearcuts))*16*16/100/100
    outputStats[1,ti,r_noi] <- sum(st_area(tmp[which(as.numeric(dam_year)==dam_years[ti]),]))
    outputStats[2,ti,r_noi] <- sum(st_area(tmp[which(as.numeric(dam_year)==dam_years[ti] & tmp$cuttingrealizationpractice%in%clearcuts),]))
    
    #ni <- which(tmp$forestdamagequalifier==dam_indexs[3])
    #outputStats[3,ti,r_noi] <- length(ni)*16*16/100/100
    #outputStats[4,ti,r_noi] <- length(which(tmp$cuttingrealizationpractice[ni]%in%clearcuts))*16*16/100/100
    outputStats[3,ti,r_noi] <- sum(st_area(tmp[which(tmp$forestdamagequalifier==dam_indexs[3] & as.numeric(dam_year)==dam_years[ti]),]))
    outputStats[4,ti,r_noi] <- sum(st_area(tmp[which(tmp$forestdamagequalifier==dam_indexs[3] & as.numeric(dam_year)==dam_years[ti] & tmp$cuttingrealizationpractice%in%clearcuts),]))
    
    #ni <- which(tmp$forestdamagequalifier==dam_indexs[4])
    #outputStats[5,ti,r_noi] <- length(ni)*16*16/100/100
    #outputStats[6,ti,r_noi] <- length(which(tmp$cuttingrealizationpractice[ni]%in%clearcuts))*16*16/100/100
    outputStats[5,ti,r_noi] <- sum(st_area(tmp[which(tmp$forestdamagequalifier==dam_indexs[4] & as.numeric(dam_year)==dam_years[ti]),]))
    outputStats[6,ti,r_noi] <- sum(st_area(tmp[which(tmp$forestdamagequalifier==dam_indexs[4] & as.numeric(dam_year)==dam_years[ti] & tmp$cuttingrealizationpractice%in%clearcuts),]))
    
    #ni <- which(tmp$forestdamagequalifier==dam_indexs[5])
    #outputStats[7,ti,r_noi] <- length(ni)*16*16/100/100
    #outputStats[8,ti,r_noi] <- length(which(tmp$cuttingrealizationpractice[ni]%in%clearcuts))*16*16/100/100
    outputStats[7,ti,r_noi] <- sum(st_area(tmp[which(tmp$forestdamagequalifier==dam_indexs[5] & as.numeric(dam_year)==dam_years[ti]),]))
    outputStats[8,ti,r_noi] <- sum(st_area(tmp[which(tmp$forestdamagequalifier==dam_indexs[5] & as.numeric(dam_year)==dam_years[ti] & tmp$cuttingrealizationpractice%in%clearcuts),]))
    
  }
  outputStats[,,r_noi] <-   outputStats[,,r_noi]/100/100
  print(round(outputStats[,,r_noi],2))
  
  
  if(DeclToRaster){
    # read region level segments used in PREBAS runs
    #CSCrun<-T
    load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/data.all_maakunta_",r_no,".rdata"))
    data.all$segID <- data.all$maakuntaID
    data.all <- data.all[segID!=0,c("segID","age","ba","pine","spruce","birch")]
    setkey(data.all,segID)
    gc()
    
    setwd("/scratch/project_2000994/PREBASruns/PREBAStesting/")
    load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
    data.IDs$segID <- data.IDs$maakuntaID
    data.IDs <- data.IDs[segID!=0]
    setkey(data.IDs,segID)
    gc()
    
    tabX <- merge(data.IDs,data.all) # coords of the segments
    rm(list=c("data.IDs","data.all"))
    gc()
    colnames(tabX)[colnames(tabX)=="x.x"]<-"x"
    colnames(tabX)[colnames(tabX)=="y.x"]<-"y"
    crsX <- ("+proj=utm +zone=35 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    bbox<-st_bbox(tmp)
    tabX<-tabX[(tabX$x>=bbox[1] & tabX$x<=bbox[3] & tabX$y>=bbox[2] &tabX$y <=bbox[4]),]
    gc()
    pts <- vect(cbind(tabX$x,tabX$y),crs = crsX)
    #values(pts) <- data.table(segID=tabX$segID, id=1:nrow(tabX),x=tabX$x,y=tabX$y)
    values(pts) <- tabX
    rm(tabX)
    gc()
    
    
    # pick only damage year 2016,2017,... data
    #dam_year <- substr(tmp$declarationarrivaldate,1,4)
    #dam <- dam[as.numeric(dam_year)>=2016-5]
    #dam_year <- dam_year[dam]#substr(tmp$declarationarrivaldate[dam],1,4)
    v <- vect(tmp)
    #v <- st_as_sf(tmp)
    ntmp <- ncol(tmp)
    values(v)[ntmp-1] <- dam_year
    names(v)[ntmp-1] <- "dam_year"
    rm(list=c("tmp","dam_year"))
    gc()
    crs(v) <- crs(pts)
    values(v)[ntmp] <- data.table(dam_id = 1:nrow(v))
    gc()
    rm(list="xy")
    gc()
    
    #XYdamages <- as.data.table(intersect(v[which(values(v)[,"dam_year"]=="2019"),], pts))
    #XYdamages <- as.data.table(intersect(v[which(values(v)[,"dam_year"]%in%c(2016:2024)),], pts))
    XYdamages <- as.data.table(intersect(v, pts))
    print(head(XYdamages))
    #XYdamages <- XYdamages[which(XYdamages$maingroup%in%c("1","2")),] # forest land and poorly productive forest land only
    
    rm(list=c("v","pts"))
    gc()
    
    #brks <- (min(as.numeric(XYdamages$dam_year),na.rm = T)-.5):
    #  (max(as.numeric(XYdamages$dam_year),na.rm = T)+.5)
    #par(mfrow = c(1,2))
    #a<-hist(as.numeric(XYdamages$dam_year),
    #        breaks=brks,plot=F)
    #ylims <- c(0,max(a$counts*16*16/100/100))
    #barplot(a$counts*16*16/100/100,names.arg = a$mids,main=paste(regnames[r_noi],"/",dam_names[1]),
    #            xlab="year",ylab="ha", ylim=ylims)
    
    #a<-hist(as.numeric(XYdamages$dam_year[which(XYdamages$cuttingrealizationpractice%in%clearcuts)]),
    #        breaks=brks,plot=F)
    #barplot(a$counts*16*16/100/100,names.arg = a$mids,main=paste(regnames[r_noi],"/ clearcuts / ",dam_names[1]),
    #        xlab="year",ylab="ha", ylim=ylims)
    
    
    #inds <- 3      
    byDamages <- F
    if(byDamages){
      for(inds in 3:5){ # go though these forest damage type indexes from dam_indexs 
        toMem2 <- ls()
        if(dam_indexs[inds]=="0"){ # all damages
          dam<-which(XYdamages$forestdamagequalifier>=1500) 
        } else { # pick only the damage type of interest
          dam<-which(XYdamages$forestdamagequalifier==dam_indexs[inds]) 
        }
        
        if(length(dam)>0){ # if there is any damage data of the given type, draw and write in memory
          #tmp <- XYdamages[dam,]
          par(mfrow = c(1,2))
          a<-hist(as.numeric(XYdamages$dam_year[dam]),
                  breaks=brks,plot=F)
          ylims <- c(0,max(a$counts*16*16/100/100))
          barplot(a$counts*16*16/100/100,names.arg = a$mids,main=paste(regnames[r_noi],"/",dam_names[inds]),
                  xlab="year",ylab="ha", ylim=ylims)
          
          a<-hist(as.numeric(XYdamages$dam_year[dam[which(XYdamages$cuttingrealizationpractice%in%clearcuts)]]),
                  breaks=brks,plot=F)
          barplot(a$counts*16*16/100/100,names.arg = a$mids,main=paste(regnames[r_noi],"/ clearcuts / ",dam_names[inds]),
                  xlab="year",ylab="ha", ylim=ylims)
          
          if(inds==30){
            tt <- XYdamages      
            ttsum = tt$pine+tt$spruce+tt$birch
            tt[,pine := pine/ttsum]
            tt[,spruce := spruce/ttsum]
            tt[,birch:= birch/ttsum]
            tt<-XYdamages[XYdamages[, .I[which.max(spruce)], by=dam_id]$V1]
            tt2 <- tt[tt$forestdamagequalifier=="1602" & 
                        tt$declarationmaintreespecies=="2" &
                        tt$cuttingrealizationpractice%in%clearcuts,]
            PAIRS <- F
            if(nrow(tt2)>0 & PAIRS){
              #pairs(tt2[,c("age","spruce","ba")],na.rm=T,
              #       upper.panel = panel.cor,
              ##       diag.panel  = panel.hist,
              #       lower.panel = panel.smooth,
              #       pch = ".",
              #       main="SBB pixels")
              par(mfrow = c(2, 2))
              for(vnams in names(tt)[c(10,17,18,20)]){
                a <- hist(as.numeric(as.matrix(tt[,..vnams])),plot=F)
                b <- hist(as.numeric(as.matrix(tt2[,..vnams])),breaks=a$breaks,plot=F)
                barplot(b$counts/a$counts,names.arg = a$mids,
                        main=paste(regnames[r_noi],"/","SBB from all decl"),
                        xlab=vnams,ylab="% of declarations",new=F,col="black")
              }
            }
          }
          
          if(toFile){
            
            newRes <- F
            if(newRes){ # If update the declaration data used in PREBAS
              fname <- paste0("DeclaredDamages_",dam_names[inds],"_rno",r_no,"_",regnames[r_noi],".rdata")
              save(XYdamages,file=paste0(savepath,"/",fname))
              print(paste("File",fname))
              print(paste("saved to folder",savepath))
            }
          }
        }
        rm(list=setdiff(ls(),toMem2))
        gc()
      }
    }}
  #if(toFile) dev.off()
  
  #  rm(tabX);
  rm(list=setdiff(ls(),toMem))
  gc()
  #return(outputStats)
}

#library(parallel)
if(!toFile){
  regs <- c(1,3:5)
#} else {
#  library(parallel)
#  regs <- c(1,3:20)
}
#regs <- 6
#regs <- c(19,20)
#out <- lapply(regs, function(jx) {
#  disturbance_extract(jx)
#  gc()
#})    

#names(out) <- regnames[regs]
#  out <- mclapply(regs, function(jx) {
#    disturbance_extract(jx)
#  }, mc.cores = 4,mc.silent=FALSE)
#}
outputStats[,,dim(outputStats)[3]] <- apply(outputStats,1:2,sum)
save(outputStats, file = paste0(savepath,"/outputStats.rdata"))
outputStatsAll <- apply(outputStats,1:2,sum)

if(dev.interactive()) dev.off()
pdf(paste0(savepath,"outputStats.pdf"))
if(!DeclToRaster){
  par(mfrow=c(4,2))
  for(r_noi in 1:dim(outputStats)[3]){
    ylims <- array(0,c(dim(outputStats)[1],2))
    ylims[,2] <- apply(outputStats[,,r_noi],1,max)
    for(ii in 1:dim(outputStats)[1]){
      #brks <- (min(as.numeric(XYdamages$dam_year),na.rm = T)-.5):
      #  (max(as.numeric(XYdamages$dam_year),na.rm = T)+.5)
      #par(mfrow = c(1,2))
      #a<-hist(as.numeric(outputStats[ii,,r_noi]),
      #        breaks=brks,plot=F)
      #ylims <- c(0,max(a$counts*16*16/100/100))
      if(ii%%2!=0){
        barplot(outputStats[ii,-1,r_noi],names.arg = dam_years[-1],
                main=paste(dimnames(outputStats)[[3]][r_noi],"/",dimnames(outputStats)[[1]][ii]),
                  xlab="year",ylab="ha", ylim=ylims[ii,])
      } else {
        hf <- barplot(outputStats[ii,-1,r_noi]/outputStats[ii-1,-1,r_noi]*100,
                names.arg = dam_years[-1],
                main=paste(dimnames(outputStats)[[1]][ii],"/ FI ave ",
                           round(sum(outputStatsAll[ii,-1])/sum(outputStatsAll[ii-1,-1])*100,1),"%"),
                xlab="year",ylab="%", ylim=c(0,105))
        if(sum(outputStats[ii-1,,r_noi])>0){
          lines(c(hf[1],hf[length(hf)]),c(1,1)*sum(outputStats[ii,-1,r_noi])/sum(outputStats[ii-1,-1,r_noi])*100,
                col="black")
        }
          lines(c(hf[1],hf[length(hf)]),c(1,1)*sum(outputStatsAll[ii,-1])/sum(outputStatsAll[ii-1,-1])*100,
              col="red")
        
      }
    }
  }
}
dev.off()
