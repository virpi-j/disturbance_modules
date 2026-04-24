neighborsAll <- function(ij, dataS, declData, dataSS, clctDist = 150, KUVA = F){
  if(KUVA) print(paste0("rno",r_noi,": id ",ij))
 # print(paste0("rno",r_noi,": id ",ij))
  sysT <- Sys.time()
  if(!exists("clctDist")) clctDist <- 150
  if(ij==1) print(paste("clctDist =",clctDist))
  if((ij)%%min(1000,round(nSegs/4))==0){
    timeT2 <- Sys.time()
    tperiter<- timeT2-timeT
    print(paste0("neighborinfo rno",r_noi,": ",ij,"/",nSegs,
                 ", time ",round(tperiter,2)))
  }
  damSegIDij <- dataS$damSegID[ij]
  # forestdamagequalifier in declarations: a number or NA, in sample 0
  iks <- which(dataSS$damSegID==damSegIDij)
  dam_yeari <- as.numeric(dataSS$dam_year[iks])[1]
  xi <- unlist(dataSS$x[iks])
  yi <- unlist(dataSS$y[iks])
  damidi <- as.character(unlist(dataSS$dam_id[iks]))[1]
  if(length(iks)==0) print(paste("error",ij))
  # select pixels within the given range from the segment i
    
  ntmpi <- which(declData$yyd<=(max(yi)+clctDist) & 
                   declData$yyd>=(min(yi)-clctDist) & 
                   declData$xxd<=(max(xi)+clctDist) & 
                   declData$xxd>=(min(xi)-clctDist))
  # Remove pixels from the same decl as segment i
  if(length(ntmpi)>0) ntmpi <- ntmpi[which(!paste0(declData$xxd[ntmpi],"_",declData$yyd[ntmpi])%in%paste0(xi,"_",yi) &
                                             !declData$dam_idd[ntmpi]%in%damidi)]
  
  n_neigh <- length(ntmpi)
  #gc()
  #print(n_neigh)
  dclct <- dclct_south <- dSBB <- dWind <- c(1e12,NA,0)
  if(length(ntmpi)>0){
    declData <- declData[ntmpi,]
    gc()
    if(KUVA){
      xm <- mean(xi)
      ym <- mean(yi)
      par(mfrow=c(2,2))
      plot(c(declData$xxd,xi)-xm,c(declData$yyd,yi)-ym,
           xlab="x",ylab="y",asp=1,col="black")
      points(declData$xxd[which(declData$damCCutInt==1)]-xm,
             declData$yyd[which(declData$damCCutInt==1)]-ym,col="blue",pch=19)
      points(xi-xm,yi-ym,pch=19,col="green",cex=2)
    }
    n_neigh <- length(ntmpi)
    ijsd <- NULL  
    if(is.na(dam_yeari)){
      print(paste("id",ij,"NA dam year"))
    } else {
      distances <- 1e12
      if(n_neigh>0){
        # rows for all neighbors, cols for all pixels in segm
        distances <- array(0,c(n_neigh,length(xi)),
                           dimnames = list(paste0("neighpix_id",1:n_neigh),
                                           paste0("segment_i_id",1:length(xi))))
        for(ijx in 1:length(xi)){
          distances[,ijx] <- rowSums(cbind(declData$xxd-xi[ijx],
                                           declData$yyd-yi[ijx])**2)
        }
        distances[distances==0] <- 1e12 # not in the segment ij
      }
      ijsd <- NULL
      # choose all the cells close enough and from the suitable time period
      if(min(distances)<= clctDist**2){ ## any cells near enough
        damYInt <- as.integer(declData$dam_yeard %in% 
                                c((dam_yeari-5):(dam_yeari-1)))
        damYInt[damYInt==0] <- 1e12
        damYInt <- array(rep(damYInt,length(xi)), dim=c(n_neigh,length(xi)))
        mins <- apply(damYInt*distances,2,min) # for each segm.pixel, min dist to neighbor
        if(min(mins)<= clctDist**2){
          ijsd <- which(apply(damYInt*distances,1,min)<clctDist**2) # which decl data is closest
        }
        declData <- declData[ijsd,]
        #gc()
        if(KUVA){
          print(paste("segm year",dam_yeari[1]))
          print("neighbor years")
          print(unique(declData$dam_yeard))
        }
        damYInt <- damYInt[ijsd,]
        distances <- distances[ijsd,]
        n_neigh <- length(ijsd)
        #gc()
        if(KUVA & n_neigh>0){
          points(declData$xxd-xm,declData$yyd-ym,pch=19,col="red")
        }
      }
      ijs <- ijsSBB <- ijsWind <- NULL
      
      # decl segments close enough
      if(length(ijsd)>0){ 
        
        # clearcuts
        ntmp <- array(rep(declData$damCCutInt,length(xi)), dim=c(n_neigh,length(xi)))*distances
        # closest neigh with cc for each pixel in segm0
        if(n_neigh>1 & length(xi)>1){ # more than one neigh and more than one segm pixel
          mins <- apply(ntmp,2,min)
        } else {
          mins <- ntmp 
        }
        nclct <- length(which(mins<(16^2+16^2))) # how many segm pixels right next to clct pixel
        if(min(mins)<= clctDist**2){
          # find the closest, recent clct pixel 
          if(length(xi)>1 & n_neigh>1){
            ijs <- unique(apply(ntmp,2,which.min)[which(mins==min(mins))])
            ijs <- ijs[which.min(dam_yeari-declData$dam_yeard[ijs])]
          } else if(length(xi)==1) {
            ijs <- which(mins==min(mins))
            ijs <- ijs[which.min(dam_yeari-declData$dam_yeard[ijs])]
          } else {
            ijs <- 1 
          }
          dd <- sqrt(min(ntmp[ijs,])) 
          dyr <- declData$dam_yeard[ijs]-dam_yeari
          dclct <- c(dd,dyr,nclct)
          if(is.na(dyr)&dd>0){
            print(ij)
            eaeta <- eaatg
          }
          if(KUVA){
            plot(c(declData$xxd,xi)-xm,c(declData$yyd,yi)-ym,
                 xlab="x",ylab="y",asp=1,col="black", 
                 main=paste0("neigh clct",ij," d=",round(dd,2),", yr=",dyr,", n=",nclct))
            points(declData$xxd[which(declData$damCCutInt==1)]-xm,
                   declData$yyd[which(declData$damCCutInt==1)]-ym,col="blue",pch=19)
            points(xi-xm,yi-ym,pch=19,col="green",cex=2)
            points(declData$xxd[ijs]-xm,declData$yyd[ijs]-ym,col="red",pch=4,cex=2)
          }
        }
        
        # clct in south
        damSouth <- damYInt*0
        for(icol in 1:length(xi)){
          damSouth[which(abs(declData$xxd-xi[icol])<sqrt((2*16)**2) & 
                           declData$yyd<yi[icol])]<-1
        }
        damSouth[damSouth==0] <- 1e12
        damSouth <- array(rep(damSouth,length(xi)), dim=c(n_neigh,length(xi)))
        damCCutInt <- array(rep(declData$damCCutInt,length(xi)), dim=c(n_neigh,length(xi)))
        
        ntmp <- damCCutInt*damYInt*distances*damSouth
        #gc()
        if(length(ijsd)>1 & length(xi)>1){
          mins <- apply(ntmp,2,min)
        } else {
          mins <- ntmp 
        }
        nclct_south <- length(which(mins<(16^2+16^2)))
        if(min(mins)<= clctDist**2){
          if(length(xi)>1 & n_neigh>1){
            ijs <- unique(apply(ntmp,2,which.min)[which(mins==min(mins))])
            ijs <- ijs[which.min(dam_yeari-declData$dam_yeard[ijs])]
          } else if(length(xi)==1) {
            ijs <- which(mins==min(mins))
            ijs <- ijs[which.min(dam_yeari-declData$dam_yeard[ijs])]
          } else {
            ijs <- 1
          }
          dd <- sqrt(min(ntmp[ijs,])) 
#          dd <- sqrt(min(mins)) 
          dyr <- declData$dam_yeard[ijs]-dam_yeari
          dclct_south <- c(dd,dyr,nclct_south)
          if(KUVA){
            plot(c(declData$xxd,xi)-xm,c(declData$yyd,yi)-ym,
                 xlab="x",ylab="y",col="black",asp=1, main=paste("neigh south d=",round(dd,2),"n=",nclct_south))
            points(declData$xxd[which(declData$damCCutInt==1)]-xm,
                   declData$yyd[which(declData$damCCutInt==1)]-ym,col="blue",pch=19)
            points(xi-xm,yi-ym,pch=19,col="green",cex=2)
            points(declData$xxd[ijs]-xm,declData$yyd[ijs]-ym,col="red",pch=4,cex=2)
          }
        }
        
        # SBB in neighborhood
        
        damBBInt <- array(rep(declData$damBBInt,length(xi)), dim=c(n_neigh,length(xi)))
        ntmp <- damBBInt*damYInt*distances
        #gc()
        if(length(ijsd)>1 & length(xi)>1){
          mins <- apply(ntmp,2,min)
        } else {
          mins <- ntmp 
        }
        nbb <- length(which(mins<(16^2+16^2)))
        if(min(mins)<= clctDist**2){
          if(length(xi)>1 & n_neigh>1){
            ijs <- unique(apply(ntmp,2,which.min)[which(mins==min(mins))])
            ijs <- ijs[which.min(dam_yeari-declData$dam_yeard[ijs])]
          } else if(length(xi)==1) {
            ijs <- which(mins==min(mins))
            ijs <- ijs[which.min(dam_yeari-declData$dam_yeard[ijs])]
          } else {
            ijs <- 1 
          }
          dd <- sqrt(min(ntmp[ijs,])) 
          #dd <- sqrt(min(mins)) 
          dyr <- declData$dam_yeard[ijs]-dam_yeari
          dSBB <- c(dd,dyr,nbb)
          if(KUVA){
            plot(c(declData$xxd,xi)-xm,c(declData$yyd,yi)-ym,
                 xlab="x",ylab="y",col="black",asp=1, 
                 main=paste("neigh bb d=",round(dd,2),"n=",nbb))
            points(declData$xxd[which(declData$damBBInt==1)]-xm,
                   declData$yyd[which(declData$damBBInt==1)]-ym,col="blue",pch=19)
            points(xi-xm,yi-ym,pch=19,col="green",cex=2)
            points(declData$xxd[ijs]-xm,declData$yyd[ijs]-ym,col="red",pch=4,cex=2)
          }
        }
        
        
        # Wind damage in neighborhood
        damWInt <- array(rep(declData$damWInt,length(xi)), dim=c(n_neigh,length(xi)))
        ntmp <- damWInt*damYInt*distances
        #gc()
        if(length(ijsd)>1 & length(xi)>1){
          mins <- apply(ntmp,2,min)
        } else {
          mins <- ntmp 
        }
        nwind <- length(which(mins<(16^2+16^2)))
        if(is.na(nwind)) break()
        if(min(mins)<= clctDist**2){
          if(length(xi)>1 & n_neigh>1){
            ijs <- unique(apply(ntmp,2,which.min)[which(mins==min(mins))])
            ijs <- ijs[which.min(dam_yeari-declData$dam_yeard[ijs])]
          } else if(length(xi)==1) {
            ijs <- which(mins==min(mins))
            ijs <- ijs[which.min(dam_yeari-declData$dam_yeard[ijs])]
          } else {
            ijs <- ijsd 
          }
          dd <- sqrt(min(mins)) 
          dyr <- declData$dam_yeard[ijs]-dam_yeari
          dWind <- c(dd,dyr,nwind)
          if(KUVA){
            plot(c(declData$xxd,xi)-xm,c(declData$yyd,yi)-ym,
                 xlab="x",ylab="y",col="black",asp=1, 
                 main=paste("neigh wind d=",round(dd,2),"n=",nwind))
            points(declData$xxd[which(declData$damWInt==1)]-xm,
                   declData$yyd[which(declData$damWInt==1)]-ym,col="blue",pch=19)
            points(xi-xm,yi-ym,pch=19,col="green",cex=2)
            points(declData$xxd[ijs]-xm,declData$yyd[ijs]-ym,col="red",pch=4,cex=2)
          }
        }
      }
    }
    
  }    
    
    ############################################################
#  if(exists("distances")) rm("distances")
#  if(exists("ntmp")) rm("ntmp")
#gc()  
  out <- c(dclct, dSBB, dWind, dclct_south)
  #print(paste(ij,Sys.time()-sysT))
  #print(out)
  return(out)
}

