neighborsAll <- function(ij, dataS, declData, clctDist = 16*70, KUVA = F){
  #print(paste0("rno",r_noi,": id ",ij))
  if(!exists("clctDist")) clctDist <- 16*70
  if((ij)%%min(1000,round(nSegs/4))==0){
    timeT2 <- Sys.time()
    tperiter<- timeT2-timeT
    print(paste0("neighborinfo rno",r_noi,": ",ij,"/",nSegs,
                 ", time ",round(tperiter,2)))
  }
  damSegIDij <- dataS$damSegID[ij]
  if(is.na(dataS$forestdamagequalifier[ij]) | dataS$forestdamagequalifier[ij]>0){ # in declarations
    iks <- which(XYdamages$damSegID==damSegIDij)
    dam_yeari <- as.numeric(XYdamages$dam_year[iks[1]]) # damage year of the dam_id
    xi <- XYdamages$x[iks]
    yi <- XYdamages$y[iks]
    damidi <- XYdamages$dam_id[iks]
  } else if(dataS$forestdamagequalifier[ij]==0){          # in sample
    iks <- which(dataS$damSegID==damSegIDij) # all segments in the dam_id
    dam_yeari <- as.numeric(dataS$dam_year[iks[1]]) # damage year of the dam_id
    xi <- dataS$x[iks]
    yi <- dataS$y[iks]
    damidi <- dataS$dam_id[iks]
  } else { print("error")}
  #if(is.na(dam_yeari)) print(paste(ij,dam_indd[iks],dam_yeari))
  # choose declarations from years dam_yeari-4:dam_yeari-1 and thus also not in declaration ij
   # possible range to be checked clctDist
  #UUSI <- T
  # select pixels within the given range from the segment i
  ntmpi <- which(declData$yyd<=(max(yi)+clctDist) & 
                   declData$yyd>=(min(yi)-clctDist) & 
                   declData$xxd<=(max(xi)+clctDist) & 
                   declData$xxd>=(min(xi)-clctDist))
  # Remove pixels from the same decl as segment i
  if(length(ntmpi)>0) ntmpi <- ntmpi[which(!paste0(declData$xxd[ntmpi],"_",declData$yyd[ntmpi])%in%paste0(xi,"_",yi) &
                                             !declData$dam_idd[ntmpi]%in%damidi)]
  
  n_neigh <- length(ntmpi)
  #print(n_neigh)
  dclct <- dclct_south <- dSBB <- dWind <- c(1e12,NA,0)
  if(length(ntmpi)>0){
    declData <- declData[ntmpi,]

    n_neigh <- length(ntmpi)
    ijsd <- NULL  
    if(is.na(dam_yeari)){
      print(paste("id",ij,"NA dam year"))
    } else {
      distances <- 1e12
      if(n_neigh>0){
        #print(ij)
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
      if(min(distances)<= clctDist**2){ ## any cells near enough
        damYInt <- as.integer(declData$dam_yeard %in% 
                                c((dam_yeari-15):(dam_yeari)))
        damYInt[damYInt==0] <- 1e12
        mins <- apply(damYInt*distances,2,min) # for each segm.pixel, min dist to neighbor
        if(min(mins)<= clctDist**2){
          ijsd <- which(apply(damYInt*distances,1,min)<clctDist**2) # which decl data is closest
        }
        declData <- declData[ijsd,]
        damYInt <- damYInt[ijsd]
        distances <- distances[ijsd,]
        n_neigh <- length(ijsd)
        if(KUVA & n_neigh>0){
          xm <- mean(xi)
          ym <- mean(yi)
          par(mfrow=c(2,2))
          plot(c(declData$xxd,xi)-xm,c(declData$yyd,yi)-ym,
               xlab="x",ylab="y",col="black",asp=1, main=paste("neigh pixs, id=",ij))
          points(xi-xm,yi-ym,pch=19,col="green")
        }
      }
      ijs <- ijsSBB <- ijsWind <- NULL
      
      # decl segments close enough
      if(length(ijsd)>0){ 
        
        # clearcuts
        ntmp <- declData$damCCutInt*damYInt*distances
        if(n_neigh>1 & length(xi)>1){
          mins <- apply(ntmp,2,min)
          #    } else if(ncol(distances)==1) {
          #      mins <- ntmp
        } else {
          mins <- ntmp 
        }
        nclct <- length(which(mins<23^2))
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
          dclct <- c(dd,dyr,nclct)
          if(KUVA){
            plot(c(declData$xxd,xi)-xm,c(declData$yyd,yi)-ym,
                 xlab="x",ylab="y",asp=1,col="black", main=paste("neigh clct d=",round(dd,2),"n=",nclct))
            points(declData$xxd[which(declData$damCCutInt==1)]-xm,
                   declData$yyd[which(declData$damCCutInt==1)]-ym,col="blue",pch=19)
            points(xi-xm,yi-ym,pch=19,col="green",cex=2)
            points(declData$xxd[ijs]-xm,declData$yyd[ijs]-ym,col="red",pch=4,cex=2)
          }
        }
        
        # clct in south
        damSouth <- damYInt*0
        for(icol in 1:length(xi)){
          damSouth[which(abs(declData$xxd-xi[icol])<32 & declData$yyd<yi[icol])]<-1
        }
        damSouth[damSouth==0] <- 1e12
        
        ntmp <- declData$damCCutInt*damYInt*(distances*damSouth)
        if(length(ijsd)>1 & length(xi)>1){
          mins <- apply(ntmp,2,min)
          #    } else if(length(xi)==1) {
          #      mins <- ntmp
        } else {
          mins <- ntmp 
        }
        nclct_south <- length(which(mins<23^2))
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
        ntmp <- declData$damBBInt*damYInt*distances
        if(length(ijsd)>1 & length(xi)>1){
          mins <- apply(ntmp,2,min)
          #    } else if(length(xi)==1) {
          #      mins <- ntmp
        } else {
          mins <- ntmp 
        }
        nbb <- length(which(mins<23^2))
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
        ntmp <- declData$damWInt*damYInt*distances
        if(length(ijsd)>1 & length(xi)>1){
          mins <- apply(ntmp,2,min)
          #    } else if(length(xi)==1) {
          #      mins <- ntmp
        } else {
          mins <- ntmp 
        }
        nwind <- length(which(mins<23^2))
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
       # if(KUVA){
      #    out <- round(c(ij, dclct, dSBB, dWind, dclct_south),2)
       #   print(out)
       # }
      }
    }
    
  }    
    
    ############################################################
  
  out <- c(dclct, dSBB, dWind, dclct_south)
  #print(out)
  return(out)
}

neighborsAll_old <- function(ij, dataS, damCCutInt, damBBInt, damWInt){
  #print(paste0("rno",r_noi,": id",ij))
  if((ij)%%min(1000,round(nSegs/4))==0){
    timeT2 <- Sys.time()
    tperiter<- timeT2-timeT
    print(paste0("neighborinfo rno",r_noi,": ",ij,"/",nSegs,
                 ", time ",round(tperiter,2)))
    #       timeT <- timeT2
  }
  damSegIDij <- dataS$damSegID[ij]
  if(is.na(dataS$forestdamagequalifier[ij]) | dataS$forestdamagequalifier[ij]>0){ # in declarations
    iks <- which(XYdamages$damSegID==damSegIDij)
    dam_yeari <- as.numeric(XYdamages$dam_year[iks[1]]) # damage year of the dam_id
    xi <- XYdamages$x[iks]
    yi <- XYdamages$y[iks]
    damidi <- XYdamages$dam_id[iks]
  } else if(dataS$forestdamagequalifier[ij]==0){          # in sample
    iks <- which(dataS$damSegID==damSegIDij) # all segments in the dam_id
    dam_yeari <- as.numeric(dataS$dam_year[iks[1]]) # damage year of the dam_id
    xi <- dataS$x[iks]
    yi <- dataS$y[iks]
    damidi <- dataS$dam_id[iks]
  } else { print("error")}
  #if(is.na(dam_yeari)) print(paste(ij,dam_indd[iks],dam_yeari))
  # choose declarations from years dam_yeari-4:dam_yeari-1 and thus also not in declaration ij
  KUVA <- F
  clctDist <- 16*7 # possible range to be checked
  #UUSI <- T
  ntmpi <- which(yyd<=(max(yi)+clctDist) & yyd>=(min(yi)-clctDist) & 
                   xxd<=(max(xi)+clctDist) & xxd>=(min(xi)-clctDist))
  if(length(ntmpi)>0) ntmpi <- ntmpi[which(!paste0(xxd[ntmpi],"_",yyd[ntmpi])%in%paste0(xi,"_",yi) &
                                             !dam_idd[ntmpi]%in%damidi)]
  ijsd <- NULL  
  if(is.na(dam_yeari)){
    ijsd <- NULL  
    print(paste("id",ij,"NA dam year"))
  } else {
    #    if(UUSI){
    distances <- 1e12
    if(length(ntmpi)>0){
      #print(ij)
      distances <- array(0,c(length(xxd[ntmpi]),length(xi)))
      #rowdist <- function(ijx) array(rowSums(cbind(xxd-max(xi[ijx]),yyd-max(yi[ijx]))**2),c(length(xxd),1))
      for(ijx in 1:length(xi)){
        distances[,ijx] <- rowSums(cbind(xxd[ntmpi]-max(xi[ijx]),yyd[ntmpi]-max(yi[ijx]))**2)
      }
      distances[distances==0] <- 1e12 # not in the segment ij
    }
    #distances <- array(0,c(length(xxd),4))
    #distances[,1] <-  rowSums(cbind(xxd-max(xi),yyd-max(yi))**2)
    #distances[,2] <-  rowSums(cbind(xxd-min(xi),yyd-max(yi))**2)
    #distances[,3] <-  rowSums(cbind(xxd-max(xi),yyd-min(yi))**2)
    #distances[,4] <-  rowSums(cbind(xxd-min(xi),yyd-min(yi))**2)
    if(min(distances)<= clctDist**2){ ## any cells near enough
      damYInt <- as.integer(dam_yeard[ntmpi] %in% c((dam_yeari-15):(dam_yeari)))
      damYInt[damYInt==0] <- 1e12
      mins <- apply(damYInt*distances,2,min)
      if(min(mins)<= clctDist**2){
        ijsd <- which(apply(damYInt*distances,1,min)<clctDist**2)
        #ijsd <- unique(apply(damYInt*distances,2,which.min)[which(mins==min(mins))])
      } #ijsd <- 1
    }
    
  }
  ijs <- ijsSBB <- ijsWind <- NULL
  if(length(ijsd)>0){ # clearcuts in neighborhood
    if(length(ijsd)>1 & ncol(distances)>1){
      ntmp <- damCCutInt[ntmpi[ijsd]]*damYInt[ijsd]*distances[ijsd,]
      mins <- apply(ntmp,2,min)
    } else if(ncol(distances)==1) {
      ntmp <- damCCutInt[ijsd]*damYInt[ijsd]*distances[ijsd,1]
      mins <- ntmp
    } else {
      mins <- ntmp <- damCCutInt[ijsd]*damYInt[ijsd]*distances[ijsd,]
    }
    nclct <- length(which(mins<23^2))
    if(min(mins)<= clctDist**2){
      if(length(xi)>1 & length(ijsd)>1){
        ijs <- ijsd[unique(apply(ntmp,2,which.min)[which(mins==min(mins))])]
        ijs <- ijs[which.min(dam_yeari-dam_yeard[ntmpi[ijs]])]
      } else if(ncol(distances)==1) {
        ijs <- ijsd[which(mins==min(mins))]
        ijs <- ijs[which.min(dam_yeari-dam_yeard[ijs])]
      } else {
        ijs <- ijsd 
      }
      #print(paste0(r_no," CC neighbors found, ij=",ij))
      if(length(ijs)>0 & KUVA){
        print(paste0("CC neighbors found, ij=",ij))
        par(mfrow=c(1,1))
        plot(c(xxd[ntmpi[ijsd]],xi)-xi[1],c(yyd[ntmpi[ijsd]],yi)-yi[1],pch=19,
             main=paste0("CC, row=",ij,", damSegid=",damSegIDij))
        points(xi-xi[1],yi-yi[1],col="red",pch=19)
        points(xxd[ntmpi[ijs]]-xi[1],yyd[ntmpi[ijs]]-yi[1],col="green",pch=4,cex=1.6)
      }
    }
    
    # SBB in neighborhood
    if(ncol(distances)>1 & length(ijsd)>1){
      ntmp <- damBBInt[ijsd]*damYInt[ijsd]*distances[ijsd,]
      mins <- apply(ntmp,2,min)
    } else if(ncol(distances)==1) {
      ntmp <- damBBInt[ijsd]*damYInt[ijsd]*distances[ijsd,1]
      mins <- ntmp
    } else {
      mins <- ntmp <- damBBInt[ijsd]*damYInt[ijsd]*distances[ijsd,]
    }
    nsbb <- length(which(mins<23^2))
    #  ntmp <- min(damCCutInt*damYInt*distances)
    if(1%in%unique(damBBInt) & min(mins)<= clctDist**2){
      if(length(xi)>1 & length(ijsd)>1){
        ijsSBB <- ijsd[unique(apply(ntmp,2,which.min)[which(mins==min(mins))])]
        ijsSBB <- ijsSBB[which.min(dam_yeari-dam_yeard[ijsSBB])]
      } else if(ncol(distances)==1) {
        ijsSBB <- ijsd[which(mins==min(mins))]
        ijsSBB <- ijsSBB[which.min(dam_yeari-dam_yeard[ijsSBB])]
      } else {
        ijsSBB <- ijsd 
      }
      
      #print(paste(1%in%unique(damBBInt),"/",sqrt(min(mins))))
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
    if(ncol(distances)>1 & length(ijsd)>1){
      ntmp <- damWInt[ijsd]*damYInt[ijsd]*distances[ijsd,]
      mins <- apply(ntmp,2,min)
    } else if(ncol(distances)==1) {
      ntmp <- damWInt[ijsd]*damYInt[ijsd]*distances[ijsd,1]
      mins <- ntmp
    } else {
      mins <- ntmp <- damWInt[ijsd]*damYInt[ijsd]*distances[ijsd,]
    }
    nclct <- length(which(mins<23^2))
    #  ntmp <- min(damCCutInt*damYInt*distances)
    if(1%in%unique(damWInt) & min(mins)<= clctDist**2){
      if(length(xi)>1 & length(ijsd)>1){
        ijsWind <- ijsd[unique(apply(ntmp,2,which.min)[which(mins==min(mins))])]
        ijsWind <- ijsWind[which.min(dam_yeari-dam_yeard[ijsWind])]
      } else if(ncol(distances)==1) {
        ijsWind <- ijsd[which(mins==min(mins))]
        ijsWind <- ijsWind[which.min(dam_yeari-dam_yeard[ijsWind])]
      } else {
        ijsWind <- ijsd 
      }
      
      #print(paste0(r_no, " Wind neighbors found, ij=",ij))
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
    
    xx <- xxd[ntmpi[ijs]]
    yy <- yyd[ntmpi[ijs]]
    dam_year <- dam_yeard[ntmpi[ijs]]
    
    neighbors <-function(ik,clctDist){
      dclct <- dclct_south <- c(1e12,NA)
      dx <- sqrt((xx-xi[ik])^2+(yy-yi[ik])^2)
      if(dx[order(dx)[1]]<dclct[1] & dx[order(dx)[1]]<clctDist & dx[order(dx)[1]]>0){
        dclct <- c(dx[order(dx)[1]],dam_year[order(dx)[1]]-dam_yeari)  
      }
      # close south neighbours with clearcut: distance and time from clearcut
      dy <- yy[which(xx<(xi[ik]+32) & xx>(xi[ik]-32))]
      dy <- yi[ik]-dy[dy<yi[ik]]
      if(length(dy)>0){
        if(dy[order(dy)[1]]<dclct_south[1] & dy[order(dy)[1]]<clctDist){#dy[order(dy)[1]] 
          dclct_south <- c(dx[order(dy)[1]],dam_year[order(dy)[1]]-dam_yeari)
        }
      }
      return(c(dclct,dclct_south))
    }
    tmp <- apply(data.table(1:length(iks)),1:2,neighbors,clctDist=clctDist)
    dclct <- tmp[1:2,which.min(tmp[1,,1]),1]
    dclct_south <- tmp[3:4,which.min(tmp[3,,1]),1]
    #if(dclct_south[1]<=clctDist) print(paste0(r_no, " south cc neighbors found, ij=",ij))
    
  }
  if(length(ijsSBB)>0){
    ik<-1
    # close neighbours with SBB: distance and time from clearcut
    xx <- xxd[ntmpi[ijsSBB]]
    yy <- yyd[ntmpi[ijsSBB]]
    dam_year <- dam_yeard[ntmpi[ijsSBB]]
    dx <- sqrt((xx-xi[ik])^2+(yy-yi[ik])^2)
    
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
  }
  if(length(ijsWind)>0){
    ik<-1
    xx <- xxd[ntmpi[ijsWind]]
    yy <- yyd[ntmpi[ijsWind]]
    dam_year <- dam_yeard[ntmpi[ijsWind]]
    
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
    
  }    
  out <- c(dclct, dclct_SBB, dclct_Wind, dclct_south)
  #if(min(out[c(1,3,5,7)])<1e12){
  #  print(ij)
  #  print(out)
  #} 
  return(out)
}

