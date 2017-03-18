setwd('..')
source("./2_Scripts/pccModel.R")


#reading distance matrix from csv file
distMat=as.matrix(read.csv("./1_Data/DistanceMatrix.csv", header=FALSE))
eventData=read.csv("./1_Data/eventData.csv",header=TRUE)

selEvents=c(1)
allInt=NULL


##############Analysis:1-dependancy on timescale###############

for (k in selEvents) {
  
  #reading from csv files and extracting data for the event date and time
  fileName=paste(eventData[k,2])
  eventStart=as.POSIXct(eventData[k,3],tz="",format = "%d/%m/%Y %H:%M")
  eventEnd=as.POSIXct(eventData[k,4],tz="",format = "%d/%m/%Y %H:%M")
  eventDur=as.numeric(difftime(eventEnd, eventStart, units="hours"))
  
  #extracting event data
  tsEvent=eventTS("./1_Data/",fileName,eventStart,eventEnd)
  
  #calculating average intensity
  eventInt=round(mean((coredata(tsEvent[index(tsEvent)==eventEnd])
                       -coredata(tsEvent[index(tsEvent)==eventStart]))/eventDur),2)
  allInt=rbind(allInt,eventInt)
  
  #plot (rainfall timeseries)
  png(paste("./3_Results/",k, ".png"),
      width=480*1.5,height=480*2)
  par(mfrow=c(2,1))
  plot(NULL,NULL,ylim=c(0,(max(tsEvent)-min(tsEvent))),xlim=c(eventStart,eventEnd),
       xlab=paste("Time (hrs) :", format(eventStart,"%d/%m/%Y"),"to",
                  format(eventEnd,"%d/%m/%Y")),
       ylab="Cum. Rainfall(mm)",xaxt='n',main=paste("Event No -",k))
  legend("topleft",legend=c(paste("Duration = ",eventDur,"hrs"),
                            paste("Avr.Intensty = ",eventInt,"mm/hr")),
         bty="n",cex=1.5)
  axis.POSIXct(1, at = seq(eventStart, eventEnd, by = "hour"),format="%H:%M")
  #grid(nx=25,ny=NULL)
  abline(v=seq(as.POSIXct(eventStart),as.POSIXct(eventEnd),"1 hour"),col="gray",
         lty=2)
  
  for (i in 1:16){
    lines((tsEvent[,i]-min(tsEvent[,i])),col=i)
  }
  
  #desired timescale(s) in minutes
  timeScaleList=c("180 min","360 min")
  timeScaleList=c("2 min","5 min","15 min","30 min","60 min") 
  
  #Correaltion Coefficient Vs Seperation distance
  plot(NULL,NULL,xlim=c(0,500),ylim=c(0.5,1),xlab="Seperation Distance (m)", 
       ylab="Correlation Coefficient")
  grid(nx=NULL,ny=NULL)
  clip(50,1000,-100,100) #to set abline boundaries 
  PCCAll=NULL
  for (j in timeScaleList) {
    
    #Extracting time series for a particular time scale
    tsScaled=extractTS(tsEvent,j)
    tsStcum=cumToStAcc(tsScaled)
    tsInt=stAccToInt(tsStcum,j)
    
    #averaging of two rainguages
    avrIntAll=NULL
    rgNo=seq(1,15,2)
    
    for (i in rgNo) {
      avrInt=apply(tsInt[,i:(i+1)],1,mean)
      avrIntAll=cbind(avrIntAll,avrInt)
    }
    
    #     #omitting zeros
    #     intThres=0
    #     rowMean=rowMeans(avrIntAll)
    #     avrIntAll=avrIntAll[which(rowMean>intThres),]
    
    #generating pearson correlation matrix
    PCCMat=PCCMatrix(avrIntAll)
    PCCAll=abind(PCCAll,PCCMat,along=3)
    
    #Removing diagonal elements to avoid the plot starting from (0,1)
    diag(PCCMat)=NA
    diag(distMat)=NA
    
    #plotting
    corrVec=c(PCCMat)
    distVec=c(distMat)
    linMod=lm(corrVec~distVec)
    
    colList=c("blue","red","green","orange","violet","black")
    #colList=c("blue","red","green","orange")
    points(distVec,corrVec,col=colList[which(timeScaleList==j)],
           pch=which(timeScaleList==j))
    abline(a=linMod[[1]][1],b=linMod[[1]][2],col=colList[which(timeScaleList==j)])
  }
  
  legend("bottom",legend=timeScaleList,pch=c(1:6),bty="n",ncol=6,lty=1,col=colList)
  dev.off()
}


##############Analysis:2-dependancy on intensity###############

intThres=12 #threshold intensity to be used

for (k in selEvents) {
  
  #reading from csv files and extracting data for the event date and time
  fileName=paste(eventData[k,2])
  eventStart=as.POSIXct(eventData[k,3],tz="",format = "%d/%m/%Y %H:%M")
  eventEnd=as.POSIXct(eventData[k,4],tz="",format = "%d/%m/%Y %H:%M")
  eventDur=as.numeric(difftime(eventEnd, eventStart, units="hours"))
  tsEvent=eventTS("./1_Data/",fileName,eventStart,eventEnd)
  eventInt=round(mean((coredata(tsEvent[index(tsEvent)==eventEnd])-
                         coredata(tsEvent[index(tsEvent)==eventStart]))/eventDur),2)
  allInt=rbind(allInt,eventInt)
  
  png(paste("./3_Results/",k, ".png"),
      width=480*1.5,height=480)
  par(mfrow=c(1,1))
  
  #desired timescale(s) in minutes
  timeScaleList=c("15 min")
  #timeScaleList=c("2 min","5 min","15 min","30 min","60 min","180 min") 
  
  #Correaltion Coefficient Vs Seperation distance
  plot(x=NULL,y=NULL,xlim=c(0,400),ylim=c(0.75,1),
       xlab="Seperation Distance (m)", ylab="Correlation Coefficient",cex=2, 
       main="Time scale = 15 min",font.lab=2)
  axis(side = 1, font = 2, cex=2)
  axis(side = 2, font = 2, cex=2)
  grid(nx=NULL,ny=NULL)
  clip(20,1000,0,1) #to set abline boundaries
  
  for (j in timeScaleList) {
    
    #Extracting time series for a particular time scale
    tsScaled=extractTS(tsEvent,j)
    tsStcum=cumToStAcc(tsScaled)
    tsInt=stAccToInt(tsStcum,j)
    
    #averaging of two rainguages
    avrIntAll=NULL
    rgNo=seq(1,15,2)
    
    for (i in rgNo) {
      avrInt=apply(tsInt[,i:(i+1)],1,mean)
      avrIntAll=cbind(avrIntAll,avrInt)
    }
    
    #setting a intensity threshold
    rowMean=rowMeans(avrIntAll)
    avrIntThresHigh=avrIntAll[which(rowMean>intThres),]
    avrIntThresLow=avrIntAll[which(rowMean<intThres),]
    
    
    #generating pearson correlation matrix
    PCCMatHigh=PCCMatrix(avrIntThresHigh)
    PCCMatLow=PCCMatrix(avrIntThresLow)
    
    #Removing diagonal elements
    diag(PCCMatHigh)=NA
    diag(PCCMatLow)=NA
    diag(distMat)=NA
    
    #plotting
    corrVecHigh=c(PCCMatHigh)
    distVecHigh=c(distMat)
    
    corrVecLow=c(PCCMatLow)
    distVecLow=c(distMat)
    
    #linear Model
    linModHigh=lm(corrVecHigh~distVecHigh)
    linModLow=lm(corrVecLow~distVecLow)
    
    #Polynomial model
    #     offset=rep(1,length(corrVec))
    #     polMod=lm(corrVec~(-1+distVec+I(distVec^2)+offset(offset)))
    #     polFun=function(x)(polMod[[1]][[2]]*(x^2)+polMod[[1]][[1]]*x+1)
    
    #exponential model
    #     expMod=lm(corrVec~log(distVec+0.01))
    #     expFun=function(x)(expMod[[1]][[2]]*log(x)+expMod[[1]][[1]])
    
    colListH=c("blue","red","green","orange","violet","black")
    colListL=c("black","red","green","orange","violet","black")
    
    #colList=c("blue","red","green","orange")
    points(distVecHigh,corrVecHigh,col=colListH[which(timeScaleList==j)],pch=1)
    points(distVecLow,corrVecLow,col=colListL[which(timeScaleList==j)],pch=2)
    
    abline(a=linModHigh[[1]][1],b=linModHigh[[1]][2],
           col=colListH[which(timeScaleList==j)])
    abline(a=linModLow[[1]][1],b=linModLow[[1]][2],
           col=colListL[which(timeScaleList==j)])
    #abline(a=linMod[[1]][1],b=linMod[[1]][2],col=colList[which(timeScaleList==j)],
    #      untf=TRUE)
    #plot(polMod)
    #curve(polFun,add=TRUE,col=colList[which(timeScaleList==j)])
    #curve(expFun,add=TRUE,col=colList[which(timeScaleList==j)])
  }
  #   legend("bottom",legend=timeScaleList,pch=c(1),bty="n",ncol=5,lty=1,
  #          col=colListH)
  legend("bottom",legend=c("> 12mm/hr","< 12mm/hr"),pch=c(1,2),bty="n",ncol=2,lty=1,
         cex=1.2,col=c("blue","black"))
  dev.off()
}



#bin
# to get the required average
# a=index(precApr25)
# b=cut(a,breaks="15 min")
# c=aggregate(precApr25,b,mean)

#precApr=read.csv(paste0(directoryData,"apr_min.csv"),header = TRUE)
#precApr[,1]<-as.POSIXct(precApr[,1], format="%d/%m/%Y %H:%M")
#precAprxts=as.xts(precApr)
#precApr25=precAprxts[.indexDate(25)]
