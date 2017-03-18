library(ggplot2)

setwd('..')
source("./2_Scripts/spaInterpol.R")

#reading distance matrix,coordinates and event data from corrosponding csv files
distMat=as.matrix(read.csv("./1_Data/DistanceMatrix.csv",header=FALSE))
rainGauPos=as.matrix(read.csv("./1_Data/loc.csv",header=TRUE))
eventData=read.csv("./1_Data/eventData.csv",header=TRUE)

#event(s)
selEvents=c(8)

#gauge number
gauNum=6

#desired timescale(s) in minutes
timeScaleList=c("180 min","360 min")
timeScaleList=c("5 min","15 min","30 min","60 min")
timeScaleList=c("5 min")

#raingauges
gauMean=TRUE #if the mean value of pair raingauge should be used or not
rgNo=seq(1,15,2) #raingauge numbers  (keep as it is if true is selected above) 

#saving plots
#png(paste("C:/QUICS/4. Conferences&Workshops/ECAM2015/3.Plots/",gauNum, ".png"),
#width=480*3,height=480)
#par(mfrow=c(1,2))

##############Kriging###############

for (k in selEvents) {
  
  #reading from csv files and extracting data for the event date and time
  fileName=paste(eventData[k,2])
  eventStart=as.POSIXct(eventData[k,3],tz="",format = "%d/%m/%Y %H:%M")
  eventEnd=as.POSIXct(eventData[k,4],tz="",format = "%d/%m/%Y %H:%M")
  eventDur=as.numeric(difftime(eventEnd, eventStart, units="hours"))
  
  
  for (j in timeScaleList) {
    
    #Extracting time series for a particular time scale
    tsInt=csvToInten("./1_Data/",fileName,eventStart,eventEnd,timeScale=j)
    
    
    if(gauMean==TRUE) {
      #averaging of two rainguages
      avrIntAll=NULL
      for (i in rgNo) {
        avrInt=apply(tsInt[,i:(i+1)],1,mean)
        avrIntAll=cbind(avrIntAll,avrInt)
      } 
    } else {
      #using single rain gauges
      avrIntAll=tsInt[,rgNo]
    }
    
    colnames(avrIntAll)=1:ncol(avrIntAll)
    #allTsInt=list(allTsInt,avrIntAll)
    
    #setting a intensity threshold
    intThres=0
    rowMean=rowMeans(avrIntAll)
    avrIntThresHigh=avrIntAll[which(rowMean>intThres),]
    avrIntThresLow=avrIntAll[which(rowMean<intThres),]
    
    #using kriging function
    krigInt=kriging_Single(rainGauPos,avrIntThresHigh,gauNum)
    measuredInt=krigInt[[3]]
    predictedInt=krigInt[[2]]
    
    #plotting
    print(qplot(measuredInt,predictedInt,
         xlab="Measured (mm/hr)",main=paste("Station - ",gauNum,":",
                                            "Time scale - ",j,"(Kriging)"),
         ylab="predicted (mm/hr)",xlim=c(0,20),ylim=c(0,20)))
         #cex.axis=1.5,cex.lab=1.5,cex.main=2)
    #textxy(avrIntThresHigh[,gauNu][-52],predictAll[-52],labs=seq(1,nrow(predictAll),1))
    #abline(a=0,b=1)
    
  }
}

##############idw###############


for (k in selEvents) {
  
  #reading from csv files and extracting data for the event date and time
  fileName=paste(eventData[k,2])
  eventStart=as.POSIXct(eventData[k,3],tz="",format = "%d/%m/%Y %H:%M")
  eventEnd=as.POSIXct(eventData[k,4],tz="",format = "%d/%m/%Y %H:%M")
  eventDur=as.numeric(difftime(eventEnd, eventStart, units="hours"))
  
  allTsInt=NULL
  for (j in timeScaleList) {
    
    #Extracting time series for a particular time scale
    tsInt=csvToInten("./1_Data/",fileName,eventStart,eventEnd,timeScale=j)
    
    
    if(gauMean==TRUE) {
      #averaging of two rainguages
      avrIntAll=NULL
      for (i in rgNo) {
        avrInt=apply(tsInt[,i:(i+1)],1,mean)
        avrIntAll=cbind(avrIntAll,avrInt)
      } 
    } else {
      #using single rain gauges
      avrIntAll=tsInt[,rgNo]
    }
    
    colnames(avrIntAll)=1:ncol(avrIntAll)
    allTsInt=list(allTsInt,avrIntAll)
    
    #setting a intensity threshold
    intThres=0
    rowMean=rowMeans(avrIntAll)
    avrIntThresHigh=avrIntAll[which(rowMean>intThres),]
    avrIntThresLow=avrIntAll[which(rowMean<intThres),]
    
    #using kriging function
    idwInt=idwInterpolation(distMat,avrIntThresHigh,gauNum)
    measuredInt=idwInt[[3]]
    predictedInt=idwInt[[2]]
    
    #plotting
    qplot(measuredInt,predictedInt,
         xlab="Measured (mm/hr)",main=paste("Station - ",gauNum,":",
                                            "Time scale - ",j,"(IDW)"),
         ylab="predicted (mm/hr)",xlim=c(0,20),ylim=c(0,20),
         cex.axis=1.5,cex.lab=1.5,cex.main=2)
    #textxy(avrIntThresHigh[,gauNu][-52],predictAll[-52],labs=seq(1,nrow(predictAll),1))
    abline(a=0,b=1)
    
  }
}



#plotting average prediction
#     meanIntFin=rowMeans(intFin)
#     plot(avrIntThresHigh[,gauNu],meanIntFin, xlab="Measured (mm/hr)",
#          main=paste("Station - ",gauNu,":","Time scale - ",j,"(Mean)"),
#          ylab="predicted (mm/hr)",xlim=c(0,20),ylim=c(0,20),
#          cex.axis=1.5,cex.lab=1.5,cex.main=2)
#     abline(a=0,b=1)

#generating grid
#     xcord=seq(415300,415750,5)
#     ycord=seq(432700,432900,5)
#     gridList=lapply(xcord, function(x) matrix(c(rep(x[1],length(ycord)),ycord),ncol=2))
#     gridMatrix=do.call(rbind,gridList)
#





