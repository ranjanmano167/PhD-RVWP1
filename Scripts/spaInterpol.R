library(geoR)
library(hydroGOF)

source("./2_Scripts/gen_util.R")


##########Kriging##############
kriging_Single=function(coordinates,intensity,NewCord,method="ordinary",
                        uvect="default",covModel="exponential") {
  
  coordinatesKrig=coordinates
  intensityKrig=intensity
    
  data=cbind(coordinatesKrig,intensityKrig)
  geodata=as.geodata(data,header=TRUE)
  variogram=variog(geodata,coords=geodata$coords,data=geodata$data,
                   uvec=uvect,messages=FALSE)
  vplot=plot(variogram,main=i)
  varioPara=variofit(variogram,cov.model=covModel,messages=FALSE)
  lines(varioPara)
  prediction=ksline(geodata,cov.model=varioPara[[3]],
                    cov.pars=c(varioPara[[2]][1],varioPara[[2]][2]),
                    nugget=varioPara[[1]],locations=coordinates[gauNo,],
                    messages=FALSE)
  
  return(prediction)
}


kriging_Com=function(coordinates,intensity,gauNo=6,method="ordinary",
                 uvect="default",covModel="exponential") {
  
  coordinatesKrig=coordinates[-gauNo,]
  intensityKrig=intensity[,-gauNo]
  measuredInt=matrix(coredata(intensity[,gauNo]),ncol=1)
  dateAndTime=data.frame(index(intensity))
  
  predictedInt=NULL
  
  for (i in 1:nrow(intensityKrig)) {
    data=cbind(coordinatesKrig,intensityKrig[i,])
    geodata=as.geodata(data,header=TRUE)
    variogram=variog(geodata,coords=geodata$coords,data=geodata$data,
                     uvec=uvect,messages=FALSE)
    vplot=plot(variogram,main=i)
    varioPara=variofit(variogram,cov.model=covModel,messages=FALSE)
    lines(varioPara)
    prediction=ksline(geodata,cov.model=varioPara[[3]],
                      cov.pars=c(varioPara[[2]][1],varioPara[[2]][2]),
                      nugget=varioPara[[1]],locations=coordinates[gauNo,],
                      messages=FALSE)
    predictedInt=rbind(predictedInt,prediction[[1]])
    
  }
  
  #plot(measuredInt,predictedInt)
  #abline(a=0,b=1)
  gof=gof(predictedInt,measuredInt)
  ResultKriging=list(dateAndTime,predictedInt,measuredInt,gof)
  return(ResultKriging)
}



##########inverse distance weighting##############  


idwInterpolation=function(distMatrix,intensity,gauNo=6) {
  
  distMat=matrix(distMatrix[,gauNo][-gauNo],ncol=1)
  intensityIDW=coredata(intensity[,-gauNo])
  measuredInt=matrix(coredata(intensity[,gauNo]),ncol=1)
  dateAndTime=data.frame(index(intensity))
  
  predictedInt=NULL
  
  for (i in 1:nrow(intensityIDW)) {
    invDistMat=1/distMat
    prediction=sum(intensityIDW[i,]*invDistMat)/sum(invDistMat)
    predictedInt=rbind(predictedInt,prediction[])    
  }
  
  #plot(measuredInt,predictedInt)
  #abline(a=0,b=1)
  gof=gof(predictedInt,measuredInt)
  ResultKriging=list(dateAndTime,predictedInt,measuredInt,gof)
  return(ResultKriging)
}

