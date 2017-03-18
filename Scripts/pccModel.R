source("./2_Scripts/gen_util.R")

###pearson correlation matrix### 

#TSIntensity: Time series of rainfall intensity for a particular time scale
PCCMatrix=function(TSIntensity){
  n=seq_len(ncol(TSIntensity))
  ff=function (a,b) cor.test(TSIntensity[,a], TSIntensity[,b])[[4]]
  corrMat=outer(n,n,Vectorize(ff))
  return(corrMat)
}
