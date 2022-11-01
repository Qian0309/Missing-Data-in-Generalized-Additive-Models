## gap of range()
dif <- function(x) {
  dif <- range(x)[2]-range(x)[1]
  return(dif)
}

## coefficients of functions
cof <- function(x) {
  C <- c(rep(0,ncol(x)-1))
  for (i in 2:ncol(x)){
    C[i-1] <- dif((x[,1])/dif(x[,i]))
  }
  return(C)
}

# test multivariate normal distribution
assessment <- function(Data, pop_dis, pro_mis) {
  #select data
  Y <- Data[,1]
  n <- ncol(Data)
  xdata <- Data[,2:n]
  #drop data
  Data.mis <- prodNA(xdata, noNA = pro_mis)
  
  ##1
  ms_completeData <- as.data.frame(Data.mis)
  for(i in 1:ncol(ms_completeData)) {
    ms_completeData[ , i][is.na(ms_completeData[ , i])] <- mean(ms_completeData[ , i], na.rm = TRUE)
  }
  ms_completeData <- ms_completeData
  as.matrix(ms_completeData)
  ##2
  imputed_Data <- mice(Data.mis, m=1, method = 'pmm', seed = 500)
  mi_completeData <- as.data.frame(complete(imputed_Data))
  as.matrix(mi_completeData)
  ##3
  rf_completeData <- rfImpute(Y ~ ., Data.mis)[,2:4]
  as.matrix(rf_completeData)
  
  #H0 (null): The variables not follow a multivariate normal distribution.
  #Ha (alternative): The variables follow a multivariate normal distribution.
  result <- mvn(xdata,mvnTest = "hz")
  result1 <- mvn(ms_completeData,mvnTest = "hz")
  result2 <- mvn(mi_completeData,mvnTest = "hz")
  result3 <- mvn(rf_completeData,mvnTest = "hz")
  
  #hoteling t square distribution test
  #H0 mu same
  #Ha mu not same
  t1 <- hotelling.stat(xdata,as.matrix(ms_completeData))$statistic
  t2 <- hotelling.stat(xdata,as.matrix(mi_completeData))$statistic
  t3 <- hotelling.stat(xdata,as.matrix(rf_completeData))$statistic
  
  Data <- as.data.frame(Data)
  O <- gam(Y~s(X1)+s(X2)+s(X3),data=Data,method="REML")
  pO <- predict.gam(O,as.data.frame(xdata))
  Me <- gam(Y~s(X1)+s(X2)+s(X3),data=ms_completeData,method="REML")
  pMe <- predict.gam(Me,as.data.frame(xdata))
  M <- gam(Y~s(X1)+s(X2)+s(X3),data=mi_completeData,method="REML")
  pM <- predict.gam(M,as.data.frame(xdata))
  rf <- gam(Y~s(X1)+s(X2)+s(X3),data=rf_completeData,method="REML")
  prf <- predict.gam(rf,as.data.frame(xdata))
  rss <- c(sum((pO-Y)^2),sum((pMe-Y)^2),sum((pM-Y)^2),sum((prf-Y)^2))
  
  res <- c(result$multivariateNormality$MVN,result1$multivariateNormality$MVN,result2$multivariateNormality$MVN,
           result3$multivariateNormality$MVN,t1,t2,t3,rss)
  return(res)
}
