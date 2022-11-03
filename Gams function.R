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
  Data.miss <- prodNA(xdata, noNA = pro_mis)
  
  ##1
  ms_completeData <- as.data.frame(Data.miss)
  for(i in 1:ncol(ms_completeData)) {
    ms_completeData[ , i][is.na(ms_completeData[ , i])] <- mean(ms_completeData[ , i], na.rm = TRUE)
  }
  ms_completeData <- ms_completeData
  as.matrix(ms_completeData)
  ##2
  imputed_Data <- mice(Data.miss, m=1, method = 'pmm', seed = 500)
  mi_completeData <- as.data.frame(complete(imputed_Data))
  as.matrix(mi_completeData)
  ##3
  rf_completeData <- rfImpute(Y ~ ., Data.miss)[,2:4]
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





Test_mse <- function(Data,EY, method, prop_miss, loop) {
  #select data
  Y <- Data[,1]
  n <- ncol(Data)
  xdata <- Data[,2:n]
  #output
  df_mean <- data.frame(matrix(ncol = n+1, nrow = loop))
  df_mice <- data.frame(matrix(ncol = n+1, nrow = loop))
  df_rf <- data.frame(matrix(ncol = n+1, nrow = loop))
  
  for (i in 1:loop){
    #drop data
    Data.miss <- prodNA(xdata, noNA = prop_miss)
    completedcase <- nrow(na.omit(Data.miss))
    ## mean substitute
    mean_imputedData <- as.data.frame(Data.miss)
    for(j in 1:ncol(mean_imputedData)) {
      mean_imputedData[ , j][is.na(mean_imputedData[ , j])] <- mean(mean_imputedData[ , j], na.rm = TRUE)
    }
    as.matrix(mean_imputedData)
    ## Mice
    #You may skip removal of collinear variables at initialization by the remove_collinear flag.
    imputed_Data <- mice(Data.miss, method = 'pmm', diagnostics = FALSE , remove_collinear = FALSE, seed = 500)
    mice_imputedData <- as.data.frame(complete(imputed_Data))
    as.matrix(mice_imputedData)
    ## Random Forest
    rf_imputedData <- rfImpute(Y ~ ., Data.miss)[,2:n]
    as.matrix(rf_imputedData)
    #
    Data <- as.data.frame(Data)
    O <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=Data,method="REML")
    TO <- predict.gam(O,type="terms",as.data.frame(xdata))
    pO <- predict.gam(O,as.data.frame(xdata))
    #
    Me <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=mean_imputedData,method="REML")
    TMe <- predict.gam(Me,type="terms",as.data.frame(xdata))
    DTMe <- colSums((TMe-TO)^2)
    pMe <- predict.gam(Me,as.data.frame(xdata))
    #
    M <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=mice_imputedData,method="REML")
    TM <- predict.gam(M,type="terms",as.data.frame(xdata))
    DTM <- colSums((TM-TO)^2)
    pM <- predict.gam(M,as.data.frame(xdata))
    #
    rf <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=rf_imputedData,method="REML")
    Trf <- predict.gam(rf,type="terms",as.data.frame(xdata))
    DTrf <- colSums((Trf-TO)^2)
    prf <- predict.gam(rf,as.data.frame(xdata))
    #
    mse_Y <- c(sum((pO-EY)^2),sum((pMe-EY)^2),sum((pM-EY)^2),sum((prf-EY)^2))
    df_mean[i,] <- c(mse_Y[2],DTMe,completedcase)
    df_mice[i,] <- c(mse_Y[3],DTM,completedcase)
    df_rf[i,] <- c(mse_Y[4],DTrf,completedcase)
  } 
  if (method == "mean"){
    return(df_mean)
  }
  else if(method == "mice"){
    return(df_mice)
  }
  else if(method == "rf"){
    return(df_rf)
  }
  else{
    return()
  }
}

Test_mse1 <- function(Data,EY, method, prop_miss, loop) {
  #select data
  Y <- Data[,1]
  n <- ncol(Data)
  xdata <- Data[,2:n]
  #output
  df_mean <- data.frame(matrix(ncol = n+1, nrow = loop))
  df_mice <- data.frame(matrix(ncol = n+1, nrow = loop))
  df_rf <- data.frame(matrix(ncol = n+1, nrow = loop))
  
  for (i in 1:loop){
    #drop data
    Data.miss <- prodNA(xdata, noNA = prop_miss)
    completedcase <- nrow(na.omit(Data.miss))
    #
    Data <- as.data.frame(Data)
    O <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=Data,method="REML")
    TO <- predict.gam(O,type="terms",as.data.frame(xdata))
    pO <- predict.gam(O,as.data.frame(xdata))
    
    if (method == "mean"){
      ## mean substitute
      mean_imputedData <- as.data.frame(Data.miss)
      for(j in 1:ncol(mean_imputedData)) {
        mean_imputedData[ , j][is.na(mean_imputedData[ , j])] <- mean(mean_imputedData[ , j], na.rm = TRUE)
      }
      as.matrix(mean_imputedData)
      Me <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=mean_imputedData,method="REML")
      TMe <- predict.gam(Me,type="terms",as.data.frame(xdata))
      DTMe <- colSums((TMe-TO)^2)
      pMe <- predict.gam(Me,as.data.frame(xdata))
      mse_Y <- c(sum((pMe-EY)^2))
      df_mean[i,] <- c(mse_Y,DTMe,completedcase)
      return(df_mean)
    }
    else if(method == "mice"){
      ## Mice
      #You may skip removal of collinear variables at initialization by the remove_collinear flag.
      imputed_Data <- mice(Data.miss, method = 'pmm', diagnostics = FALSE , remove_collinear = FALSE, seed = 500)
      mice_imputedData <- as.data.frame(complete(imputed_Data))
      as.matrix(mice_imputedData)
      #
      M <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=mice_imputedData,method="REML")
      TM <- predict.gam(M,type="terms",as.data.frame(xdata))
      DTM <- colSums((TM-TO)^2)
      pM <- predict.gam(M,as.data.frame(xdata))
      mse_Y <- c(sum((pM-EY)^2))
      df_mice[i,] <- c(mse_Y,DTM,completedcase)
      return(df_mice)
    }
    else if(method == "rf"){
      ## Random Forest
      rf_imputedData <- rfImpute(Y ~ ., Data.miss)[,2:n]
      as.matrix(rf_imputedData)
      #
      rf <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=rf_imputedData,method="REML")
      Trf <- predict.gam(rf,type="terms",as.data.frame(xdata))
      DTrf <- colSums((Trf-TO)^2)
      prf <- predict.gam(rf,as.data.frame(xdata))
      mse_Y <- c(sum((prf-EY)^2))
      df_rf[i,] <- c(mse_Y,DTrf,completedcase)
      return(df_rf)
    }
    else{
      return()
    }
  }
}

