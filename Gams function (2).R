## gap of range()
diff <- function(x) {
  diff <- range(x)[2]-range(x)[1]
  return(diff)
}

## rescale ratio of f 
rescale <- function(fx,Term_fx) {
  re_f <- data.frame(matrix(ncol = ncol(fx), nrow = nrow(fx)))
  for (i in 1:ncol(fx)){
    re_f[,i] <- fx[,i]/diff(fx[,i])*diff(Term_fx[,i])
  }
  return(re_f)
}

# test multivariate normal distribution

#f4 <- f4-sum(f4)/300
#c2 <- c(f4/diff(f4)*diff(TMean[,4]))
#DTMean <- sqrt((sum(TMean[,4]-c2)^2)/300)
#matplot(xdata[,4], cbind(c2,TMean[,4]),col=c("red","green"),lty=c(1,1))

Test_rmse <- function(Data, EY, method, prop_miss, loop) {
  
  # Define Y
  Y <- Data[,1]
  # Number of columns
  # 1 Y, n xi's and n fx_i's (2n+1 columns in total)
  N <- ncol(Data)
  # Number of observations of each variable
  r <- nrow(Data)
  # locate the column of fx_1
  col_f1 <- (N-1)/2+2
  # Define all columns that containing fx_i's data
  fx_Data <- Data[,col_f1:N]
  
  # Define a dataset containing Y and x_i's, excluding all fx_i's columns
  Data <- Data[,1:col_f1-1]
  # Number of coulmns in dataframe "Data"
  n <- ncol(Data)
  # Define a dataset only containing x_i's, excluding all Y and fx_i's columns
  xdata <- Data[,2:ncol(Data)]
  
  # Define matrix for storing Y_rmse, x_rmse and number of complete case
  df_mean <- data.frame(matrix(ncol = n+1, nrow = loop))
  df_mice <- data.frame(matrix(ncol = n+1, nrow = loop))
  df_rf <- data.frame(matrix(ncol = n+1, nrow = loop))
  
  # Apply GAMs to the non-dropping full dataframe
  Data <- as.data.frame(Data)
  O <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=Data,method="REML")
  TO <- predict.gam(O,type="terms",as.data.frame(xdata))
  pO <- predict.gam(O,as.data.frame(xdata))
  
  # Compute root mean squared error for y variable and each x variables after imputing the missing data
  # If the "Mean" imputation method was chosen
  if (method == "mean"){
    for (i in 1:loop){
      
      # Drop data
      Data.miss <- prodNA(xdata, noNA = prop_miss)
      completedcase <- nrow(na.omit(Data.miss))
    
      ## Apply "Mean" imputation method on dataframe with missing data
      mean_imputedData <- as.data.frame(Data.miss)
      for(j in 1:ncol(mean_imputedData)) {
        mean_imputedData[ , j][is.na(mean_imputedData[ , j])] <- mean(mean_imputedData[ , j], na.rm = TRUE)
      }
      as.matrix(mean_imputedData)
      
      # Apply GAMs to the full dataframe after mean imputation
      Mean_mod <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=mean_imputedData,method="REML")
     
      # GAMs prediction for each xi
      x_pred <- predict.gam(Mean_mod,type="terms",as.data.frame(xdata))  
      
      # Rescale fx_i's
      re_fx <- rescale(fx_Data, x_pred)
      for (m in 1:ncol(re_fx)){
        re_fx[, m] <- c(re_fx[, m] - rep(colMeans(re_fx)[m], r))
      }
      
      # Compute root mean squared error (rmse) for each xi's
      x_rmse <- sqrt(colSums((x_pred-re_fx)^2)/r)
      
      # GAMs prediction for y
      Y_pred <- predict.gam(Mean_mod,as.data.frame(xdata))
      # Compute root mean squared error (rmse) for y
      Y_rmse <- c(sqrt(sum((Y_pred-EY)^2)/r))
      
      # Add the above computed x_rmse, y_rmse and the number of complete cases after dropping data
      # into the matrix previously created (eg. loop 1 results go into the first row of the matrix)
      df_mean[i,] <- c(Y_rmse,x_rmse,completedcase)

    }
    return(df_mean)
  }
  
  # Compute root mean squared error for y variable and each x variables after imputing the missing data
  # If the "Mice" imputation method was chosen
  else if(method == "mice"){
    for (i in 1:loop){
      
      # Drop data
      Data.miss <- prodNA(xdata, noNA = prop_miss)
      completedcase <- nrow(na.omit(Data.miss))
      
      ## Apply "Mice" imputation method on dataframe with missing data
      #You may skip removal of collinear variables at initialization by the remove_collinear flag.
      imputed_Data <- mice(Data.miss, method = 'pmm', diagnostics = FALSE , remove_collinear = FALSE, seed = 500)
      mice_imputedData <- as.data.frame(complete(imputed_Data))
      as.matrix(mice_imputedData)
     
      # Apply GAMs to the full dataframe after mice imputation
      Mice_mod <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=mice_imputedData,method="REML")
     
      # GAMs prediction for each xi
      x_pred <- predict.gam(Mice_mod,type="terms",as.data.frame(xdata))
      
      # Rescale fx_i's
      re_fx <- rescale(fx_Data, x_pred)
      for (m in 1:ncol(re_fx)){
        re_fx[, m] <- c(re_fx[, m] - rep(colMeans(re_fx)[m], r))
      }
      
      # Compute root mean squared error (rmse) for each xi's
      x_rmse <- sqrt(colSums((x_pred-re_fx)^2)/r)
      
      # GAMs prediction for y
      Y_pred <- predict.gam(Mice_mod,as.data.frame(xdata))
      # Compute root mean squared error (rmse) for y
      Y_rmse <- c(sqrt(sum((Y_pred-EY)^2)/r))
      
      # Add the above computed x_rmse, y_rmse and the number of complete cases after dropping data
      # into the matrix previously created (eg. loop 1 results go into the first row of the matrix)
      df_mice[i,] <- c(Y_rmse, x_rmse, completedcase)
    }
    
    # Add column names for each column in the dataframe "df_mice"
    colnames(df_mice) <- c("Y_rmse", "x1_rmse", "x2_rmse", "x3_rmse", "x4_rmse", "x5_rmse", "x6_rmse", "# Complete case")
    
    return(df_mice)
  }
  
  # Compute root mean squared error for y variable and each x variables after imputing the missing data
  # If the "rain forest" imputation method was chosen
  else if(method == "rf"){
    for (i in 1:loop){
      # drop data
      Data.miss <- prodNA(xdata, noNA = prop_miss)
      completedcase <- nrow(na.omit(Data.miss))
      
      ## Apply "Random Forest" imputation method on dataframe with missing data
      rf_imputedData <- rfImpute(Y ~ ., Data.miss)[,2:n]
      as.matrix(rf_imputedData)
      
      # Apply GAMs to the full dataframe after mice imputation
      rf_mod <- gam(Y~s(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6),data=rf_imputedData,method="REML")

      # GAMs prediction for each xi
      x_pred <- predict.gam(rf_mod,type="terms",as.data.frame(xdata))
      
      # Rescale fx_i's
      re_fx <- rescale(fx_Data, x_pred)
      for (m in 1:ncol(re_fx)){
        re_fx[, m] <- c(re_fx[, m] - rep(colMeans(re_fx)[m], r))
      }
      
      # Compute root mean squared error (rmse) for each xi's
      x_rmse <- sqrt(colSums((x_pred-re_fx)^2)/r)
      
      # GAMs prediction for y
      Y_pred <- predict.gam(rf_mod,as.data.frame(xdata))
      # Compute root mean squared error (rmse) for y
      Y_rmse <- c(sqrt(sum((Y_pred-EY)^2)/r))
      
      # Add the above computed x_rmse, y_rmse and the number of complete cases after dropping data
      # into the matrix previously created (eg. loop 1 results go into the first row of the matrix)
      df_rf[i,] <- c(Y_rmse, x_rmse,completedcase)
    }
    return(df_rf)
  }
  else{
    return()
  }
}

