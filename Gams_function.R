#simulation
gsim <- function (eg = 1, n = 400, dist = "normal", scale = 2) {
  if (eg == 1) {
    x1 <- runif(n, 0, 1)
    x2 <- x1 * 0.7 + runif(n, 0, 0.3)
    x3 <- runif(n, 0, 1)
    x4 <- 0.5*x3^2+ runif(n, 0, 0.5)
    x5 <- runif(n, 0, 1)
    x6 <- 0.6*x5^4+ runif(n, 0, 0.4)
    
    f1 <- function(x) x
    f2 <- function(x) 2 * sin(pi * x)
    f3 <- function(x) exp(2 * x)
    f4 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
      (10 * x)^3 * (1 - x)^10
    f5 <- function(x) x^2
    f6 <- function(x) x^3
    f <- f1(x1) + f2(x2) + f3(x3) + f4(x4) + f5(x5) + f6(x6)
    
    if (dist == "normal") {
      e <- rnorm(n, 0, scale)
      y <- f + e
    }
    else if (dist == "poisson") {
      g <- exp(f * scale)
      f <- log(g)
      y <- rpois(rep(1, n), g)
    }
    else if (dist == "binary") {
      f <- (f - 5) * scale
      g <- binomial()$linkinv(f)
      y <- rbinom(g, 1, g)
    }
    else stop("dist not recognised")
    
    data <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6, 
                       f = f, f1 = f1(x1), f2 = f2(x2), f3 = f3(x3), f4 = f4(x4), f5 = f5(x5), f6 = f6(x6))
    return(data)
  }
  
  else if (eg == 2) {
    x1 <- rnorm(n, 0.5, 0.1)
    x2 <- x1 * 0.2+rnorm(n, 0.2,0)*2
    x3 <- rnorm(n, 0.7, 0.1)-0.1
    x4 <- 0.5*(x3-2)^2 + rnorm(n, 5, 3)*0.01 -0.6
    x5 <- rnorm(n, 0.7, 0.1)-0.1
    x6 <- 0.6*x5^4+ rnorm(n, 1.5, 0.4)*0.1
    
    f1 <- function(x) x
    f2 <- function(x) 2 * sin(pi * x)
    f3 <- function(x) exp(2 * x)
    f4 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
      (10 * x)^3 * (1 - x)^10
    f5 <- function(x) x^2
    f6 <- function(x) x^3
    f <- f1(x1) + f2(x2) + f3(x3) + f4(x4) + f5(x5) + f6(x6)
    
    if (dist == "normal") {
      e <- rnorm(n, 0, scale)
      y <- f + e
    }
    else if (dist == "poisson") {
      g <- exp(f * scale)
      f <- log(g)
      y <- rpois(rep(1, n), g)
    }
    else if (dist == "binary") {
      f <- (f - 5) * scale
      g <- binomial()$linkinv(f)
      y <- rbinom(g, 1, g)
    }
    else stop("dist not recognised")
    
    data <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6, 
                       f = f, f1 = f1(x1), f2 = f2(x2), f3 = f3(x3), f4 = f4(x4), f5 = f5(x5), f6 = f6(x6))
    return(data)
  }
  
  else if (eg == 3) {
    x1 <- runif(n, 0, 1)
    x2 <- x1 * 0.7 + runif(n, 0, 0.3)
    x3 <- sin(pi * x2)* 0.7+ runif(n, 0, 0.3)
    x4 <- 0.5*(cos(pi * (2.5*x1))+1)* 0.7+ runif(n, 0, 0.3)
    x5 <- 0.5*(cos(pi * (1.5*x2))+1)* 0.7+ runif(n, 0, 0.3)
    x6 <- rnorm(n, 0.7, 0.1)-0.1
    
    f1 <- function(x) 10*(x-0.5)^2
    f2 <- function(x) 2.5*x
    f3 <- function(x) (0.2 * x^11 * (10 * (1 - x))^6)*2.5/3
    f4 <- function(x) 1/4*(0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10)
    f5 <- function(x) 2*sin(pi * x)
    f6 <- function(x) 4*x^3
    f <- f1(x1) + f2(x2) + f3(x3) + f4(x4) + f5(x5) + f6(x6)
    
    if (dist == "normal") {
      e <- rnorm(n, 0, scale)
      y <- f + e
    }
    else if (dist == "poisson") {
      g <- exp(f * scale)
      f <- log(g)
      y <- rpois(rep(1, n), g)
    }
    else if (dist == "binary") {
      f <- (f - 5) * scale
      g <- binomial()$linkinv(f)
      y <- rbinom(g, 1, g)
    }
    else stop("dist not recognised")
    
    data <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6, 
                       f = f, f1 = f1(x1), f2 = f2(x2), f3 = f3(x3), f4 = f4(x4), f5 = f5(x5), f6 = f6(x6))
    return(data)
  }
  else{
    
  }
}

#f4 <- f4-sum(f4)/300
#c2 <- c(f4/diff(f4)*diff(TMean[,4]))
#DTMean <- sqrt((sum(TMean[,4]-c2)^2)/300)
#matplot(xdata[,4], cbind(c2,TMean[,4]),col=c("red","green"),lty=c(1,1))

Test_rmse <- function(Data, method, prop_miss, loop, graph = "T") {
  
  # Define Y
  y <- Data[,1]
  # Number of columns
  # 1 Y, n xi's and n fx_i's (2n+1 columns in total)
  N <- ncol(Data)
  # Number of observations of each variable
  r <- nrow(Data)
  # XDATA
  w <-N/2
  xdata <- Data[,2:w]
  # Define all columns that containing fx_i's data
  w1 <- w+2
  fx_Data <- Data[,w1:N]
  w2 <- w+1
  Ey <- Data[,w2]
  # Number of coulmns in dataframe "Data"
  n <- ncol(xdata)
  #x data with y
  data <- Data[,1:w]
  #gams formual
  p<-n
  form <- as.formula(paste0("y~",paste0("s(x",1:p,")",collapse="+")))
  
  # Define matrix for storing Y_rmse, x_rmse and number of complete case
  df_result <- data.frame(matrix(ncol = n+2, nrow = loop))
  
  #graph data from loop
  sx_imputed_pred1 <- data.frame(matrix(0,ncol = ncol(fx_Data), nrow = nrow(fx_Data)))
  Y_imputed_pred1 <- rep(0,length(y))
  
  
  # Apply GAMs to the full dataframe after mean imputation
  mod <- gam(formula=form,data=as.data.frame(data),method="REML")
  
  # Apply GAMs to the non-dropping full dataframe
  Y_pred <- predict.gam(mod)
  
  # GAMs prediction for each xi
  sx_pred <- predict.gam(mod,type="terms")  
  
  # Rescale fx_i's,x_pred
  for (m1 in 1:ncol(fx_Data)){
    fx_Data[, m1] <- fx_Data[, m1] - colMeans(fx_Data)[m1]
  }
  for (m2 in 1:ncol(sx_pred)){
    sx_pred[, m2] <- sx_pred[, m2] - colMeans(sx_pred)[m2]
  }
  
  if (graph == "T"){
    df1 <- data.frame(matrix(0,ncol = 2, nrow = r))
    df2 <- data.frame(matrix(0,ncol = 2, nrow = r))
    par(mfrow=c(2,3))
    for (k in 1:ncol(fx_Data)){
      df1[,1] <- xdata[,k]
      df1[,2] <- sx_pred[,k]
      df2[,1] <- xdata[,k]
      df2[,2] <- fx_Data[,k]
      print(ggplot() +
              geom_line(data = df1, aes(x = X1, y = X2, color = "predicted")) +
              geom_line(data = df2, aes(x = X1, y = X2, color = "ture")) +
              scale_color_manual(name = "Lines",
                                 values = c("predicted" = "blue", "ture" = "red")))
    }
  }
  else{
  }
  YO_rmse <- sqrt(sum((Y_pred-Ey)^2)/r)
  
  
  
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
      b <- cbind(y,mean_imputedData)
      
      # Apply GAMs to the full dataframe after mean imputation
      Mean_mod <- gam(formula=form,data=as.data.frame(b),method="REML")
      
      # GAMs prediction for y
      Y_imputed_pred <- predict.gam(Mean_mod)
      
      # GAMs prediction for each xi
      sx_imputed_pred <- predict.gam(Mean_mod,type="terms")  
      
      # Rescale fx_i's,x_pred
      
      for (m in 1:ncol(sx_imputed_pred)){
        sx_imputed_pred[, m] <- sx_imputed_pred[, m] - colMeans(sx_imputed_pred)[m]
      }
      
      # Compute root mean squared error (rmse) for each xi's
      f1 <- function(x) 10*(x-0.5)^2
      f2 <- function(x) 2.5*x
      f3 <- function(x) (0.2 * x^11 * (10 * (1 - x))^6)*2.5/3
      f4 <- function(x) 1/4*(0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10)
      f5 <- function(x) 2*sin(pi * x)
      f6 <- function(x) 4*x^3
      f_imputed_xdata <- data.frame(matrix(0,ncol = ncol(fx_Data), nrow = r))
      
      #for (x in 1:ncol(fx_Data)){
        #f_imputed_xdata[,x] <- fxn[x](mean_imputedData[,x])
      #}
      f_imputed_xdata[,1] <- f1(mean_imputedData[,1])
      f_imputed_xdata[,2] <- f2(mean_imputedData[,2])
      f_imputed_xdata[,3] <- f3(mean_imputedData[,3])
      f_imputed_xdata[,4] <- f4(mean_imputedData[,4])
      f_imputed_xdata[,5] <- f5(mean_imputedData[,5])
      f_imputed_xdata[,6] <- f6(mean_imputedData[,6])
      
      for (m3 in 1:ncol(f_imputed_xdata)){
        f_imputed_xdata[, m3] <- f_imputed_xdata[, m3] - colMeans(f_imputed_xdata)[m3]
      }
      
      x_rmse <- sqrt(colSums((sx_imputed_pred-f_imputed_xdata)^2)/r)
      
      # Compute root mean squared error (rmse) for y
      Y_rmse <- sqrt(sum((Y_imputed_pred-Ey)^2)/r)
      
      # Add the above computed x_rmse, y_rmse and the number of complete cases after dropping data
      # into the matrix previously created (eg. loop 1 results go into the first row of the matrix)
      df_result[i,] <- c(Y_rmse,x_rmse,completedcase)
      
      # Add column names for each column in the dataframe "df_result"
      colnames(df_result) <- c("Y_rmse", "x1_rmse", "x2_rmse", "x3_rmse", "x4_rmse", "x5_rmse", "x6_rmse", "# Complete case")
      
      # in the loop graph data
      if (graph == "T"){
        df1 <- data.frame(matrix(0,ncol = 2, nrow = r))
        df2 <- data.frame(matrix(0,ncol = 2, nrow = r))
        par(mfrow=c(2,3))
        for (k in 1:ncol(fx_Data)){
          df1[,1] <- mean_imputedData[,k]
          df1[,2] <- sx_imputed_pred[,k]
          df2[,1] <- xdata[,k]
          df2[,2] <- fx_Data[,k]
          print(ggplot() +
                  geom_line(data = df1, aes(x = X1, y = X2, color = "predicted")) +
                  geom_line(data = df2, aes(x = X1, y = X2, color = "ture")) +
                  scale_color_manual(name = "Lines",
                                     values = c("predicted" = "blue", "ture" = "red")))
        }
        sx_imputed_pred1 <- sx_imputed_pred1 + sx_imputed_pred
        Y_imputed_pred1 <- Y_imputed_pred1 + Y_imputed_pred
      }
      else{
      }
    }
    
    # out of loop graph
    if (graph == "T"){
      df1 <- data.frame(matrix(0,ncol = 2, nrow = r))
      df2 <- data.frame(matrix(0,ncol = 2, nrow = r))
      sx_imputed_pred1 <- sx_imputed_pred1/loop
      par(mfrow=c(2,3))
      for (k in 1:ncol(fx_Data)){
        df1[,1] <- mean_imputedData[,k]
        df1[,2] <- sx_imputed_pred1[,k]
        df2[,1] <- xdata[,k]
        df2[,2] <- fx_Data[,k]
        print(ggplot() +
                geom_line(data = df1, aes(x = X1, y = X2, color = "predicted")) +
                geom_line(data = df2, aes(x = X1, y = X2, color = "ture")) +
                scale_color_manual(name = "Lines",
                                   values = c("predicted" = "blue", "ture" = "red")))
      }
      ggplot()+
        geom_point(aes(x=Y_imputed_pred1/loop, y=Ey))
    }
    else{
    }
    return(df_result)
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
      b <- cbind(y,mice_imputedData)
      
      # Apply GAMs to the full dataframe after mean imputation
      Mice_mod <- gam(formula=form,data=as.data.frame(b),method="REML")
      
      # GAMs prediction for y
      Y_imputed_pred <- predict.gam(Mice_mod)
      
      # GAMs prediction for each xi
      sx_imputed_pred <- predict.gam(Mice_mod,type="terms")  
      
      # Rescale fx_i's,x_pred
      
      for (m in 1:ncol(sx_imputed_pred)){
        sx_imputed_pred[, m] <- sx_imputed_pred[, m] - colMeans(sx_imputed_pred)[m]
      }
      
      # Compute root mean squared error (rmse) for each xi's
      f1 <- function(x) 10*(x-0.5)^2
      f2 <- function(x) 2.5*x
      f3 <- function(x) (0.2 * x^11 * (10 * (1 - x))^6)*2.5/3
      f4 <- function(x) 1/4*(0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10)
      f5 <- function(x) 2*sin(pi * x)
      f6 <- function(x) 4*x^3
      f_imputed_xdata <- data.frame(matrix(0,ncol = ncol(fx_Data), nrow = r))
      
      #for (x in 1:ncol(fx_Data)){
      #f_imputed_xdata[,x] <- fxn[x](mean_imputedData[,x])
      #}
      f_imputed_xdata[,1] <- f1(mice_imputedData[,1])
      f_imputed_xdata[,2] <- f2(mice_imputedData[,2])
      f_imputed_xdata[,3] <- f3(mice_imputedData[,3])
      f_imputed_xdata[,4] <- f4(mice_imputedData[,4])
      f_imputed_xdata[,5] <- f5(mice_imputedData[,5])
      f_imputed_xdata[,6] <- f6(mice_imputedData[,6])
      
      for (m3 in 1:ncol(f_imputed_xdata)){
        f_imputed_xdata[, m3] <- f_imputed_xdata[, m3] - colMeans(f_imputed_xdata)[m3]
      }
      
      x_rmse <- sqrt(colSums((sx_imputed_pred-f_imputed_xdata)^2)/r)
      
      # Compute root mean squared error (rmse) for y
      Y_rmse <- sqrt(sum((Y_imputed_pred-Ey)^2)/r)
      
      # Add the above computed x_rmse, y_rmse and the number of complete cases after dropping data
      # into the matrix previously created (eg. loop 1 results go into the first row of the matrix)
      df_result[i,] <- c(Y_rmse,x_rmse,completedcase)
      
      # Add column names for each column in the dataframe "df_result"
      colnames(df_result) <- c("Y_rmse", "x1_rmse", "x2_rmse", "x3_rmse", "x4_rmse", "x5_rmse", "x6_rmse", "# Complete case")
      
      # in the loop graph data
      if (graph == "T"){
        df1 <- data.frame(matrix(0,ncol = 2, nrow = r))
        df2 <- data.frame(matrix(0,ncol = 2, nrow = r))
        par(mfrow=c(2,3))
        for (k in 1:ncol(fx_Data)){
          df1[,1] <- mice_imputedData[,k]
          df1[,2] <- sx_imputed_pred[,k]
          df2[,1] <- xdata[,k]
          df2[,2] <- fx_Data[,k]
          print(ggplot() +
                  geom_line(data = df1, aes(x = X1, y = X2, color = "predicted")) +
                  geom_line(data = df2, aes(x = X1, y = X2, color = "ture")) +
                  scale_color_manual(name = "Lines",
                                     values = c("predicted" = "blue", "ture" = "red")))
        }
        sx_imputed_pred1 <- sx_imputed_pred1 + sx_imputed_pred
        Y_imputed_pred1 <- Y_imputed_pred1 + Y_imputed_pred
      }
      else{
      }
    }
    
    # out of loop graph
    if (graph == "T"){
      df1 <- data.frame(matrix(0,ncol = 2, nrow = r))
      df2 <- data.frame(matrix(0,ncol = 2, nrow = r))
      sx_imputed_pred1 <- sx_imputed_pred1/loop
      par(mfrow=c(2,3))
      for (k in 1:ncol(fx_Data)){
        df1[,1] <- mice_imputedData[,k]
        df1[,2] <- sx_imputed_pred1[,k]
        df2[,1] <- xdata[,k]
        df2[,2] <- fx_Data[,k]
        print(ggplot() +
                geom_line(data = df1, aes(x = X1, y = X2, color = "predicted")) +
                geom_line(data = df2, aes(x = X1, y = X2, color = "ture")) +
                scale_color_manual(name = "Lines",
                                   values = c("predicted" = "blue", "ture" = "red")))
      }
      ggplot()+
        geom_point(aes(x=Y_imputed_pred1/loop, y=Ey))
    }
    else{
    }
    return(df_result)
  }
  
  # Compute root mean squared error for y variable and each x variables after imputing the missing data
  # If the "rain forest" imputation method was chosen
  else if(method == "rf"){
    for (i in 1:loop){
      # drop data
      Data.miss <- prodNA(xdata, noNA = prop_miss)
      completedcase <- nrow(na.omit(Data.miss))
      
      ## Apply "Random Forest" imputation method on dataframe with missing data
      rf_imputedData <- rfImpute(y ~ ., Data.miss)
      as.data.frame(rf_imputedData)
      rf_imputedData_x <- rf_imputedData[,2:(length(fx_Data)+1)]
      
      # Apply GAMs to the full dataframe after mean imputation
      rf_mod <- gam(formula=form,data=as.data.frame(rf_imputedData),method="REML")
      
      # GAMs prediction for y
      Y_imputed_pred <- predict.gam(rf_mod)
      
      # GAMs prediction for each xi
      sx_imputed_pred <- predict.gam(rf_mod,type="terms")  
      
      # Rescale fx_i's,x_pred
      
      for (m in 1:ncol(sx_imputed_pred)){
        sx_imputed_pred[, m] <- sx_imputed_pred[, m] - colMeans(sx_imputed_pred)[m]
      }
      
      # Compute root mean squared error (rmse) for each xi's
      f1 <- function(x) 10*(x-0.5)^2
      f2 <- function(x) 2.5*x
      f3 <- function(x) (0.2 * x^11 * (10 * (1 - x))^6)*2.5/3
      f4 <- function(x) 1/4*(0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10)
      f5 <- function(x) 2*sin(pi * x)
      f6 <- function(x) 4*x^3
      f_imputed_xdata <- data.frame(matrix(0,ncol = ncol(fx_Data), nrow = r))
      
      #for (x in 1:ncol(fx_Data)){
      #f_imputed_xdata[,x] <- fxn[x](mean_imputedData[,x])
      #}
      f_imputed_xdata[,1] <- f1(rf_imputedData_x[,1])
      f_imputed_xdata[,2] <- f2(rf_imputedData_x[,2])
      f_imputed_xdata[,3] <- f3(rf_imputedData_x[,3])
      f_imputed_xdata[,4] <- f4(rf_imputedData_x[,4])
      f_imputed_xdata[,5] <- f5(rf_imputedData_x[,5])
      f_imputed_xdata[,6] <- f6(rf_imputedData_x[,6])
      
      for (m3 in 1:ncol(f_imputed_xdata)){
        f_imputed_xdata[, m3] <- f_imputed_xdata[, m3] - colMeans(f_imputed_xdata)[m3]
      }
      
      x_rmse <- sqrt(colSums((sx_imputed_pred-f_imputed_xdata)^2)/r)
      
      # Compute root mean squared error (rmse) for y
      Y_rmse <- sqrt(sum((Y_imputed_pred-Ey)^2)/r)
      
      # Add the above computed x_rmse, y_rmse and the number of complete cases after dropping data
      # into the matrix previously created (eg. loop 1 results go into the first row of the matrix)
      df_result[i,] <- c(Y_rmse,x_rmse,completedcase)
      
      # Add column names for each column in the dataframe "df_result"
      colnames(df_result) <- c("Y_rmse", "x1_rmse", "x2_rmse", "x3_rmse", "x4_rmse", "x5_rmse", "x6_rmse", "# Complete case")
      
      # in the loop graph data
      if (graph == "T"){
        df1 <- data.frame(matrix(0,ncol = 2, nrow = r))
        df2 <- data.frame(matrix(0,ncol = 2, nrow = r))
        par(mfrow=c(2,3))
        for (k in 1:ncol(fx_Data)){
          df1[,1] <- rf_imputedData_x[,k]
          df1[,2] <- sx_imputed_pred[,k]
          df2[,1] <- xdata[,k]
          df2[,2] <- fx_Data[,k]
          print(ggplot() +
                  geom_line(data = df1, aes(x = X1, y = X2, color = "predicted")) +
                  geom_line(data = df2, aes(x = X1, y = X2, color = "ture")) +
                  scale_color_manual(name = "Lines",
                                     values = c("predicted" = "blue", "ture" = "red")))
        }
        sx_imputed_pred1 <- sx_imputed_pred1 + sx_imputed_pred
        Y_imputed_pred1 <- Y_imputed_pred1 + Y_imputed_pred
      }
      else{
      }
    }
    
    # out of loop graph
    if (graph == "T"){
      df1 <- data.frame(matrix(0,ncol = 2, nrow = r))
      df2 <- data.frame(matrix(0,ncol = 2, nrow = r))
      sx_imputed_pred1 <- sx_imputed_pred1/loop
      par(mfrow=c(2,3))
      for (k in 1:ncol(fx_Data)){
        df1[,1] <- rf_imputedData_x[,k]
        df1[,2] <- sx_imputed_pred1[,k]
        df2[,1] <- xdata[,k]
        df2[,2] <- fx_Data[,k]
        print(ggplot() +
                geom_line(data = df1, aes(x = X1, y = X2, color = "predicted")) +
                geom_line(data = df2, aes(x = X1, y = X2, color = "ture")) +
                scale_color_manual(name = "Lines",
                                   values = c("predicted" = "blue", "ture" = "red")))
      }
      ggplot()+
        geom_point(aes(x=Y_imputed_pred1/loop, y=Ey))
    }
    else{
    }
    return(df_result)
  }
  
  # else
  else{
    
  }
}

