library(readr);library(glasso);library(glmnet);library(CVglasso);library(pracma)
#library(xlsx)


##### Five Functions Required for Multivariate Adaptive Lasso Based on Lasso###############

##(1)Coordinate Descent Function———————————————————————————————————————
coor_discent=function(x,y,w,lambda,beta){
  
  ##   This function performs coordinate descent optimization for multivariate regression.
  ##   
  ##   Args:
  ##     x: The nxp matrix of predictor variables.
  ##     y: The nxs matrix of response variables.
  ##     w: The weight matrix.
  ##     lambda: The regularization parameter.
  ##     beta: The initial coefficient matrix.
  ##
  ##   Returns:
  ##     Updated coefficient matrix after performing coordinate descent optimization.

  
  
  p=ncol(x)
  s=ncol(y)
  a = rep(0,nrow(x))
  b = rep(0,nrow(x))
  
  for(m in 1:p){
    for(n in 1:s){
      for(i in 1:nrow(x)){
        a[i] = x[i, m]*(y[i,n]-x[i, ]%*%beta[, n]+x[i, m]*beta[m, n])
        b[i] = x[i, m]*x[i, m]
      }
      z = (sum(a))/nrow(x)
      r = lambda*w[m,n]
      q = sum(b)/nrow(x)
      if (r<abs(z)){
        if(z>0) beta[m, n] = (z-r)/q
        else if(z<0) beta[m, n] = (z+r)/q
        else beta[m, n] = 0
      }
      else beta[m, n] = 0
    }
  }
  return(beta)
}


###(2) Coordinate Descent Function with Diagonal Elements of Precision Matrix——————————————————————————————————————————
diag_coordinate = function(x,y,w,lambda,beta,precision_diag){
  ###   This function performs coordinate descent optimization for multivariate regression with diagonal elements of the precision matrix.
  ##   
  ###   Args:
  ###     x: The nxp matrix of predictor variables.
  ###     y: The nxs matrix of response variables.
  ###     w: The weight matrix.
  ###     lambda: The regularization parameter.
  ###     beta: The initial coefficient matrix.
  ###     precision_diag: A vector containing diagonal elements of the precision matrix.
  ##
  ###   Returns:
  ###     Updated coefficient matrix after performing coordinate descent optimization.
  ###     
  ###   Notes:
  ###     This function is similar to the coordinate descent function but considers diagonal elements of the precision matrix
  ###     in calculating the regularization threshold.
  
  p=ncol(x)
  s=ncol(y)
  a = rep(0,nrow(x))
  b = rep(0,nrow(x))
  
  for(m in 1:p){
    for(n in 1:s){
      for(i in 1:nrow(x)){
        a[i] = (x[i, m]*(y[i,n]-x[i, ]%*%beta[, n]+x[i, m]*beta[m, n]))
        b[i] = (x[i, m]*x[i, m])
      }
      z = (sum(a))
      r = lambda*w[m,n]*sqrt(nrow(x))/(2*precision_diag[n])
      q = sum(b)
      if (r<abs(z)){
        if(z>0) beta[m, n] = (z-r)/q
        else if(z<0) beta[m, n] = (z+r)/q
        else beta[m, n] = 0
      }
      else beta[m, n] = 0
    }
  }
  return(beta)
}


##(3)Multivariate BIC Function——————————————————————————————————————————————
multi_bic = function(x,y,beta){ 
  # x is an n*p matrix, y is an n*s matrix, beta is a p*s matrix
  
  s=ncol(y)
  if(is.null(s)){s = 1}
  bic=rep(0,s)
  n = nrow(x)
  for (i in 1:s){  
    sk=sum(beta[,i]!=0) # Count of non-zero elements in the ith column of beta
    RSS=sum((y[,i]-x%*%beta[,i])^2) # Residual Sum of Squares (RSS)
    bic[i]=n*log(RSS)+log(n)*sk # BIC calculation for the ith response variable
  }
  bic=sum(bic)
  return(bic)
}


##(4)Multivariate Adaptive Lasso Function——————————————————————————————————————————
mAdaLas = function(x,y,lambdamax=10,lambdamin=0.0001)
{
  ##   This function performs Multivariate Adaptive Lasso based on a sequence of lambda values.
  ##
  ##   Args:
  ##     x: The nxp matrix of predictor variables.
  ##     y: The nxs matrix of response variables.
  ##     lambdamax: The maximum value of the regularization parameter lambda.
  ##     lambdamin: The minimum value of the regularization parameter lambda.
  ##
  ##   Returns:
  ##     A list containing various information including BIC values, estimated coefficients, and timings.
  ##     
  ##   Notes:
  ##     The function performs the following steps:
  ##     1. Calculate Lasso estimates using cross-validation (`cv.glmnet`) for each response variable.
  ##     2. Form a weight matrix based on the absolute values of the Lasso estimates.
  ##     3. Generate a sequence of lambda values and iterate through them.
  ##     4. Perform coordinate descent optimization (`coor_discent`) for each lambda value.
  ##     5. Calculate the BIC for each estimate and choose the lambda value with the lowest BIC.
  
  time1 = Sys.time()
  n = nrow(x); p = ncol(x);s = ncol(y) # Matrix dimensions
  s2 = cv.glmnet(x, y,family = "mgaussian")  # Cross-validation for glmnet
  time2 = Sys.time()
  
  beta = matrix(0,p,s)
  for (i in 1:s) {
    beta[,i] = s2$glmnet.fit$beta[[i]][,which.min(s2$cvm)]
  } 
  
  
  betalasso=beta  # Lasso estimates from the first step
  
  # Forming the weight matrix
  w=abs(beta)
  beta0=rep(0,1);betanot0=rep(0,1)
  for (i in 1:p){
    for (j in 1:s){
      if (w[i,j]==0)
        w[i,j]=1/(n)
    }
  }
  w=1/w
  
  time3=Sys.time()
  
  # Generate a sequence of lambda values
  tmax=log(lambdamax);tmin=log(lambdamin)
  tseq=seq(tmax,tmin,-0.6)
  lambdatest=exp(tseq)
  nn=length(lambdatest)
  betaall=array(0,c(p,s,nn))
  
  # Calculate BIC for each lambda and find the best one
  bic=rep(0,nn)
  
  for (i in 1:nn){
    betaall[,,i]=coor_discent(x,y,w,lambdatest[i],betalasso)  # Using coordinate descent
    bic[i]=multi_bic(x,y,betaall[,,i]) # Calculate BIC for current estimate
  }
  bestlambda=which.min(bic)
  final_beta=betaall[,,bestlambda]
  time4=Sys.time()
  #}
  out = list(i=i,bic=bic,betalasso = betalasso,betafinal = final_beta,
             cordinateTime = time2 - time1,cordinateTime2 = time4 - time3,
             lambada =lambdatest[bestlambda],nn=nn)
}


##(5)Multivariate Weighted Lasso Function with Diagonal Precision Matrix Elements————————————————————————————————————
d_m_ada_lasso = function(x,y,precision_diag,lambdamax=100,lambdamin=0.001)
{
  ## d_m_ada_lasso <- function(x, y, precision_diag, lambdamax = 100, lambdamin = 0.001) {
  ##   This function performs calibrated Multivariate Weighted Lasso(MWL) with diagonal elements of the precision matrix.
  ##
  ##   Args:
  ##     x: The nxp matrix of predictor variables.
  ##     y: The nxs matrix of response variables.
  ##     precision_diag: A vector containing diagonal elements of the precision matrix.
  ##     lambdamax: The maximum value of the regularization parameter lambda.
  ##     lambdamin: The minimum value of the regularization parameter lambda.
  ##
  ##   Returns:
  ##     A list containing various information including BIC values, estimated coefficients, and the chosen lambda.
  ##     
  ##   Notes:
  ##     The function performs the following steps:
  ##     1. Calculate Lasso estimates using cross-validation (`cv.glmnet`) for each response variable.
  ##     2. Form a weight matrix based on the absolute values of the Lasso estimates.
  ##     3. Generate a sequence of lambda values and iterate through them.
  ##     4. Perform calibrated coordinate descent optimization (`diag_coordinate`) with diagonal precision elements for each lambda value.
  ##     5. Calculate the BIC for each estimate and choose the lambda value with the lowest BIC.
  
  
  n = nrow(x); p = ncol(x);s = ncol(y) # Matrix dimensions
  s2 = cv.glmnet(x, y,family = "mgaussian") # Cross-validation for glmnet
  
  
  beta = matrix(0,p,s)
  for (i in 1:s) {
    beta[,i] = s2$glmnet.fit$beta[[i]][,which.min(s2$cvm)] 
  }
  betalasso=beta # Lasso estimates from the first step
  
  # Forming the weight matrix
  w=abs(beta)
  beta0=rep(0,1);betanot0=rep(0,1)
  for (i in 1:p){
    for (j in 1:s){
      if (w[i,j]==0)
        w[i,j]=1/(n)
    }
  }
  w=1/w
  
  # Generate a sequence of lambda values
  tmax=log(lambdamax);tmin=log(lambdamin)
  tseq=seq(tmax,tmin,-0.6)
  lambdatest=exp(tseq)
  nn=length(lambdatest)
  betaall=array(0,c(p,s,nn))
  
  # Calculate BIC for each lambda and find the best one
  bic=rep(0,nn)
  
  
  for (i in 1:nn){
    betaall[,,i]=diag_coordinate(x,y,w,lambdatest[i],betalasso,precision_diag)
    bic[i]=multi_bic(x,y,betaall[,,i])
  }
  bestlambda=which.min(bic)
  final_beta=betaall[,,bestlambda]
  
  out = list(bic=bic,betalasso = betalasso,betafinal = final_beta,
             lambada =lambdatest[bestlambda],nn=nn)
}




