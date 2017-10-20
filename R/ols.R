#Written by: Radojevic, Marco; Besch, Johannes
rm(list=ls())
ols <- function(X,Y,data){
  rownames <- as.matrix(colnames(data))
  beta.hat = solve(t(X) %*% X) %*% t(X) %*% Y
  betah <- t(beta.hat)
  yhat <- betah %*% t(X)
  yhat <- t(yhat)
  resi <- Y-yhat
  # calculate residual quartils
  resq = quantile(resi,p=seq(0,1,length.out=5))
  names(resq) = c("Min", "1Q", "Median", "3Q", "Max")
  n = nrow(Y)
  k = ncol(X)
  df = n-k
  vcov = 1 / (n-k) * as.numeric(t(resi) %*% resi) * solve(t(X) %*% X)   # calculate variance-covariance-matrix
  s.e. = sqrt(diag(vcov))   # calculate the standard errors
  t = as.vector(betah) / s.e.   # calculate t value and add them to the result
  p = 2*pt(abs(t), df=n-k,lower.tail= FALSE)# calculate p value and add them to the result
  stars = ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,"."," "))))
  residual.error = sqrt(sum((yhat-Y)^2)/(n-k))
  multiplersq = 1 - sum((resi)^2) / sum((Y-mean(Y))^2) 
  adjustedr = 1 - (1 - multiplersq)*(n-1) / (n-k)  # calculate adjusted rsquared
  # calculate F-statistics and its p-value
  f =  multiplersq/(k - 1) / ((1 - multiplersq) / (n-k))
  f.p = pf(f, k-1, n-k, lower.tail=F)
  o <- k-1
  s.e. <- as.matrix(s.e.)
  t <- as.matrix(t)
  p <- as.matrix(p)
  rownames[1] <- c("Intercept") 
  result = as.data.frame(cbind(rownames,beta.hat,s.e.,t,p,stars))#  result.names
  colnames(result) = c("Coeff.","Estimate","Std. Error","t","Pr(>|t|)","")
  
  writeLines("Residuals:")
  print(resq)
  writeLines(" ")
  print(result,row.names=F)
  cat(paste("---", 
            "Signif. codes:   0 *** 0.001 ** 0.01 * 0.05 . 0.1  1",
            sep="\n"))
  cat('\n\n', 'Residual standard error:', residual.error, 'on', df, 'degrees of freedom', '\n',
      'Multiple R-squared:', multiplersq, 'Adjusted R-squared:', adjustedr, '\n',
      'F-statistic:', f, 'on', df, 'DF',  'p-value:', f.p)

}
