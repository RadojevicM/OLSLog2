logit <- function(fml,dat) {
  loglik<- function(theta, y, x) {
    y <- y
    x <- as.matrix(x)
    beta   <- theta[1:ncol(x)]               
    ll <- sum(-y*log(1+exp(-(x%*%beta)))-(1-y)*log(1+exp(x%*%beta)))
    return(-ll)
  }
  out <- rownames(attr(terms(fml),"factors"))[1]
  tempdf <- model.frame(dat)
  x <- as.matrix(model.matrix(fml,data=tempdf))
  y <- as.numeric(as.matrix(dat[,match(out,colnames(dat))]))
  theta.ini <- rep(0,(dim(x)[2]))
  names(theta.ini) <- colnames(x)
  mle <- optim(theta.ini,loglik,x=x,y=y,hessian=T,method = "BFGS")
  output <- list(beta=mle$par,s.e.=sqrt(diag(solve(mle$hessian))),ll=2*mle$value)
  return(output)
}



