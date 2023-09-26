summarize <- function(x){
  y=rbind(x[[1]],x[[2]],x[[3]])
  z=data.frame(cbind(apply(y,2,mean,na.rm=TRUE),apply(y,2,sd,na.rm=TRUE),t(apply(y,2,quantile,c(.5,.025,.975),na.rm=TRUE))))
  names(z)=c("mean","sd","median","lower","upper")
  return(z)
}

beta.moments <- function(mu,sigma){
	alpha = (mu^2-mu^3-mu*sigma^2)/sigma^2
	beta = (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/(sigma^2)
	return(list(alpha = alpha, beta = beta))
}

gamma.moments <- function(mu,sigma){
	alpha <- (mu^2)/(sigma^2)
	beta <- mu/(sigma^2)
	return(list(alpha=alpha,beta=beta))
}