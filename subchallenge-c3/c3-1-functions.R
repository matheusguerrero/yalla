require(mvnfast)
require(this.path)
setwd(this.path::here())

# CDF, quantile function, Jacobian for transformation to delta-Laplace margins.

#------------------------------------------------------------------------------------------------------------------
# pdlaplace: univariate distribution function


pdlaplace<-function(z,mu,sigma,delta,lower.T=T)
{
  
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  nz<-length(z)
  result<-numeric(nz)
  if(length(delta)==1){delta<-rep(delta,nz)}
  if(lower.T==T){
    result[z<mu]<- 0.5*pgamma((((mu-z)/sigma)[z<mu])^delta[z<mu],shape=1/delta[z<mu],scale=1,lower.tail = F)
    result[z>=mu]<- 0.5+0.5*pgamma((((z-mu)/sigma)[z>=mu])^delta[z>=mu],shape=1/delta[z>=mu],scale=1)
  }else{
    
    result[z<mu]<- 1-0.5*pgamma((((mu-z)/sigma)[z<mu])^delta[z<mu],shape=1/delta[z<mu],scale=1,lower.tail = F)
    result[z>=mu]<- 0.5*pgamma((((z-mu)/sigma)[z>=mu])^delta[z>=mu],shape=1/delta[z>=mu],scale=1,lower.tail=F)
  }
  
  return(result)
}

#------------------------------------------------------------------------------------------------------------------
# qdlaplace: univariate quantile function

qdlaplace<-function(q,mu,sigma,delta)
{
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  nq<-length(q)
  result<-numeric(nq)
  if(length(mu)==1){mu<-rep(mu,nq)}
  if(length(sigma)==1){sigma<-rep(sigma,nq)}
  if(length(delta)==1){delta<-rep(delta,nq)}
  result[q<0.5]<- mu[q<0.5] - sigma[q<0.5]*qgamma(1-2*q[q<0.5],shape=1/delta[q<0.5],scale=1)^(1/delta[q<0.5])
  result[q>=0.5]<- mu[q>=0.5] + sigma[q>=0.5]*qgamma(2*q[q>=0.5]-1,shape=1/delta[q>=0.5],scale=1)^(1/delta[q>=0.5])
  return(result)
}


#------------------------------------------------------------------------------------------------------------------
# rdlaplace: univariate random number generation

rdlaplace<-function(n,mu,sigma,delta)
{
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  return(qdlaplace(runif(n),mu=mu,delta=delta,sigma=sigma))
}


#------------------------------------------------------------------------------------------------------------------
# ddlaplace: univariate density function

ddlaplace<-function(z,mu,sigma,delta,log=FALSE)
{
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  ld<- -abs((z-mu)/(sigma))^delta +log(delta)-log(2*sigma)-lgamma(1/delta)
  if(log==TRUE){return(ld)}
  else{return(exp(ld))}
}

#------------------------------------------------------------------------------------------------------------------
# dmvdlaplace: multivariate density function with Gaussian copula


# with faster MVN computation: requires matrix Sigma AND its cholesky factorization, chol(Sigma), which results in faster
# computation when used in apply()

dmvdlaplace<-function(z,mu,sigmad,SigmaChol,Sigma,delta,log=FALSE)
{
  zn<-qnorm(pdlaplace(z,mu=mu,sigma=sigmad,delta=delta,lower.T=F),mean=mu,sd=sqrt(diag(Sigma)),lower.tail =F)
  ld1<-dmvn(zn,mu=mu,sigma = SigmaChol,log=TRUE,isChol=TRUE)
  ld2<-sum(ddlaplace(z,mu=mu,sigma=sigmad,delta=delta,log=TRUE))-sum(dnorm(zn, mean=mu, sd=sqrt(diag(Sigma)), log=TRUE))
  if(log==TRUE){return(ld1+ld2)}
  else{return(exp(ld1+ld2))}
}

#------------------------------------------------------------------------------------------------------------------
# rmvdlaplace: multivariate random number generation with Gaussian copula

rmvdlaplace<-function(n,dim,mu,sigmad,Sigma,delta)
{
  if(dim(Sigma)[1]!=dim){stop("dim and Sigma do not match")}
  NU<-t(apply(rmvn(n,mu=rep(0,length=dim),sigma=Sigma),1,pnorm,sd=sqrt(diag(Sigma))))
  return(t(apply(NU,1,qdlaplace,mu=mu,sigma=sigmad,delta=delta)))
}

# Heffernan and Tawn nll

nll=function(par,X0,X1X2){
  
 alpha.vec=par[1:2]
 beta.vec=par[3:4]
 sig.vec=par[5:6]
 mu.vec=par[7:8]
 delta.vec=par[9:10]
 rho=par[11]
  
    #Add in any constraints here
 if(any(alpha.vec< -1) | any(alpha.vec > 1 ) | any(beta.vec> 1 ) | any(beta.vec < 0) | any(sig.vec<0) | any(delta.vec < 0) )  return(1e10)
 
 cov.mat=diag(sig.vec)%*%matrix(c(1,rho,rho,1),2,2)%*%diag(sig.vec)

      nll<-apply(cbind(X0,X1X2),1,function(x){
        mean=alpha.vec*x[1]+(x[1]^beta.vec)*mu.vec
        var=diag(c(x[1]^beta.vec))%*%cov.mat%*%diag(c(x[1]^beta.vec))
        dmvdlaplace(x[2:3]-mean,mu= c(0,0),sigmad=sqrt(diag(var)),SigmaChol=chol(var),Sigma=var,
                    delta=delta.vec,log=T)
      })
    
   
  
  if(is.finite(sum(nll))){
    return(-sum(nll))
  }else{return(1e10)}
  
  
}


data<-read.csv("../Data/Coputopia.csv")

X=as.matrix(data[,3:5])


X0=X[,1]

X1X2=X[,2:3]

plot(X0,X1X2[,1])

u=quantile(X0,0.9)

X0.exceed=X0[X0>u]
X1X2.exceed=X1X2[X0>u,]

points(X0.exceed,X1X2.exceed[,1],col="red")

init.par=rep(0.5,11)
fit1<-optim(par=init.par,fn=nll,X0=X0.exceed,X1X2=X1X2.exceed)

