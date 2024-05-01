


library(rmutil)
u=qlaplace(0.9)

betas=read.csv("Data/all_beta.csv")
input=read.csv("Data/Model1_alpha_res.csv")
input2=read.csv("Data/Model2_alpha_res.csv")
input3=read.csv("Data/Model3_alpha_res.csv")

pars=vector("list",3)
pars[[1]]=list("beta"=t(as.matrix(betas$Model1)),"alpha.vec"=cbind(input$Model1_a1,input$Model1_a2),"resid"=cbind(input$res_z1,input$res_z2))
pars[[2]]=list("beta"=t(as.matrix(betas$Model2)),"alpha.vec"=cbind(input2$Model2_a1,input2$Model2_a2),"resid"=cbind(input2$res_z1,input2$res_z2))
pars[[3]]=list("beta"=t(as.matrix(betas$Model3)),"alpha.vec"=cbind(input3$Model3_a1,input3$Model2_a2),"resid"=cbind(input3$res_z1,input3$res_z2))

rsim=function(beta,alphas,resid,v,cond.ind){
  
  X0=v+rexp(1)
  
  alpha.ind=sample(1:nrow(alphas),1)
  resid.ind=sample(1:nrow(resid),1)
  
  X1X2=(alphas[alpha.ind,])*X0+X0^beta*resid[resid.ind,]
  
  if(cond.ind==1) return(c(X0,X1X2))
  if(cond.ind==2) return(c(X1X2[1],X0,X1X2[2]))
  if(cond.ind==3) return(c(X1X2,X0))
  
}


rsim_all=function(pars,n,B,v){
  
  cond.inds=as.matrix(sample(1:3,B,replace=T))
  
  out<-t(apply(cond.inds,1,function(x){
    rsim(pars[[x]]$beta,pars[[x]]$alpha.vec,pars[[x]]$resid,v,x)
    
  }))
  
import.prob=apply(out,1,function(x){
    1/sum(x>v)
  })
import.prob=import.prob/sum(import.prob)

samp.inds = sample(1:length(import.prob),size=n,prob=import.prob)

  return(out[samp.inds,])
}


data<-read.csv("Data/Coputopia.csv")

X=as.matrix(data[,3:5])


X_L=qlaplace(exp(-exp(-X)))

maxs=apply(X_L,1,max)
max.exceed.inds=which(maxs>u)
max.nonexceed.inds=which(maxs<=u)

n=50000

boo<-rbinom(n=n,size=1,prob=length(max.exceed.inds)/nrow(X_L)  ) 

sims.exceeds<-rsim_all(pars,n=sum(boo),B=5*sum(boo),v=u)

sims=matrix(nrow=n,ncol=3)

sims[1:sum(boo),]=sims.exceeds
orig.inds=sample(max.nonexceed.inds,size=sum(boo==0),replace=T)
sims[-(1:sum(boo)),]=X_L[orig.inds,]

sims=-log(-log(plaplace(sims)))


plot(X[,1],X[,2])
abline(v=-log(-log(0.9)),col="red")
abline(h=-log(-log(0.9)),col="red")
points(sims[,1],sims[,2],col="blue",pch=20)

plot(X[,2],X[,3])
abline(v=-log(-log(0.9)),col="red")
abline(h=-log(-log(0.9)),col="red")
points(sims[,2],sims[,3],col="blue",pch=20)


plot(X[,1],X[,3])
abline(v=-log(-log(0.9)),col="red")
abline(h=-log(-log(0.9)),col="red")
points(sims[,1],sims[,3],col="blue",pch=20)

