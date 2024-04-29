setwd(this.path::here())
data<-read.csv("../Data/Amaurot.csv")

library(evd)


Y=data$Y


library(evd)


ncores <- 40 # For parallel computation
K <- 350000 # Size of training data
n.samples <- 4e7 # n* number of observations to caluclate "theoretical" extreme quantile
n.boot <- 750 # Number of bootstrap samples

sigs<-xis<-us<-matrix(ncol=n.boot,nrow=21000)
for(i in 1:n.boot){
  try(load(paste0("../subchallenge-c1/GPD/0.6_GPD_predictions_test",i,".Rdata")),silent=T)
  
  us[,i]=u
  xis[,i]=preds[,2]
  sigs[,i]=preds[,1]
  
}

# Remove some problematic values
us[us>100]=NA
sigs[sigs>200]=NA


library(parallel)

print(paste("ncores=",ncores))

cl <- makeCluster(ncores,port=11300)
setDefaultCluster(cl=cl)
invisible(clusterEvalQ (cl , library("collapse")))
clusterExport(cl , "Y")


simulate <- function(xis,us,sigs, K,n.samples, seeds) {
  

  
  parLapply(seeds,cl=cl, fun=function(i) {
    set.seed(i)
    us=us[!is.na(us)]
    u.vec=sample(us[us > min(Y)],size=21000,replace=T)
    xi.vec=sample(xis[!is.na(xis)],size=21000,replace=T)
    
    sig.vec=sample(sigs[!is.na(sigs)],size=21000,replace=T)
    
    u.vec=sample(u.vec,size=n.samples,replace=T)
    xi.vec=sample(xi.vec,size=n.samples,replace=T)
    sig.vec=sample(sig.vec,size=n.samples,replace=T)
    
    Z<-rep(0, n.samples)
    boo=sample(c(0,1),prob=c(0.6,0.4),size=n.samples,replace=T)
    unique.u=unique(u.vec[boo==0])
    
    tab<-as.numeric(qtab(u.vec[boo==0]))
    
    Z[boo==0]=  unlist(apply(cbind(sort(unique.u),tab),1,function(x) sample(Y[Y < x[1]],size=x[2],replace=T)))
    
    Z[boo==1] <- ((1-runif(sum(boo==1),0,1))^(-xi.vec[boo==1])-1)*sig.vec[boo==1]/xi.vec[boo==1]+u.vec[boo==1]
    
    
    q=quantile(Z,prob=1-1/(300*200))
    
    Z<-c(sample(Z,21000,replace=T))
    
    dim(Z)=c(21000,1)
    #dim(Z)=c(1,21000)
    
    list(Z,q)
  })
  
  
}

clusterExport(cl , "xis")
clusterExport(cl , "us")
clusterExport(cl , "sigs")
set.seed(1)

print(n.samples)
clusterExport(cl , "n.samples")
invisible(clusterEvalQ (cl , set.seed(1)))

seeds=1:1000
clusterExport(cl , "seeds")

tmp <- simulate(xis,us,sigs,K=1000,n.samples,seeds)
Z_test<-lapply(tmp,"[[", 1)
theta_test<-t(as.matrix(unlist(lapply(tmp,"[[", 2))))
print("Test data simulated")

save(Z_test,theta_test,
     file=paste0("NBE/replicates/train_N",K,".Rdata"))
invisible(clusterEvalQ (cl , set.seed(2)))

seeds=1001:(1000+K/5)
clusterExport(cl , "seeds")

tmp <- simulate(xis,us,sigs,K=K/5,n.samples,seeds)
Z_val<-lapply(tmp,"[[", 1)
theta_val<-t(as.matrix(unlist(lapply(tmp,"[[", 2))))
print("Validation data simulated")

invisible(clusterEvalQ (cl , set.seed(3)))

seeds=(1001+K/5):(1000+K/5+K)
clusterExport(cl , "seeds")
tmp <- simulate(xis,us,sigs,K=K,n.samples, seeds)
Z_train<-lapply(tmp,"[[", 1)
theta_train<-t(as.matrix(unlist(lapply(tmp,"[[", 2))))

print("Training data simulated")
save(Z_train,Z_val,Z_test,theta_train,theta_test,theta_val,
     file=paste0("NBE/replicates/train_N",K,".Rdata"))


stopCluster(cl)

