rm(list=ls())

##### Loading needed packages
require(this.path)
require(evd)
require(heatmaply)
require(pheatmap)
require(ggplot2)
require(maxLik)

##################################
##### User-defined functions #####
##################################

##### Vector norm
mynorm <- function(x){
  sqrt(sum(x^2))
}

##### Log-likelihood
posloglik <- function(theta,x,thres = 0.2){
  sum(log(theta[1] + (theta[2]*(1 - theta[1]*thres)*x^(theta[2]-1))/(thres^theta[2])),na.rm = TRUE)
}

##### Nonparametric tail probability estimation
tail.prob <- function(u, q, t = rep(1,ncol(u)),thres = 0.2,init.par.up = c(20,5),num.rep = 10){
  
  N <- nrow(u)
  
  u2 <- array(NA, dim = dim(u))
  for(i in 1:ncol(u)){
    u2[,i] <- (1-u[,i])/t[i]
  }
  u_max <- apply(u2,MARGIN = 1,FUN = max,na.rm = TRUE)
  u_new <- u_max[u_max < thres]
  n_T <- length(u_new)
  
  theta.comp <- data.frame(k1 = rep(NA,num.rep),eta = rep(NA,num.rep),pll = rep(NA,num.rep))
  
  for(b in 1:num.rep){
    
    init.k1 <- runif(1,0,init.par.up[1])
    init.eta <- runif(1,1,init.par.up[2])
    A <- diag(c(1,1))
    B <- matrix(c(0,-1),nrow = 2,ncol = 1)
    theta.b <- maxLik(posloglik,start = c(init.k1,init.eta),x = u_new,thres = thres,
                      constraints = list(ineqA = A,ineqB = B))
    theta.comp$k1[b] <- theta.b$estimate[1]
    theta.comp$eta[b] <- theta.b$estimate[2]
    theta.comp$pll[b] <- theta.b$maximum
  }
  
  theta.est <- theta.comp[which.max(theta.comp$pll),1:2]
  
  prob.est <- n_T*(theta.est[1]*q + ((1 - theta.est[1]*thres)*q^theta.est[2])/(thres^theta.est[2]))/N
  
  return(as.numeric(prob.est))
}

##############################################################
##### Separating variables into tail-dependent subgroups #####
##############################################################

##### Importing data
U1 = read.csv("../Data/UtopulaU1.csv")
U2 = read.csv("../Data/UtopulaU2.csv")
U = cbind(U1, U2)

##### Calculating TPDM
n = dim(U)[1]; p = dim(U)[2]
U.trans = U
for(i in 1:p){
  u = rank(U[,i], ties.method = "random")/(n+1)
  U.trans[,i] = sqrt(-1/log(u))
}
TPDM <- matrix(nrow = p, ncol = p)
for(i in 1:p){
  for(j in i:p){
    print(c(i,j))
    R = sapply(1:n, function(t) mynorm(c(U.trans[t,i],U.trans[t,j])))
    r = quantile(R,0.99)
    idx.pair = which(R>r)
    sum_one = U.trans[idx.pair,i]*U.trans[idx.pair,j]/R[idx.pair]^2
    TPDM[i,j] <- TPDM[j,i] <- sum(sum_one)/length(idx.pair)
  }
}
# saveRDS(TPDM, "../subchallenge-c4/results/TPDM.rds")
# TPDM <- readRDS("../subchallenge-c4/results/TPDM.rds")

##### Selecting cut-off for defining subgroups
chi.mat <- TPDM
plot(sort(chi.mat)); abline(h = 0.12, col = "red")
## The selected cut-off is 0.12.

##### Visualizing tail dependent subgroups based on TPDM
heatmaply(TPDM,ylab="",xlab="",k_row=5,k_col=5,
          dendrogram = T, showticklabels = F,show_dendrogram = c(TRUE,FALSE),
          colors = (blues9),key.title="EDM", limits=c(0,0.5),
          height=6,width=8)

##### Separating variables into five subgroups
pairs = which(chi.mat>=0.12, arr.ind = T)
pairs = data.frame(pairs)

idx.bi = c()
pairs.bi = list()
for (i in 1:(nrow(pairs)-1)){
  for (j in (i+1):nrow(pairs)) {
    if (pairs[i,1] == pairs[j,2] & pairs[i,2] == pairs[j,1]) {
      pairs.bi = c(pairs.bi, list(pairs[i,]))
      idx.bi = c(idx.bi, i, j)
    }}}
pairs.bi = data.frame(do.call("rbind", pairs.bi))

isolated.points = (1:50)[-unique(as.numeric(do.call("rbind", pairs)))]
if (length(isolated.points) != 0) {
  groups0 = c(lapply(1:length(isolated.points), function(i) isolated.points[i]), 
              list(as.numeric(pairs[1,])))
} else {
  groups0 = list(unique(as.numeric(pairs[1,])))
}
for (i in 2:nrow(pairs)) {
  print(i)
  print(length(groups0))
  group.logic = 0
  for (j in 1:length(groups0)) {
    pair.logic = as.numeric(pairs[i,]) %in% groups0[[j]]
    group.logic = group.logic + sum(pair.logic)
    if (pair.logic[1] == F & pair.logic[2] == T) {
      groups0[[j]] = c(groups0[[j]], pairs[i,1])
    } else if (pair.logic[1] == T & pair.logic[2] == F) {
      groups0[[j]] = c(groups0[[j]], pairs[i,2])
    }
  }
  if (group.logic == 0) {
    groups0 = c(groups0, list(unique(as.numeric(pairs[i,]))))
  }
}

groups = groups0
while (TRUE) {
  groups.old = groups
  gpair.include = sapply(1:length(groups.old), function(i) {
    sapply(1:length(groups.old), function(j) {
      ifelse(sum(groups.old[[i]] %in% groups.old[[j]])>0, T, F)
    }) })
  if (sum(gpair.include) - length(groups.old) == 0) {break}
  
  groups = list()
  keep.merge = TRUE
  merge.id = 1; merged.gp = rep(FALSE, ncol(gpair.include))
  while (keep.merge) {
    print(merge.id)
    if (!merged.gp[merge.id]) {
      newgp = sort(unique(unlist(groups.old[which(gpair.include[merge.id,])])))
      groups = c(groups, list(newgp))
      merged.gp[gpair.include[merge.id,]] = TRUE
      # merged.gp = merged.gp + sum(gpair.include[merge.id,])
      # print(merged.gp)
      if (sum(merged.gp) == length(groups.old)) { keep.merge = FALSE }
    } 
    merge.id = merge.id + 1
  }
}
groups # obtaining variables belonging in each subgroup


##########################################################
##### Nonparametric estimation of tail probabilities #####
##########################################################

##### Setting function parameters
q <- 1/300 # quantile level
t1 <- c(rep(1,25),rep(12,25)) # constant multiplier for each variable for p3
t2 <- c(rep(1,50)) # constant multiplier for each variable for p3

##### Transformed data based on the Gumbel margins
U.unif <- exp(-exp(-U))

##### Selecting threshold for estimating tail probabilities
ps1_comp <- c()
ps2_comp <- c()

set.seed(1) # seed number for reproducibility

thres.all <- seq(0.1,0.4,by = 0.005) # grid of thresholds

for(ind in 1:length(thres.all)){

  ## Variable groupings
  g.num <- length(groups)

  ## Probabilities per cluster
  p.est1 <- vector("numeric",length(groups))
  p.est2 <- vector("numeric",length(groups))


  for(g in 1:g.num){

    var.list <- groups[[g]]
    t1.list <- t1[var.list]
    t2.list <- t2[var.list]

    p.est1[g] <- tail.prob(u = U.unif[,var.list],q = q,t = t1.list,thres = thres.all[ind],init.par.up = c(2.5,1.5),num.rep = 200)
    p.est2[g] <- tail.prob(u = U.unif[,var.list],q = q,t = t2.list,thres = thres.all[ind],init.par.up = c(2.5,1.5),num.rep = 200)

    cat(paste("Threshold = ",thres.all[ind],", group = ",g,sep = ""),"\r")

  }

  ps1_comp <- c(ps1_comp,cumprod(p.est1)[length(p.est1)])
  ps2_comp <- c(ps2_comp,cumprod(p.est2)[length(p.est2)])

}

# saveRDS(ps1_comp,"../subchallenge-c4/results/ps1_comp.rds")
# saveRDS(ps2_comp,"../subchallenge-c4/results/ps2_comp.rds")

##### Plotting estimated tail probabilities for each threshold
# ps1_comp <- readRDS("../subchallenge-c4/results/ps1_comp.rds")
# ps2_comp <- readRDS("../subchallenge-c4/results/ps2_comp.rds")

### Estimated tail probabilities per threshold for p3
p3 <- ggplot(data=data.frame(x=thres.all,
                         out=log(ps2_comp)), aes(x=thres.all, y=log(ps2_comp), group=1)) +
  geom_point()+theme(
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20))+
  geom_vline(xintercept=0.25,col="red")+
  geom_point()+xlab( expression(phi["*"]))+
  ylab( expression(log(hat(p[3]))))
p3
ggsave(p3, filename = "../Figures/p3.pdf", height = 3.5, width = 7, bg = "transparent" )

### Estimated tail probabilities per threshold for p4
p4 <- ggplot(data=data.frame(x=thres.all,
                            out=log(ps1_comp)), aes(x=thres.all, y=log(ps1_comp), group=1)) +
  geom_point()+theme(
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20))+
  geom_vline(xintercept=0.25,col="red")+
  geom_point()+xlab( expression(phi["*"]))+
  ylab( expression(log(hat(p[4]))))
p4
ggsave(p4, filename = "../Figures/p4.pdf", height = 3.5, width = 7, bg = "transparent" )

## The selected threshold is 0.25.

##### Estimates for threshold = 0.25
ps1_comp[31] # estimate for p3
ps2_comp[31] # estimate for p4

##### Bootstrap for thres = 0.25
ps1_bs <- c()
ps2_bs <- c()

set.seed(1) # seed number for reproducibility
B <- 1000 # number of bootstrap replicates
for(bs in 1:B){

  samp_ind <- sample(1:nrow(U.unif),nrow(U.unif),replace = TRUE)

  U.unif_bs <- U.unif[samp_ind,]

  ## Variable groupings
  g.num <- length(groups)
  
  ## Probabilities per cluster
  p.est1 <- vector("numeric",length(groups))
  p.est2 <- vector("numeric",length(groups))

  for(g in 1:g.num){

    var.list <- groups[[g]]
    t1.list <- t1[var.list]
    t2.list <- t2[var.list]

    p.est1[g] <- tail.prob(u = U.unif_bs[,var.list],q = q,t = t1.list,thres = 0.25,init.par.up = c(2.5,1.5),num.rep = 200)
    p.est2[g] <- tail.prob(u = U.unif_bs[,var.list],q = q,t = t2.list,thres = 0.25,init.par.up = c(2.5,1.5),num.rep = 200)

    cat(paste("BS number = ",bs,", group = ",g,sep = ""),"\r")

  }

  ps1_bs <- c(ps1_bs,cumprod(p.est1)[length(p.est1)])
  ps2_bs <- c(ps2_bs,cumprod(p.est2)[length(p.est2)])

}

# saveRDS(ps1_bs,"../subchallenge-c4/results/ps1_bs.rds")
# saveRDS(ps2_bs,"../subchallenge-c4/results/ps2_bs.rds")

##### Empirical distribution of estimates at threshold = 0.25
# ps1_bs <- readRDS("../subchallenge-c4/results/ps1_bs.rds")
# ps2_bs <- readRDS("../subchallenge-c4/results/ps2_bs.rds")
par(mfrow = c(1,2))
boxplot(log(ps1_bs,base = 10))
boxplot(log(ps2_bs,base = 10))

#### Final answer for C4 is the median of the bootstrap estimates
c4_answer <- c(median(ps1_bs),median(ps2_bs))
# saveRDS(c4_answer,"../subchallenge-c4/results/C4_Answer.rds")
# c4_answer <- readRDS("../subchallenge-c4/results/C4_Answer.rds")
c4_answer

