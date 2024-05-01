# Clear all previous R environment data and load necessary libraries
rm(list=ls())
require(evgam)
require(qgam)
require(ggeffects) 
require(gratia)
require(parallel)
require(psych)
require(mgcv)
require(this.path)
setwd(this.path::here())


# Load and inspect the test dataset containing different covariate combinations
data.test <- read.csv("../Data/AmaurotTestSet.csv")
colnames(data.test)

# Load and inspect the main dataset for Amaurot
data <- read.csv("../Data/Amaurot.csv")
colnames(data)
class(data)

# Convert 'Season' from numerical or character to a factor (categorical variable)
data$Season <- as.factor(data$Season)

# Prepare the dataset by combining test and training data, handling missing values
n0 <- dim(data)[1]
n1 <- dim(data.test)[1]
data.testY <- cbind(rep(NA, n1), data.test)
colnames(data.testY) <- colnames(data)

# Combine test and training data
data2 <- rbind(data.testY, data)

# Handle missing data by imputing mean values and creating indicators for missingness
dname <- names(data)[-1]
dat1 <- data2
n = nrow(data2)
for (i in 1:8) {
  by.name <- paste("m", dname[i], sep = "")
  dat1[[by.name]] <- is.na(dat1[[dname[i]]])
  if (sum(is.na(dat1[[dname[i]]])) > 0) {
    dat1[[dname[i]]][dat1[[by.name]]] <- mean(dat1[[dname[i]]], na.rm = TRUE)
    lev <- rep(1, n)
    lev[dat1[[by.name]]] <- 1:sum(dat1[[by.name]])
    id.name <- paste("id", dname[i], sep = "")
    dat1[[id.name]] <- factor(lev)
    dat1[[by.name]] <- as.numeric(dat1[[by.name]])
  } else {
    print(i)
  }
}

# Transform variable V1 using natural logarithm for better modeling
dat2 <- dat1
dat2$V1 <- log(dat1$V1)

# Define the quantile to be estimated and fit a quantile regression model using qgam package

# Set lower q0 and n.boots values for use in subchallenge C2
q0 = 0.972 # q0 <- 0.6
n.boots <- 2500 # n.boots <- 750




fit <- qgam(Y ~ 1 + Season + WindSpeed + s(V2, by = ordered(!mV2)) +
              s(V3, by = ordered(!mV3)) + s(V4, by = ordered(!mV4)) +
              s(WindDirection, by = ordered(!mWindDirection)),
            data = dat2[-(1:100),], qu = q0)

# Predict exceedance thresholds
u = predict(fit, newdata = dat2[-(1:100),])

# Calculate residuals where Y exceeds the threshold u
excess <- dat2$Y[-(1:100)] - u
dat3 <- cbind(excess, dat2[-(1:100),])
is.na(dat3$excess[dat3$excess <= 0]) <- TRUE

# Fit generalized Pareto distribution models to model exceedances
fmla_gpd1 <- list(excess ~ Season + WindSpeed + s(V2, by = ordered(!mV2)) + s(V1, by = ordered(!mV1)) +
                    s(V3, by = ordered(!mV3)) + s(V4, by = ordered(!mV4)) +
                    s(WindDirection, by = ordered(!mWindDirection)),
                  ~ Season + WindSpeed)
m_gpd1 <- evgam(fmla_gpd1, dat3, family = "gpd")
summary(m_gpd1)

fmla_gpd2 <- list(excess ~ WindSpeed + s(V2, by = ordered(!mV2)) +
                    s(V3, by = ordered(!mV3)) + s(V4, by = ordered(!mV4)) +
                    s(WindDirection, by = ordered(!mWindDirection)),
                  ~ Season + WindSpeed)
m_gpd2 <- evgam(fmla_gpd2, dat3, family = "gpd")
summary(m_gpd2)

# Compare model performance using Bayesian Information Criterion (BIC)
BIC(m_gpd1)
BIC(m_gpd2)

# Bootstrap the quantile estimation to provide confidence intervals
quants.9999 <- mclapply(1:n.boots, FUN = function(i) {
  dat.2 <- dat2[-(1:100),]
  it = i + 526
  set.seed(it)
  boot.inds = sample(1:n0, size = n0, replace = T)
  
  fit <- qgam(Y ~ 1 + WindSpeed + s(V2, by = ordered(!mV2)) +
                s(V3, by = ordered(!mV3)) + s(V4, by = ordered(!mV4)) +
                s(WindDirection, by = ordered(!mWindDirection)),
              data = dat.2[boot.inds,], qu = q0)
  
  # Predict exceedance thresholds for the bootstrap sample and test data
  u.boot = predict(fit, newdata = dat.2[boot.inds,])
  u = predict(fit, newdata = dat.2)
  u.test = predict(fit, newdata = dat2[1:100,])
  
  excess <- dat.2$Y[boot.inds] - u.boot
  dat3.boot <- cbind(excess, dat.2[boot.inds,])
  is.na(dat3.boot$excess[dat3.boot$excess <= 0]) <- TRUE
  
  # Fit the GPD model to the bootstrap sample
  fmla_gpd1 <- list(excess ~ WindSpeed + s(V2, by = ordered(!mV2)) +
                      s(V3, by = ordered(!mV3)) + s(V4, by = ordered(!mV4)) +
                      s(WindDirection, by = ordered(!mWindDirection)),
                    ~ Season + WindSpeed)
  m_gpd_1 <- evgam(fmla_gpd1, dat3.boot, family = "gpd")
  
  summary(m_gpd_1)
  
  # Predict quantiles for the test data
  quants <- predict(m_gpd_1, newdata = dat2[(1:100),], prob = (0.9999 - q0) / (1 - q0)) + u.test
  
  # If q0 = 0.6, save predicted GPD fits for use in subchallenge C2
  if( q0 == 0.6){
   preds<-predict(m_gpd_1, newdata=dat3, type="response")

   save(preds, u,q0, file=paste0("GPD/",q0,"_GPD_predictions_test",i,".Rdata"))
  
  }
  return(quants)
}, mc.cores = 40)

# Combine quantile estimates from all bootstraps
q.9999 <- predict(m_gpd2, newdata = dat2[(1:100),], prob = (0.9999 - q0) / (1 - q0)) + predict(fit, newdata = dat2[(1:100),])
quantile.9999 <- matrix(unlist(quants.9999), ncol = n.boots, nrow = 100)

# Calculate confidence intervals for the quantile estimates
CI.q.9999 <- apply(quantile.9999, FUN = quantile, probs = c(.5, .25, .75), MARGIN = 1)

# Prepare the final output with quantile estimates and confidence intervals
AnswerC1 = cbind(q.9999, t(CI.q.9999))

# Save the results in a file named 'C1.Rdata'
save(AnswerC1, file = "C1.Rdata")
