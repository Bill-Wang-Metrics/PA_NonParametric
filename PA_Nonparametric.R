library(Rlab)
library(gplm)
library(stats4)
library(plotrix)
library(MASS)
library(grf)
set.seed(42)

logit <- function(x) exp(2*x)/(1+exp(2*x))
#Pa <- function(x) 0.7 - ((0.5*(sin(pi*(2*x-1)))+.5)^2*logit((2*x-1))-0.1*(5*x-2))
#Pa <- function(x) logit(x)
#Pa <- function(x) rep(0.3,length(x))
Pa <- function(x) 0.5 - (.5*(cos(pi*(2*x-1)))+.5)*logit((2*x-1))^2 + (x-0.5)^2
plotTRUE <- function(){
  x <- seq(0,1,0.01)
  par(mar=c(5,5,2,1)+.1)
  plot(x,Pa(x),type = "l",lwd = 3, ylab = "Attention Probability",xlab = "X",cex.lab=2,
  xlim = c(0,1), ylim = c(0,1),yaxt="n",xaxt="n")
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
}
par(mfrow = c(1,1))
plotTRUE()


P10 <- function(x) Pa(x)*logit(x)
P11 <- function(x) logit(x)

# observe k households for T months
T = 24
k = 300
x <- runif(k)
P <- P10(x)
D0 <- rbern(k,P)

Takeup_generator <- function(D0){
  prob_takeup <- ifelse(D0==1,P11(x),P10(x))
  Dt <- rbern(k,prob_takeup)
}

Dall <- t(matrix(D0))
for (t in 2:(T+1)){
  D0 <- Takeup_generator(Dall[t-1,])
  Dall <- rbind(Dall,D0)
}


colnames(Dall) <- NULL
rownames(Dall) <- NULL
#Dall <- Dall[-1,]
#View(Dall)


sum(Dall[-1,])/k/T # 0.5304167 when k = 100, T =24



# Initialize empty matrices for SA and FA
SA <- matrix(nrow = 0, ncol = 2)
FA <- matrix(nrow = 0, ncol = 2)

# Iterate over time periods and households
for (t in 2:(T + 1)) {  # Start from t = 2 because t = 1 is the initial period
  for (i in 1:k) {
    # Check participation status in the previous period
    if (Dall[t - 1, i] == 0) {
      # If not in the program last period, add to SA
      SA <- rbind(SA, c(Dall[t, i], x[i]))
    } else if (Dall[t - 1, i] == 1) {
      # If in the program last period, add to FA
      FA <- rbind(FA, c(Dall[t, i], x[i]))
    }
  }
}

# Add column names for clarity
colnames(SA) <- c("Participation", "x")
colnames(FA) <- c("Participation", "x")
x = seq(0,1,0.05)

par(mfrow = c(2,2))
#par(mfrow = c(1,1))

######## constant PA
nll <- function(a,c,d){
  PA <- a
  PDSA <- exp(c*SA[,2]+d)/(1+exp(c*SA[,2]+d))
  PDFA <- exp(c*FA[,2]+d)/(1+exp(c*FA[,2]+d))
  -sum(SA[,1]*log(PA*PDSA) - PA*PDSA) - sum(FA[,1]*log(PDFA) - PDFA)
}

est <- stats4::mle(minuslog=nll, start=list(a=2,c=2,d=0))
plotTRUE()
#x = seq(0.05,0.95,0.1)
PA <- rep(coef(est)[1],length(x))
interval <- rep(1.96 * sqrt(vcov(est)[1,1]),length(x))
par(new = TRUE)
plotCI(x,PA,
       li = ifelse(PA - interval>0,PA - interval,0),
       ui = ifelse(PA + interval<1,PA + interval,1),  
       ylab = "",xlab = "",yaxt="n",xaxt="n",
       xlim = c(0,1), ylim = c(0,1),col = "darkgreen", pch=19)

######## logit PA
nll <- function(a,b,c,d){
  PA <- exp(a*SA[,2]+b)/(1+exp(a*SA[,2]+b))
  PDSA <- exp(c*SA[,2]+d)/(1+exp(c*SA[,2]+d))
  PDFA <- exp(c*FA[,2]+d)/(1+exp(c*FA[,2]+d))
  -sum(SA[,1]*log(PA*PDSA) - PA*PDSA) - sum(FA[,1]*log(PDFA) - PDFA)
}

est <- stats4::mle(minuslog=nll, start=list(a=2,b=0,c=2,d=0))
#x = seq(0.05,0.95,0.1)
u <- coef(est)[1]*x + coef(est)[2]
PA <- exp(u)/(1+exp(u))

plotTRUE()
par(new = TRUE)
first <- cbind(x*PA^2*exp(-u),exp(-u)*PA^2)
second <- vcov(est)[1:2,1:2]
interval <- 1.96 * sqrt(rowSums(first %*% second * first))
plotCI(x,PA,
       li = ifelse(PA - interval>0,PA - interval,0),
       ui = ifelse(PA + interval<1,PA + interval,1), 
       ylab = "",xlab = "",yaxt="n",xaxt="n",
       xlim = c(0,1), ylim = c(0,1),col = "blue", pch=19)


######## PA SIEVE
logit_FA <- glm(Participation ~ x, family = binomial(link = "logit"), data = as.data.frame(FA))
basis_m <- floor(log(nrow(SA)))
SIEVE_SA <- lm(Participation ~ poly(x,basis_m, raw = TRUE), data = as.data.frame(SA))
newdata <- as.data.frame(cbind(rep(1,length(x)),x))
SIEVE_SAx <- predict(SIEVE_SA,newdata = newdata)
logit_SA <- predict(logit_FA,newdata = newdata,type = "response")

plotTRUE()
par(new = TRUE)
Z = cbind(1,poly(SA[,2],basis_m, raw = TRUE))
Zx = cbind(1,poly(x,basis_m, raw = TRUE))
ZZinv = ginv(t(Z)%*%Z) # same as solve(t(Z)%*%Z) when t(Z)%*%Z is invertible
u_hat <- SIEVE_SA$residuals
Zu <- sweep(t(Z),2,u_hat,"*")
Vn = Zx %*% ZZinv %*% Zu %*% t(Zu) %*% ZZinv %*% t(Zx)
An <- 1/sqrt(diag(Vn))
interval <- 1.96 * sqrt((1/logit_SA/An)^2) 
plotCI(x,SIEVE_SAx/logit_SA,
     li = ifelse(SIEVE_SAx/logit_SA - interval>0,SIEVE_SAx/logit_SA - interval,0),
     ui = ifelse(SIEVE_SAx/logit_SA + interval<1,SIEVE_SAx/logit_SA + interval,1), 
     ylab = "",xlab = "",yaxt="n",xaxt="n",
     xlim = c(0,1), ylim = c(0,1),col = "red", pch=19)


######## PA RF
RF_SA <- regression_forest(matrix(SA[,2]),matrix(SA[,1]),min.node.size = min(dim(SA)[1]/32,50),seed = 42)
RF_SAx <- predict(RF_SA, matrix(x),estimate.variance = TRUE)
plotTRUE()
par(new = TRUE)
interval = 1.96 * sqrt(RF_SAx$variance.estimates)/logit_SA
plotCI(x,RF_SAx$predictions/logit_SA,
       li = ifelse(RF_SAx$predictions/logit_SA - interval>0,RF_SAx$predictions/logit_SA - interval,0),
       ui = ifelse(RF_SAx$predictions/logit_SA + interval<1,RF_SAx$predictions/logit_SA + interval,1), 
       ylab = "",xlab = "",yaxt="n",xaxt="n",
       xlim = c(0,1), ylim = c(0,1),col = "darkorange", pch=19)

# Boost_SA <- boosted_regression_forest(matrix(SA[,2]),matrix(SA[,1]))
# Boost_SAx <- predict(RF_SA, matrix(x))
# plotTRUE()
# par(new = TRUE)
# plot(x,Boost_SAx$predictions/logit_SA,ylab = "Attention Probability",xlab = "X",
#      xlim = c(0,1), ylim = c(0,1),col = "darkorchid",pch=19)



