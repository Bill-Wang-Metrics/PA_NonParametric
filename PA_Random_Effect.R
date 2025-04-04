library(Rlab)
library(gplm)
library(stats4)
library(plotrix)
library(MASS)
library(grf)
library(scales)
library(dplyr)
library(EstimationTools)
library(gaussquad)
library(statmod)
set.seed(42)

u2 <- function(x1,x2,q) 5*x1 - 3*x2 + q - 1
u1 <- function(x1,x2,q) 3*x1 - 5*x2 + q - 2

logit <- function(x1,x2,q) exp(u2(x1,x2,q))/(1+exp(u2(x1,x2,q)))
Pa <- function(x1,x2,q) exp(u1(x1,x2,q))/(1+exp(u1(x1,x2,q)))
P10 <- function(x1,x2,q) Pa(x1,x2,q)*logit(x1,x2,q)
P11 <- function(x1,x2,q) logit(x1,x2,q)

# observe k households for T months
T = 24
k = 300
x1 <- runif(k)
x2 <- rbern(k,.5) - 0.5
q <- rnorm(k)
P <- P10(x1,x2,q)
D0 <- rbern(k,P)

Takeup_generator <- function(D0){
  prob_takeup <- ifelse(D0==1,P11(x1,x2,q),P10(x1,x2,q))
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


sum(Dall[-1,])/k/T 


SA <- matrix(nrow = 0, ncol = 3)
FA <- matrix(nrow = 0, ncol = 3)

# Iterate over time periods and households
for (t in 2:(T + 1)) {  # Start from t = 2 because t = 1 is the initial period
  for (i in 1:k) {
    # Check participation status in the previous period
    if (Dall[t - 1, i] == 0) {
      # If not in the program last period, add to SA
      SA <- rbind(SA, c(Dall[t, i], x[i], q[i]))
    } else if (Dall[t - 1, i] == 1) {
      # If in the program last period, add to FA
      FA <- rbind(FA, c(Dall[t, i], x[i], q[i]))
    }
  }
}

# Add column names for clarity
colnames(SA) <- c("Participation", "x","q")
colnames(FA) <- c("Participation", "x","q")

par(mfrow = c(1,1))
x_FA <- hist(FA[,2])
x_SA <- hist(SA[,2])
plot(x_SA,ylim = c(0,1000),xlab = "x",main = "",col = alpha("skyblue", 0.75))
plot(x_FA,ylim = c(0,1000),xlab = "x",main = "",col = alpha("red", 0.5),add = TRUE)

q_FA <- hist(FA[,3])
q_SA <- hist(SA[,3])
plot(q_SA,ylim = c(0,1500),xlab = "q",main = "",col = alpha("skyblue", 0.75))
plot(q_FA,ylim = c(0,1500),xlab = "q",main = "",col = alpha("red", 0.5),add = TRUE)

data_GQ <- data.frame(matrix(0,nrow = k*(T+1), ncol = 6))
colnames(data_GQ) <- c("id","month","WIC_1","X1","X2","q")
data_GQ$id <- rep(1:k,each = T+1)
data_GQ$month <- rep(seq(1,T+1),k)
data_GQ$WIC_1 <- c(Dall)
data_GQ$X1 <- rep(x1,each = T+1)
data_GQ$X2 <- rep(x2,each = T+1)
data_GQ$q <- rep(q,each = T+1)

data_GQ <- data_GQ %>%
  group_by(id) %>%  # Group by household
  mutate(
    WIC_1_lag = lag(WIC_1, default = 0),  # Create lagged WIC_1 (default to 0 for the first month)
    ID0 = ifelse(month == 1 | WIC_1_lag == 0, 1, 0)  # Define ID0
  ) %>%
  ungroup()


utility <- function(Q,a,b,c,sigma,x){
  n = dim(x)[1]
  x1 <- x$X1
  x2 <- x$X2
  q <- rep(Q, length.out = n)
  r <- sqrt(2*sigma^2) * q
  a * x1 + b * x2 + r + rep(c, length.out = n) 
}

ua <- function(Q,d,e,f,sigma,x){
  n = dim(x)[1]
  x1 <- x$X1
  x2 <- x$X2
  q <- rep(Q, length.out = n)
  r <- sqrt(2*sigma^2) * q
  d * x1 + e * x2 + r + rep(f, length.out = n) 
}

Palogit <- function(Q,a,b,c,d,e,f,sigma,data_GQid){
  u2 <- utility(Q,a,b,c,sigma,data_GQid)
  u1 <- ua(Q,d,e,f,sigma,data_GQid)
  Pit <- (exp(u1)/(1+exp(u1)))^(data_GQid$ID0 == 1) * (exp(u2)/(1+exp(u2)))
  # Pit <- exp(u)/(1+exp(u))
  prod(Pit^(data_GQid$WIC_1 == 1)*(1 - Pit)^(data_GQid$WIC_1 == 0) + 1e-10) 
} 

quadrature <- gauss.quad(n = 100, kind = "hermite")
nodes <- quadrature$nodes
weights <- quadrature$weights



nll <- function(params){
  a <- params[1]
  b <- params[2]
  c <- params[3]
  d <- params[4]
  e <- params[5]
  f <- params[6]
  sigma <- exp(params[7])
  
  neg_ll <- 0
  for(i in 1:k){
    data_GQid <- data_GQ[data_GQ$id == i,]
    
    # Evaluate the integrand at the quadrature nodes
    integrand_values <- sapply(nodes, function(Q) {
      Palogit(Q, a, b, c, d, e, f, sigma, data_GQid)
    })
    
    # Approximate the integral using Gaussian quadrature
    integral <- sum(weights * integrand_values)
    
    neg_ll <- neg_ll - log(integral)
  }
  return(neg_ll)
}


optim(par = c(2,2,2,2,2,2,2), fn = nll, method = "BFGS")