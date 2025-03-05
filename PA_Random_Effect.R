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
set.seed(42)

u <- function(x,q) 5*x+q-2
logit <- function(x,q) exp(u(x,q))/(1+exp(u(x,q)))
Pa <- function(x,q) pnorm(u(x,q))
P10 <- function(x,q) Pa(x,q)*logit(x,q)
P11 <- function(x,q) logit(x,q)

# observe k households for T months
T = 24
k = 300
x <- runif(k)
q <- rnorm(k)
P <- P10(x,q)
D0 <- rbern(k,P)

Takeup_generator <- function(D0){
  prob_takeup <- ifelse(D0==1,P11(x,q),P10(x,q))
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
plot(q_SA,ylim = c(0,1000),xlab = "q",main = "",col = alpha("skyblue", 0.75))
plot(q_FA,ylim = c(0,1000),xlab = "q",main = "",col = alpha("red", 0.5),add = TRUE)

data_GQ <- data.frame(matrix(0,nrow = k*(T+1), ncol = 5))
colnames(data_GQ) <- c("id","month","WIC_1","X","q")
data_GQ$id <- rep(1:k,each = T+1)
data_GQ$month <- rep(seq(1,T+1),k)
data_GQ$WIC_1 <- c(Dall)
data_GQ$X <- rep(x,each = T+1)
data_GQ$q <- rep(q,each = T+1)

data_GQ <- data_GQ %>%
  group_by(id) %>%  # Group by household
  mutate(
    WIC_1_lag = lag(WIC_1, default = 0),  # Create lagged WIC_1 (default to 0 for the first month)
    ID0 = ifelse(month == 1 | WIC_1_lag == 0, 1, 0)  # Define ID0
  ) %>%
  ungroup()


utility <- function(q,a,b,x) a*x + q + b 
logit <- function(q,a,b,x){
  u <- utility(q,a,b,x)
  exp(u)/(1+exp(u))
} 

Pd <- (a,b,X){
  
}



nll <- function(a,b,c,d){
  
}

