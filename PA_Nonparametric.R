library(Rlab)
library(gplm)
library(stats4)
set.seed(42)

logit <- function(x) exp(2*x)/(1+exp(2*x))
Pa <- function(x) (0.5*(sin(pi*(2*x-1)))+.5)^2*logit((2*x-1))-0.1*(2*x-2)
x <- seq(0,1,0.01)
plot(x,Pa(x),type = "l",lwd = 3, ylab = "Attention Probability",xlab = "X")


P10 <- function(x) Pa(x)*logit(x)
P11 <- function(x) logit(x)

# observe k households for T months
T = 24
k = 100
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


logit_FA <- glm(Participation ~ x, family = binomial(link = "logit"), data = as.data.frame(FA))

SIEVE_SA <- lm(Participation ~ poly(x,floor(log(nrow(SA)))), data = as.data.frame(SA))


logit_SA <- predict(logit_FA,newdata = as.data.frame(SA),type = "response")

plot(SA[,2],SIEVE_SA$fitted.values/logit_SA,
     xlim = c(0,1), ylim = c(0,1),ylab = "",xlab = "",col = "red")
par(new = TRUE)
x <- seq(0,1,0.01)
plot(x,Pa(x),type = "l",lwd = 3, ylab = "Attention Probability",xlab = "X",
     xlim = c(0,1), ylim = c(0,1))

# nll <- function(theta0,theta1) {
#   x <- Y$age[-idx]
#   y <- Y$Count[-idx]
#   mu = exp(theta0 + x*theta1)
#   -sum(y*(log(mu)) - mu)
# }

nll <- function(a,b,c,d){
  PA <- exp(a*SA[,2]+b)/(1+exp(a*SA[,2]+b))
  PDSA <- exp(c*SA[,2]+d)/(1+exp(c*SA[,2]+d))
  PDFA <- exp(c*FA[,2]+d)/(1+exp(c*FA[,2]+d))
  -sum(SA[,1]*log(PA*PDSA) - PA*PDSA) - sum(FA[,1]*log(PDFA) - PDFA)
}

est <- stats4::mle(minuslog=nll, start=list(a=2,b=0,c=2,d=0))

u <- coef(est)[1]*x + coef(est)[2]
PA <- exp(u)/(1+exp(u))
par(new = TRUE)
plot(x,PA, ylab = "Attention Probability",xlab = "X",
     xlim = c(0,1), ylim = c(0,1),col = "blue")