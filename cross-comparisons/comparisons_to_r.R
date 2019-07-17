##################################################################
# R-code for testing/comparisons I have made with zEpid
##################################################################
# Last edit; 2019/3/30

### TMLE testing ###
library(dplyr)
library(tmle)

data <- read.csv(file="sample.csv", header=TRUE, sep=",")

# Testing complete case
cc <- data[complete.cases(data$dead), ]
Y <- cc$dead
A <- cc$art
W <- select(cc, male, age0, age_rs1, age_rs2, cd40, cd4_rs1, cd4_rs2, dvl0)
result <- tmle(Y,A,W,
               gform = A ~ male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               Qform = Y ~ A + male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               gbound = 0,
               family = "binomial", 
               verbose = T)
summary(result)
result <- tmle(Y,A,W,
               gform = A ~ male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               Qform = Y ~ A + male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               gbound = c(0.025, 0.975),
               family = "binomial", 
               verbose = T)
summary(result)
result <- tmle(Y,A,W,
               gform = A ~ male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               Qform = Y ~ A + male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               gbound = c(0.025, 0.95),
               family = "binomial", 
               verbose = T)
summary(result)

Y <- cc$cd4_wk45
A <- cc$art
W <- select(cc, male, age0, age_rs1, age_rs2, cd40, cd4_rs1, cd4_rs2, dvl0)
result <- tmle(Y,A,W,
               gform = A ~ male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               Qform = Y ~ A + male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               gbound = 0,
               family = "gaussian", 
               verbose = T)
summary(result)
result <- tmle(Y,A,W,
               gform = A ~ male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               Qform = Y ~ A + male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               gbound = 0,
               family = "poisson", 
               verbose = T)
summary(result)

# Testing missing binary outcome data
Y <- data$dead
A <- data$art
W <- select(data, male, age0, age_rs1, age_rs2, cd40, cd4_rs1, cd4_rs2, dvl0)
D <- data$miss_d
result <- tmle(Y,A,W, Delta=D,
               gform = A ~ male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               Qform = Y ~ A + male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               g.Deltaform = Delta ~ A + male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               gbound = 0,
               family = "binomial", 
               verbose = T)
summary(result)
result$estimates$ATE
result$estimates$RR
result$estimates$OR

# Testing missing continuous outcome data
Y <- data$cd4_wk45
A <- data$art
W <- select(data, male, age0, age_rs1, age_rs2, cd40, cd4_rs1, cd4_rs2, dvl0)
D <- data$miss_c
result <- tmle(Y,A,W,Delta = D,
               gform = A ~ male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               Qform = Y ~ A + male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               g.Deltaform = Delta ~ A + male + age0 + age_rs1 + age_rs2 + cd40 + cd4_rs1 + cd4_rs2 + dvl0,
               gbound = 0,
               family="gaussian", verbose = T)
summary(result)
result$estimates$ATE

### Iterative Conditional g-formula testing ###
library(ltmle)

data <- read.csv(file="C:/Users/zivic/Desktop/data.csv", header=TRUE, sep=",")
ids <- data$id
drops <- c("W", "id")
data <- data[ , !(names(data) %in% drops)]
data$Y2[data$Y1 == 1] <- 1
data$Y3[data$Y2 == 1] <- 1

result <- ltmle(data, Anodes=c("A1", "A2", "A3"), Cnodes=NULL, Lnodes=c("L1", "L2", "L3"), Ynodes=c("Y1", "Y2", "Y3"), abar=c(1, 1, 1),
                Qform=c(Y1='Q.kplus1 ~ A1 + L1', Y2='Q.kplus1 ~ A2 + L2', Y3='Q.kplus1 ~ A3 + L3'), 
                gcomp=TRUE, SL.library=NULL, survivalOutcome = TRUE, id=ids)
result$estimates  # Treat-all = 0.4140978

result <- ltmle(data, Anodes=c("A1", "A2", "A3"), Cnodes=NULL, Lnodes=c("L1", "L2", "L3"), Ynodes=c("Y1", "Y2", "Y3"), abar=c(0, 0, 0),
                Qform=c(Y1='Q.kplus1 ~ A1 + L1', Y2='Q.kplus1 ~ A2 + L2', Y3='Q.kplus1 ~ A3 + L3'), 
                gcomp=TRUE, SL.library=NULL, survivalOutcome = TRUE, id=ids)
result$estimates  # Treat-all = 0.6464508

drops <- c("A3", "Y3", "L3")
data2 <- data[ , !(names(data) %in% drops)]
result <- ltmle(data2, Anodes=c("A1", "A2"), Cnodes=NULL, Lnodes=c("L1", "L2"), Ynodes=c("Y1", "Y2"), abar=c(1, 1),
                Qform=c(Y1='Q.kplus1 ~ A1 + L1', Y2='Q.kplus1 ~ A2 + L2'), 
                gcomp=TRUE, SL.library=NULL, survivalOutcome = TRUE, id=ids)
result$estimates  # Treat-all = 0.3349204

result <- ltmle(data2, Anodes=c("A1", "A2"), Cnodes=NULL, Lnodes=c("L1", "L2"), Ynodes=c("Y1", "Y2"), abar=c(0, 0),
                Qform=c(Y1='Q.kplus1 ~ A1 + L1', Y2='Q.kplus1 ~ A2 + L2'), 
                gcomp=TRUE, SL.library=NULL, survivalOutcome = TRUE, id=ids)
result$estimates  # Treat-all = 0.5122774

result <- ltmle(data, Anodes=c("A1", "A2", "A3"), Cnodes=NULL, Lnodes=c("L1", "L2", "L3"), Ynodes=c("Y1", "Y2", "Y3"), abar=c(1, 0, 1),
                Qform=c(Y1='Q.kplus1 ~ A1 + L1', Y2='Q.kplus1 ~ A2 + L2', Y3='Q.kplus1 ~ A3 + L3'), 
                gcomp=TRUE, SL.library=NULL, survivalOutcome = TRUE, id=ids)
result$estimates  # Treat-all = 0.4916937

result <- ltmle(data, Anodes=c("A1", "A2", "A3"), Cnodes=NULL, Lnodes=c("L1", "L2", "L3"), Ynodes=c("Y1", "Y2", "Y3"), abar=c(0, 1, 0),
                Qform=c(Y1='Q.kplus1 ~ A1 + L1', Y2='Q.kplus1 ~ A2 + L2', Y3='Q.kplus1 ~ A3 + L3'), 
                gcomp=TRUE, SL.library=NULL, survivalOutcome = TRUE, id=ids)
result$estimates  # Treat-all = 0.5634683

result <- ltmle(data, Anodes=c("A1", "A2", "A3"), Cnodes=NULL, Lnodes=c("L1", "L2", "L3"), Ynodes=c("Y1", "Y2", "Y3"), abar=c(1, 1, 1),
                Qform=c(Y1='Q.kplus1 ~ A1 + L1', Y2='Q.kplus1 ~ A2 + A1 + L2', Y3='Q.kplus1 ~ A3 + A2 + L3'), 
                gcomp=TRUE, SL.library=NULL, survivalOutcome = TRUE, id=ids)
result$estimates  # Treat-all = 0.4334696

result <- ltmle(data, Anodes=c("A1", "A2", "A3"), Cnodes=NULL, Lnodes=c("L1", "L2", "L3"), Ynodes=c("Y1", "Y2", "Y3"), abar=c(0, 0, 0),
                Qform=c(Y1='Q.kplus1 ~ A1 + L1', Y2='Q.kplus1 ~ A2 + A1 + L2', Y3='Q.kplus1 ~ A3 + A2 + L3'), 
                gcomp=TRUE, SL.library=NULL, survivalOutcome = TRUE, id=ids)
result$estimates  # Treat-all = 0.6282985

### G-estimation Tests ###
library(DTRreg)

data <- read.csv(file="C:/Users/zivic/Desktop/data.csv", header=TRUE, sep=",")

# Continuous: 1 parameter
snm1 <- DTRreg(cd4_wk45, list(~ 1), 
               list(art ~ male + age0 + age_sq + age_cu + cd40 + cd4_sq + cd4_cu + dvl0), 
               list(~ 1), data=data, method='gest')
summary(snm1)

# Continuous: 2 parameter
snm1 <- DTRreg(cd4_wk45, list(~ male), 
               list(art ~ male + age0 + age_sq + age_cu + cd40 + cd4_sq + cd4_cu + dvl0), 
               list(~ 1), data=data, method='gest')
summary(snm1)

# Continuous: 3 parameter
snm1 <- DTRreg(cd4_wk45, list(~ male + cd40), 
               list(art ~ male + age0 + age_sq + age_cu + cd40 + cd4_sq + cd4_cu + dvl0), 
               list(~ 1), data=data, method='gest')
summary(snm1)

# Continuous: 5 parameter
snm1 <- DTRreg(cd4_wk45, list(~ male + dvl0 + cd40 + age0), 
               list(art ~ male + age0 + age_sq + age_cu + cd40 + cd4_sq + cd4_cu + dvl0), 
               list(~ 1), data=data, method='gest')
summary(snm1)


# Binary: 1 parameter
snm1 <- DTRreg(dead, list(~ 1), 
               list(art ~ male + age0 + age_sq + age_cu + cd40 + cd4_sq + cd4_cu + dvl0), 
               list(~ 1), data=data, method='gest')
summary(snm1)

# Binary: 2 parameter
snm1 <- DTRreg(dead, list(~ male), 
               list(art ~ male + age0 + age_sq + age_cu + cd40 + cd4_sq + cd4_cu + dvl0), 
               list(~ 1), data=data, method='gest')
summary(snm1)

# Binary: 3 parameter
snm1 <- DTRreg(dead, list(~ male + cd40), 
               list(art ~ male + age0 + age_sq + age_cu + cd40 + cd4_sq + cd4_cu + dvl0), 
               list(~ 1), data=data, method='gest')
summary(snm1)

# Binary: 5 parameter
snm1 <- DTRreg(dead, list(~ male + dvl0 + cd40 + age0), 
               list(art ~ male + age0 + age_sq + age_cu + cd40 + cd4_sq + cd4_cu + dvl0), 
               list(~ 1), data=data, method='gest')
summary(snm1)


### SMD ###
library(stddiff)

treat<-c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
category<-c(1, 2, 3, 1, 2, 3, 1, 3, 2, 1, 2, 3, 2, 1)
data<-data.frame(treat,category)
stddiff.category(data=data,gcol=1,vcol=c(2,2))
