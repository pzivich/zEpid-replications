##################################################################
# R-code for testing/comparisons I have made with zEpid
##################################################################
# Last edit; 2019/2/16


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
