library(R.matlab)
getwd()
setwd("C:\\Users\\apalu\\Desktop\\mgr\\72107")
dane <-readMat("Data_DieboldLi.mat")

maturities <- dane$maturities
maturities <-matrix(maturities)

yields <- as.data.frame(dane$Data)

lambda0 <- 0.0609

X <- cbind(matrix(1,length(maturities)), 
           (1-exp(-lambda0*maturities))/(lambda0*maturities), 
           ((1-exp(-lambda0*maturities))/(lambda0*maturities)-exp(-lambda0*maturities)))
beta <- matrix(0,nrow(yields),3)

residuals <- matrix(0,nrow(yields),length(maturities))

for (i in 1:nrow(yields)){
  olsfit <- lm(unlist(yields[i,])~X-1)
  beta[i,] <- unname(olsfit$coefficients)
  residuals[i,] <- unname(olsfit$residuals)
}

EstMdlVAR <- ar(beta, order = 1, method = "ols")

A0 <-matrix(EstMdlVAR$ar,c(3,3))

Q0 <- matrix(EstMdlVAR$var.pred,c(3,3))
#Q0 <- t(chol(Q0))
Q0 <- rep(1,9)
Q0 <- c(log(Q0[1]),0,0,log(Q0[5]),0,log(Q0[9]))

H0 <- cov(residuals)
H0 <- diag(H0)
H0 <- H0[H0>0]
H0 <- sqrt(H0)

H0<-log(rep(1,17))

mu0 <- colMeans(beta)

Cov0 <- cov(beta)
Cov0 <- t(chol(Cov0))
Cov0 <- c(Cov0[1],Cov0[2],Cov0[3],Cov0[5],Cov0[6],Cov0[9])

param0 <- list(H0)
param0 <-unlist(param0)

param<-param0
yield<-yields
maturity<-maturities

library(KFAS)
Build_DieboldLi <- function(param, yields){
  yield <- yields
  maturity <-  c(3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120)
  numFactors <- 3           
  numMaturities <- length(maturity)
  
  lambda <- 0.0778
  
  C <- cbind(matrix(1,length(maturity)), 
             (1-exp(-lambda*maturity))/(lambda*maturity), 
             (((1-exp(-lambda*maturity))/(lambda*maturity))-exp(-lambda*maturity)))
  
  C<-matrix(C,17,3)
  
  A <- matrix(c(0.9944,    0.0286,   -0.0221,
  -0.0290,    0.9391,    0.0396,
  0.0253,  0.0229,    0.8415),c(3,3),byrow = T)
  
  Q <- matrix(c(0.1149,   -0.0266   ,-0.0719,
                -0.0266,    0.3943 ,   0.0140,
                -0.0719 ,   0.0140,    1.2152),
              c(3,3),byrow = T)
  
  H <- diag(exp(param[1:17]),numMaturities,numMaturities)
  #H <- exp(H) %*% H
  
  m <- matrix(c(8.0246,   -1.4423,   -0.4188),c(numFactors))
  
  intercept <- C %*% m
  deflatedYield <- yield - matrix(t(intercept)[col(yield)],348,17)
  
  return(SSModel(as.matrix(deflatedYield) ~ -1 + SSMcustom(Z = as.matrix(C), T = as.matrix(A), R =  matrix(c(1,0,0,0,1,0,0,0,1),3,3), Q = as.matrix(Q)), H=as.matrix(H)))
}

initial_model <- Build_DieboldLi(param = param0, yields = yields)
tt <- matrix(initial_model$y,c(348,17)) + (t(matrix(initial_model$Z,c(17,3)) %*% matrix(c(param0[33],param0[34],param0[35]))))[col(matrix(initial_model$y,c(348,17)))]
colSums(yields - tt)

updatefn<- function(pars, model, yields){
  new_model <- Build_DieboldLi(param = pars, yields = yields)
  model$y <- new_model$y
  model$Z <- new_model$Z
  model$T <- new_model$T
  model$Q <- new_model$Q
  model$H <- new_model$H
  return(model)
}
param0<-fitted_model$optim.out$par
fitted_model <- fitSSM(initial_model, param0, updatefn, update_args = list(yields = yields), 
                       method = "L-BFGS-B", control = list(trace = 6, maxit = 1000))

pars=param0
calc_loglik <- function(pars, yields){
  pars[33]<-pars[33]*10
  pars[34]<-pars[34]*10
  pars[35]<-pars[35]*10
  new_model <- Build_DieboldLi(param = pars, yields = yields)
  cost <- -logLik(new_model)
  #print(cost)
  if (is.na(cost)){
    return(9999)
  } else{
    return(cost)
  }
}
library(pso)
param0[15]<-0.99999999
param0[33]<-param0[33]/10
param0[34]<-param0[34]/10
param0[35]<-param0[35]/10
est_param <- psoptim(par = param0, fn = calc_loglik, control = list(trace = 6, maxit=10000), yields = yields)
est_param <- optim(par = param0, fn = calc_loglik, method = "L-BFGS-B", control = list(trace = 6, maxit = 10), yields = yields)

est_param$par[33] = est_param$par[33]*10
est_param$par[34] = est_param$par[34]*10
est_param$par[35] = est_param$par[35]*10
estimated_model <- Build_DieboldLi(param = est_param$par, yields = yields)
tt <- matrix(estimated_model$y,c(348,17)) + (t(matrix(estimated_model$Z,c(17,3)) %*% matrix(c(est_param$par[33],est_param$par[34],est_param$par[35]))))[col(matrix(estimated_model$y,c(348,17)))]


mu <- matrix(c(est_param$par[33],est_param$par[34],est_param$par[35]),c(3,1))
intercept <- matrix(estimated_model$Z,c(17,3)) %*% mu 
deflatedYields <- yields - matrix(t(intercept)[col(yield)],348,17)
deflatedStates <- KFS(estimated_model)
estimatedStates <- deflatedStates$a + matrix(t(mu)[col(deflatedStates$a)],c(349,3))

colMeans(estimatedStates)
plot(estimatedStates[,1], col="blue")
lines(beta[,1], col="red") 
mean(estimatedStates[2:349,1]-beta[,1])

plot(estimatedStates[2:349,2], col="blue", type="l")
lines(beta[,2], col="red") 
mean(estimatedStates[2:349,2]-beta[,2])

plot(estimatedStates[2:349,3], col="blue", type="l")
lines(beta[,3], col="red") 
mean(estimatedStates[2:349,3]-beta[,3])


lembda <- est_param$par[36]

tau <- 0:max(maturities)   
decay = c(lambda0, lambda)
loading <- matrix(0,length(decay),length(tau))
for (i in 1:2){
  loading[i,] <- ((1-exp(-decay[i]*tau))/(decay[i]*tau)-exp(-decay[i]*tau))
}
plot(loading[1,], col="blue", type="l")
lines(loading[2,], col="red") 


residualsSSM <- yields - t(matrix(estimated_model$Z,c(17,3)) %*% t(unname(estimatedStates)))
residuals2Step <- yields - t(X %*% t(beta))

residualMeanSSM <- 100*colMeans(residualsSSM)
residualStdSSM <- 100*apply(residualsSSM, 2, sd)
residualMean2Step <- 100*colMeans(residuals2Step)
residualStd2Step <- 100*apply(residuals2Step, 2, sd)

sum(residualMeanSSM^2)
sum(residualMean2Step^2)

library(ggplot2)
p<-ggplot(estimatedStates[,1] aes(y=custom1))

geom_line())
p


#test publ wynikow
maturity <-  c(3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120)
lambda <- 0.0778
C <- cbind(matrix(1,length(maturity)), 
           (1-exp(-lambda*maturity))/(lambda*maturity), 
           (((1-exp(-lambda*maturity))/(lambda*maturity))-exp(-lambda*maturity)))
tt <- matrix(estimated_model$y,c(348,17)) + (t(matrix(C,c(17,3)) %*% matrix(c(8.0246,-1.4423,-0.4188))))[col(matrix(estimated_model$y,c(348,17)))]
colMeans(yields - tt)




estimatedStates <-  deflatedStates$s + opttt$m0[col(deflatedStates$s)]

plot(c(estimatedStates[,3]-beta[,3]),type = "l")
