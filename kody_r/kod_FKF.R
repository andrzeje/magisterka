#install.packages('R.matlab')
library(R.matlab)
getwd()
setwd("C:\\Users\\apalu\\Desktop\\mgr\\72107")
dane <-readMat("Data_DieboldLi.mat")

maturities <- dane$maturities
maturities <-matrix(maturities)

yields <- as.data.frame(dane$Data)

lambda0 = 0.0609
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

#library(vars)
#EstMdlVAR <- VAR(beta,1)
#plot(unlist(beta[,1]))
EstMdlVAR <- ar(beta, order = 1, method = "ols")

A0 <-matrix(EstMdlVAR$ar,c(3,3))

Q0 <- matrix(EstMdlVAR$var.pred,3,3)
#Q0 <- log(diag(Q0))
Q0 <- t(chol(Q0))

Q0 <- c((Q0[1]),Q0[2],Q0[3],(Q0[5]),Q0[6],(Q0[9]))
#Q0 <- rep(0,6)
param0 <- c(0.994379534, -0.029000675, 0.025273442, 0.028596299, 0.939107357, 0.022918016, -0.022113864, 0.039582824,
            0.841468145, 0.307593089, -0.045278948, 0.142056958, 0.616976663, 0.025480271, 0.882416067, 0.268233668,
            0.075493268, 0.090294723, 0.104508294, 0.099150289, 0.0864779, 0.078629986, 0.072091124, 0.07268272,
            0.079096104, 0.102954872, 0.092607469, 0.100415228, 0.111761682, 0.106970331, 0.15070029, 0.172784562,
            8.024602633, -1.442312972, -0.418834565, 0.077764082)
H0 <- cov(residuals)
H0 <- diag(H0)
H0 <- H0[H0>0]
#H0 <- log(H0)

H0 <- sqrt(H0)
#H0 <- rep(0,17)
mu0 <- colMeans(beta)
mu0 <- c(8.0246,   -1.4423,   -0.4188)

#Cov0 <- cov(beta)
#Cov0 <- t(chol(Cov0))
#Cov0 <- c(Cov0[1],Cov0[2],Cov0[3],Cov0[5],Cov0[6],Cov0[9])

#mu0, Cov0,

param0 <- list(A0, Q0, H0, mu0, lambda0)
param0 <-unlist(param0)

param<-param0
yield<-yields
maturity<-maturities
library(FKF)
Build_DieboldLi <- function(param, yields){
  maturity <-  c(3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120)
  numFactors <- 3           
  numMaturities <- length(maturity)
  lenTseries <- nrow(yields)
  
  lambda <- tail(param,1)
  
  C <- cbind(matrix(1,length(maturity)), 
             (1-exp(-lambda*maturity))/(lambda*maturity), 
             (((1-exp(-lambda*maturity))/(lambda*maturity))-exp(-lambda*maturity)))
  
  C<-matrix(C,numMaturities,numFactors)
  
  offsetA<-numFactors*numFactors
  offsetQ <- offsetA + numFactors*(numFactors+1)/2
  offsetH <- offsetQ + numMaturities
  
  A <- matrix(param[1:offsetA],numFactors,numFactors)
  
  Qvalues <- param[(offsetA+1):offsetQ]
  
  Q <- matrix(c((Qvalues[1]),(Qvalues[2]),(Qvalues[3]),
                (0),(Qvalues[4]),(Qvalues[5]),
                (0),(0),(Qvalues[6])),numFactors,numFactors)
  Q <- Q %*% t(Q) 
  
  H <- diag((param[(offsetQ+1):offsetH]),numMaturities,numMaturities)
  H <- H %*% H
  
  m <- matrix(param[(offsetH+1):(offsetH+numFactors)],c(numFactors))
  
  intercept <- C %*% m
  #deflatedYield <- yields - matrix(t(intercept)[col(yields)],lenTseries,numMaturities)
  
  m0<-matrix(0,numFactors)
  #diag(0,numFactors,numFactors)
  C0<-matrix(c(4.33110517514536,	-0.352953505542163,	0.674691362241432,
  -0.352953505542163,	3.96117338253019,	0.794438604581112,
  0.674691362241433,	0.794438604581112,	2.95707829241291),3,3)
  return(fkf(a0=c(0,0,0),P0=C0,dt=c(0,0,0),ct=-intercept,Tt=A,Zt=C,HHt=Q,GGt=H,yt=t(as.matrix(yields)))  )
  #return(list(dlm(m0=m0, C0=C0, FF = C, V = H, GG = A, W = Q),deflatedYield))
}
aa<-fkf(a0=c(0,0,0),P0=C0,dt=c(0,0,0),ct=-intercept,Tt=A,Zt=C,HHt=Q,GGt=H,yt=t(as.matrix(yields)))

library(FKF)

opt <- function(param,yields){
  #print(param)
  cost <- tryCatch(expr = {
    result <- Build_DieboldLi(param = param0, yields = yields)
    #dlmLL(as.matrix(result[[2]]),result[[1]])
    -result$logLik
  },
  error = function(e){print(e); return(9999)},
  warning = {9999}
  )
  #print(cost)
  if(!is.na(cost)){
    return(cost)
  }else{
      return(9999)
    }
}

opt(param0,yields)


library(pso)
aa <- psoptim(par = param0, fn = opt, lower = -10, upper = 10, control = list(trace = 6, maxit=1000), yields = yields)
results <- optim(par = param0, fn = opt, method = "L-BFGS-B", yields = yields,control = list(trace = 6, maxit=1000))
aa <- maxBHHH(fn = opt, start = param0, control = list(printLevel = 6), yields = yields)
par0 <- aa$par
param0[33]<-param0[33]/10
param0[34]<-param0[34]/10
param0[35]<-param0[35]/10
param0[15]<-0.99999


param0<-par0

lplp<-dlmMLE(as.matrix(yields), param0, Example_DieboldLi)
opttt <- Example_DieboldLi(param = aa$par)
opttt[[1]]
opttt$GG
opttt$W
opttt$m0
opttt$C0


colMeans(a$a)

intercept <- opttt[[1]]$FF %*% opttt[[1]]$m0
deflatedYields <- yields - t(intercept)
a <- dlmFilter(as.matrix(deflatedYields), opttt[[1]])
deflatedStates <- dlmSmooth(deflatedYields,opttt)

estimatedStates <-  deflatedStates$s + opttt$m0[col(deflatedStates$s)]

plot(c(estimatedStates[,3]-beta[,3]),type = "l")

