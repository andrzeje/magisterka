

numFactors <- 3           
numMaturities <- length(maturities)
lenTseries <- nrow(yields)


lambda0 <- 0.0609

X <- cbind(matrix(1,length(maturities)), 
           (1-exp(-lambda0*maturities))/(lambda0*maturities), 
           ((1-exp(-lambda0*maturities))/(lambda0*maturities)-exp(-lambda0*maturities)))
beta <- matrix(0,nrow(yields),3)

residuals <- matrix(0,nrow(yields),length(maturities))

for (i in 1:nrow(yields)){
  olsfit <- lm(unlist(yields[i,])~X-1,na.action = na.exclude)
  beta[i,] <- unname(olsfit$coefficients)
  residuals[i,] <- unname(olsfit$residuals)
}

EstMdlVAR <- ar(beta, order = 1, method = "ols")

A0 <-matrix(EstMdlVAR$ar,c(3,3))

Q0 <- matrix(EstMdlVAR$var.pred,c(3,3))
#Q0 <- t(chol(Q0))
Q0 <- sqrt(c(Q0[1],0,0,Q0[5],0,Q0[9]))

H0 <- cov(residuals)
H0 <- diag(H0)
H0 <- H0[H0>0]
H0 <- sqrt(H0)

mu0 <- colMeans(beta)

param0 <- list(A0, Q0, H0, mu0, lambda0)
param0 <-unlist(param0)

numParams <- length(param0)
#param0_matlab <- c(0.990104430123422, -0.0281240136916211, 0.0517849291657113, 0.0249684220940171, 0.942561543530871,
#0.0124733219740373, -0.00229431907867095, 0.0286993873750660, 0.788078795496486, 0.338906023889202, 0, 0, 
#0.627898385454499, 0, 1.10237105011337, 0.141709398459525, 0.0728948481276874, 0.114923387137557, 0.111200084878294,
#0.0905579451691759, 0.0767207451043713, 0.0722210755302052, 0.0707643092581408, 0.0701289067342859, 0.0726736638115494,
#0.106242058067296, 0.0902962131596956, 0.103745274752660, 0.0980121532955626, 0.0912201363438591, 0.117941897636865,
#0.133544184092624, 8.34544392169689, -1.57244152304139, 0.202991858616548, 0.0609000000000000)

#param<-param0
#yield<-yields
#maturity<-maturities

library(KFAS)
Build_DieboldLi <- function(param, yields, maturities){
  yield <- yields
  maturity <-  maturities
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
  Q <- matrix(c(Qvalues[1],Qvalues[2],Qvalues[3],0,Qvalues[4],Qvalues[5],0,0,Qvalues[6]),numFactors,numFactors)
  Q <- Q %*% t(Q) 
  
  H <- diag(param[(offsetQ+1):offsetH],numMaturities,numMaturities)
  H <- H %*% H
  
  m <- matrix(param[(offsetH+1):(offsetH+numFactors)],c(numFactors))
  
  intercept <- C %*% m
  deflatedYield <- yield - matrix(t(intercept)[col(yield)],lenTseries,numMaturities)
  
  #P0<-matrix(c(4.33110517514536,	-0.352953505542163,	0.674691362241432,
  #-0.352953505542163,	3.96117338253019,	0.794438604581112,
  #0.674691362241433,	0.794438604581112,	2.95707829241291),3,3,byrow=T)
  
  P0 <- matrix(solve((diag(1,9) - A %x% A)) %*% matrix(Q,9,1),3,3) #page 138 - Time Series Analysis by State Space Methods; Second Edition, J. Durbin, S. J. Koopman
  
  #P0 <- diag(100,3)
  
  P1inf <- diag(0,3,3)
  a1 <- matrix(0,3)
  return(SSModel(as.matrix(deflatedYield) ~ -1 + SSMcustom(P1 = P0, a1 = a1, Z = as.matrix(C), T = as.matrix(A), R =  matrix(c(1,0,0,0,1,0,0,0,1),3,3), Q = as.matrix(Q)), H=as.matrix(H)))
}

initial_model <- Build_DieboldLi(param = param0, yields = yields, maturities = maturities)

updatefn<- function(pars, model, yields, maturities){
  new_model <- Build_DieboldLi(param = pars, yields = yields, maturities = maturities)
  model$y <- new_model$y
  model$Z <- new_model$Z
  model$T <- new_model$T
  model$Q <- new_model$Q
  model$H <- new_model$H
  return(model)
}
logLik(initial_model)
fitted_model <- fitSSM(initial_model, param0, updatefn, update_args = list(yields = yields, maturities = maturities), 
                       method = "BFGS", control = list(trace = 6, maxit = 100000, ndeps=rep(1e-8,numParams)))

#pars=param0
calc_loglik <- function(pars, yields,maturities){
  pars[numParams-3]<-pars[numParams-3]*10
  pars[numParams-2]<-pars[numParams-2]*10
  pars[numParams-1]<-pars[numParams-1]*10
  new_model <- Build_DieboldLi(param = pars, yields = yields, maturities = maturities)
  cost <- -logLik(new_model)
  #print(cost)
  #print(pars)
  if (is.na(cost)){
    return(9999)
  } else{
    return(cost)
  }
}
library(pso)
#param0[15]<-0.99999999
param0[numParams-3]<-param0[numParams-3]/10
param0[numParams-2]<-param0[numParams-2]/10
param0[numParams-1]<-param0[numParams-1]/10
est_param <- psoptim(par = param0, fn = calc_loglik, control = list(trace = 6, s=100), yields = yields, maturities = maturities)
 est_param <- optim(par = param0, fn = calc_loglik, method = "BFGS", control = list(trace = 6, maxit = 1000), yields = yields, maturities = maturities)

est_param$par[numParams-3] = est_param$par[numParams-3]*10
est_param$par[numParams-2] = est_param$par[numParams-2]*10
est_param$par[numParams-1] = est_param$par[numParams-1]*10
estimated_model <- Build_DieboldLi(param = est_param$par, yields = yields, maturities = maturities)
tt <- matrix(estimated_model$y,c(348,17)) + (t(matrix(estimated_model$Z,c(17,3)) %*% matrix(c(est_param$par[33],est_param$par[34],est_param$par[35]))))[col(matrix(estimated_model$y,c(348,17)))]
estimated_parameters <- est_param$par

estimated_model <- fitted_model$model
estimated_parameters <- fitted_model$optim.out$par

#matlab final loglik: 3184.553024590277  
#matlab estimated parameters:
# params <- c(0.994379534, -0.029000675, 0.025273442, 0.028596299, 0.939107357, 0.022918016, -0.022113864, 0.039582824,
#             0.841468145, 0.307593089, -0.045278948, 0.142056958, 0.616976663, 0.025480271, 0.882416067, 0.268233668,
#             0.075493268, 0.090294723, 0.104508294, 0.099150289, 0.0864779, 0.078629986, 0.072091124, 0.07268272,
#             0.079096104, 0.102954872, 0.092607469, 0.100415228, 0.111761682, 0.106970331, 0.15070029, 0.172784562,
#             8.024602633, -1.442312972, -0.418834565, 0.077764082)

mu <- matrix(c(estimated_parameters[numParams-3],estimated_parameters[numParams-2],estimated_parameters[numParams-1]),c(3,1))
intercept <- matrix(estimated_model$Z,c(numMaturities,3)) %*% mu 
deflatedYields <- yields - matrix(t(intercept)[col(yields)],lenTseries,numMaturities)
KFS_done <- KFS(estimated_model)
deflatedStates <- KFS_done$att
estimatedStates <- deflatedStates + matrix(t(mu)[col(deflatedStates)],c(lenTseries,3))

colMeans(estimatedStates)
plot(estimatedStates[,1], col="blue", type="l")
lines(beta[,1], col="red") 
mean(estimatedStates[,1,drop=F]-beta[,1,drop=F])

plot(estimatedStates[,2], col="blue", type="l")
lines(beta[,2], col="red") 
mean(estimatedStates[,2]-beta[,2])

plot(estimatedStates[,3], col="blue", type="l")
lines(beta[,3], col="red") 
mean(estimatedStates[,3]-beta[,3])


lambda <- estimated_parameters[numParams]

tau <- 0:max(maturities)   
decay = c(lambda0, lambda)
loading <- matrix(0,length(decay),length(tau))
for (i in 1:2){
  loading[i,] <- ((1-exp(-decay[i]*tau))/(decay[i]*tau)-exp(-decay[i]*tau))
}

plot(loading[1,], col="blue", type="l")
lines(loading[2,], col="red") 


residualsSSM <- yields - t(matrix(estimated_model$Z,c(numMaturities,3)) %*% t(unname(estimatedStates)))
residuals2Step <- yields - t(X %*% t(beta))

residualMeanSSM <- 100*colMeans(residualsSSM)
residualStdSSM <- 100*apply(residualsSSM, 2, sd)
residualMean2Step <- 100*colMeans(residuals2Step)
residualStd2Step <- 100*apply(residuals2Step, 2, sd)

abs(residualMean2Step)-abs(residualMeanSSM)

sum(residualsSSM^2)
sum(residuals2Step^2)

forecastDeflated <- predict(estimated_model,n.ahead = 1)
a=c()
for(i in 1:numMaturities){
  new<-forecastDeflated[[i]]
  a<-cbind(a,new)
}

plot(unname(a)[1,]+intercept,type='l')

plot(unlist(unname(yields)[lenTseries,]),type='l')

