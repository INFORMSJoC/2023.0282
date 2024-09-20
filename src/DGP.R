library(MASS) 
library(Matrix)
library(pracma)
library(glasso)
library(mvtnorm)
library(glmnet)
library(knockoff)
library(glassoFast)
library(huge)
library(BDcocolasso)
library(parallel)
library(foreach)
library(doParallel)

DGP_simulation1 = function(n , p ){
  rho = 0.5
  ini_sigma=gen_sigma(p)
  X=mvrnorm(n,rep(0,p),ini_sigma)
  #generate beta
  k= 25
  nonzero = sample(p, k)
  eff=c(1.5,3,-1.5,-3)
  xx=rep(0,p)
  betastar=rep(0,p)
  for(i in 1:p){
    xx[i]=sample(eff,1)
  }
  betastar=xx*(1:p %in% nonzero)
  y= X %*% betastar + c(rnorm(n))
  A= tau * matrix(rnorm(n*p), n, p)
  Z= X + A
  
  return(list(Z=Z,
              y=y,
              beta = betastar))
  
}

DGP_simulation2 = function(n , p ){

  data_band=huge.generator(n=n,d=p,graph = "band")
  X=data_band$data
  #generate beta
  k= 25
  nonzero = sample(p, k)
  eff=c(1.5,3,-1.5,-3)
  xx=rep(0,p)
  betastar=rep(0,p)
  for(i in 1:p){
    xx[i]=sample(eff,1)
  }
  betastar=xx*(1:p %in% nonzero)
  y= X %*% betastar + c(rnorm(n))
  A= tau * matrix(rnorm(n*p), n, p)
  Z= X + A
  
  return(list(Z=Z,
              y=y,
              beta = betastar))
  
}


DGP_simulation3 = function(n, p ){
  rho = 0.5
  ini_sigma=gen_sigma(p)
  X=rmvt(n,ini_sigma,df=10) #t-distribution
  #generate beta
  k = 25
  nonzero = sample(p, k)
  eff=c(1.5,3,-1.5,-3)
  xx=rep(0,p)
  betastar=rep(0,p)
  for(i in 1:p){
    xx[i]=sample(eff,1)
  }
  betastar=xx*(1:p %in% nonzero)
  y= X %*% betastar + c(rnorm(n))
  
  A= tau * matrix(rnorm(n*p), n, p)
  Z= X + A
  
  return(list(Z=Z,
              y=y,
              beta = betastar))
  
}


### step 1: simulate 100 data set
runtimes = 100
data_Z <- list()
data_Y <- list()
data_B <- list()

n = 200
tau = 0.25

p = 400 # 500, 600
for (i in 1:runtimes) {
  data = DGP_simulation1(n , p)
  data_Z[[i]] = data$Z
  data_Y[[i]] = data$y
  data_B[[i]] = data$beta
  print(i)
}

for (i in 1:runtimes) {
  data = DGP_simulation2(n , p)
  data_Z[[i]] = data$Z
  data_Y[[i]] = data$y
  data_B[[i]] = data$beta
  print(i)
}

pp <- c(140, 180, 220, 260, 300)
p <- pp[1]
for (i in 1:runtimes) {
  data = DGP_simulation3(n , p)
  data_Z[[i]] = data$Z
  data_Y[[i]] = data$y
  data_B[[i]] = data$beta
  print(i)
}


