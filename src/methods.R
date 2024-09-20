
DP <- function(Z, y, tau){
  n = dim(Z)[1]; p = dim(Z)[2]
  asd = CVglasso(Z, tau = tau, cores = 5)
  Omega= asd$Omega 

  Znew=gknockoffX(Z,Omega)

  s=1/eigen(Omega)$value[1]
  GZ=t(Znew)%*%Znew/n
  subome=Omega
  jiu=(diag(p)-subome%*%(s*diag(p)))
  ER11=tau*tau*diag(p)
  ER12=ER11%*%jiu
  ER21=t(jiu)%*%ER11
  ER22=t(jiu)%*%ER11%*%jiu
  ER1=cbind(ER11,ER12)
  ER2=cbind(ER21,ER22)
  ER=rbind(ER1,ER2)
  
  pena = "lasso"
  reco2.fit=coco(Znew,y,n,2*p,tau=ER, noise="additive", block = F,  K=5, center.Z = F, scale.Z = F, center.y = F, scale.y = F, penalty=pena)
  betaco1=reco2.fit$beta.opt
  
  return(betaco1)
}

RANK <- function(Z, y){
  n = dim(Z)[1]; p = dim(Z)[2]
  asd = CVglasso(Z, tau = 0, cores = 5)
  Omega= asd$Omega 
  Znew=gknockoffX(Z,Omega)
  fit.lasso = cv.glmnet(Znew, y ,nfolds=5)
  ridge1 =glmnet(Znew, y)
  beta1<-t(coef(ridge1,s=fit.lasso$lambda.min)[-1])
  return(beta1)
  
}


BHq <- function(Z, y, tau){
  #or use CVglasso(Z, tau = tau, cores = 5)
  obj<- score.nodewiselasso(Z, wantTheta = TRUE,
                            verbose = FALSE,
                            lambdaseq = "quantile",
                            parallel = FALSE,
                            oldschool = FALSE,
                            lambdatuningfactor = 1,
                            cv.verbose = FALSE,
                            do.ZnZ = TRUE)
  Theta = obj$out
  
  #debiased cocolasso
  decoco.fit=coco(Z,y,n,p,tau=tau, noise="additive", block = F, K=5, center.Z = F, scale.Z = F, center.y = F, scale.y = F, penalty="lasso")
  beta_co=decoco.fit$beta.opt
  sigw = tau^2*diag(p)
  part2b = n^(-1)*t(Z)%*%y - (n^(-1)*t(Z)%*%Z - sigw)%*%beta_co
  bwh = beta_co + Theta%*%part2b
  
  Tau = matrix(0,p,p)
  for(i in 1:n){
    Tau11 = as.vector(Z[i,])%*%(y[i] - t(as.vector(Z[i,]))%*%beta_co) + sigw%*%beta_co
    Tau = Tau + n^(-1)*Tau11%*%t(Tau11)
  }
  
  Ome_wh = Theta%*%Tau%*%Theta
  w_wh = rep(0,p)
  pval = rep(0,p)
  for(j in 1:p){
    sej <- Ome_wh[j,j]/sqrt(n)
    w_wh[j] = bwh[j]/sej
    pval[j] <- 2 * pnorm(-abs(w_wh[j]))
  }
  return(pval)
}


DP_de <- function(Z, y, tau){
  n = dim(Z)[1]; p = dim(Z)[2]
  asd = CVglasso(Z, tau = tau, cores = 5)
  Omega= asd$Omega 
  
  Znew=gknockoffX(Z,Omega)
  
  s=1/eigen(Omega)$value[1]
  GZ=t(Znew)%*%Znew/n
  subome=Omega
  jiu=(diag(p)-subome%*%(s*diag(p)))
  ER11=tau*tau*diag(p)
  ER12=ER11%*%jiu
  ER21=t(jiu)%*%ER11
  ER22=t(jiu)%*%ER11%*%jiu
  ER1=cbind(ER11,ER12)
  ER2=cbind(ER21,ER22)
  ER=rbind(ER1,ER2)
  
  #########
  pena = "lasso"
  reco2.fit=coco(Znew,y,n,2*p,tau=ER, noise="additive", block=F, K=5, center.Z = F, scale.Z = F, center.y = F, scale.y = F, penalty=pena)
  beta_co2 =reco2.fit$beta.opt
  
  ###########debiased lasso
  asd = CVglasso(Znew, tau = tau, cores = 5)
  Theta2 = asd$Omega 

  sigw = ER
  part2b = n^(-1)*t(Znew)%*%y - (n^(-1)*t(Znew)%*%Znew - sigw)%*%beta_co2
  bwh22 = beta_co2 + Theta2%*%part2b
  ##also recommand scad after correcting measurement errors
  ##scad is an unbiased method, can replace the debiased lasso in some cases
  # CV_s=cv.ncvreg(tilZ, tily, penalty="SCAD",  nfolds = 10, eps = 1e-6)
  # beta=coef(CV_s, s = "lambda.min", exact=TRUE)[-1]
  # beta = as.numeric(beta)
  

  return(bwh22)

}

