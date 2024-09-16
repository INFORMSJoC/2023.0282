##########part 1: read the raw data and use bgx to generate mean and std
if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
BiocManager::install(version = "3.13")
BiocManager::install("bgx")
BiocManager::install("hgu133plus2cdf")
BiocManager::install("affy")

library(BiocManager)
library(bgx)
library(affy)
library(hgu133plus2cdf)

bdata<-ReadAffy() #input your file path containing all the .CEL files such as "Pasient2.CEL"
bset<-bgx(bdata)
va<-assayData(bset)
mean<-va$exprs
std<-va$se.exprs

da<-t(mean)
rownames(da)=NULL

da2<-t(std)
rownames(da2)=NULL

##########part 2: use RMA and SIS to screening

bgcorrect.methods=c(bgcorrect.methods,"rma001","rma002")

b = data 
#This obtains an AffyBatch object from your CEL files
a = expresso(b, bgcorrect.method="rma",
             normalize.method="quantiles", pmcorrect.method="pmonly",
             summary.method="medianpolish")
#RMA-Mean
aa=exprs(a)
exprs(a)=2^aa
aaa=exprs(a)

##response: the total hip T-score obtained from the data website
y<-c(-4.1,-3,-2.1,0.6,-2.9,0.5,0.4,-2.4,-2.1,0.1,0.7,-2.3,-2.9,-3.2,-2.8,-0.1,1.9,1.6,-3.8,-0.6,-2.6,-2.6,
     2,-1.3,-0.9,1.1,0.2,-1.4,-3.6,-2,-4,-1.4,0.5,-1.6,0.6,-1.8,0.3,-2.7,0.1,1.8,0.3,-2.5,-3,0.1,-2.2,-0.9,-0.2,-1.7,-2.1,0.7,0.1,
     0.5,-0.6,-0.1,-0.1,0.6,-2.7,-0.7,-1.4,1.7,-2.2,-1.7,-1.7,-1.6,-1.3,-0.9,0.3,0.2,1,1.7,-0.1,0.3,-0.8,-0.8,2,-2.7,-2.8,
     -1.5,0.8,0.2,0.6,1.2,-1.1,-3.7)

tt2=as.matrix(y)
model=SIS(t(aaa), tt2, family='gaussian', penalty = "lasso", tune = "bic",  varISIS='cons', nsis=1800)
ix=model$sis.ix0

##########part 3: calculate measurement errors and final data

da=da[,ix]
da2=da2[,ix]

n=nrow(da)
p=ncol(da)


#W_bar
W_bar=rep(0,p)
for(j in 1:p){
  W_bar[j]=sum(da[ ,j])/n
}

#s
s2=rep(0,p)
for(j in 1:p){
  s2[j] = var(da[,j])*(n-1)/n
}
s=sqrt(s2)
#W
W=da
for(i in 1:n){
  for(j in 1:p){
    W[i,j]= (da[i,j]- W_bar[j] )/s[j]
  }
}
#var
#da2=da2^2 #input is se
sigma_u2=rep(0,p)
Sig_u=rep(0,p)
for(j in 1:p){
  sigma_u2[j] = sum(da2[,j])/n
  Sig_u[j] = sigma_u2[j]/s2[j]
}
index = which(sigma_u2 < 0.5*s2)
ix0=index
Z=W[,ix0]
Sigma_A= diag(Sig_u[ix0])

##obtained Z, y, Sigma_A, 

library(graphics)
library(MASS) 
library(Matrix)
library(pracma)
library(scalreg)
library(glasso)
library(mvtnorm)
library(msgps)
library(glmnet)
library(knockoff)
library(glassoFast)
library(BDcocolasso)

y0 = y
Z0 = Z
Sigma_A0 = Sigma_A

n = dim(Z0)[1]; pz = p = dim(Z0)[2]

iter = 100
relevant <- list()
error <- list()

for (ii in 1:iter) {
  n = dim(Z0)[1]
  squta = 1 : n
  observation_indices = sample(1:n, 64)
  index_prediction = squta[-observation_indices]
  
  
  X_mid <- Z0[observation_indices,]
  y_mid <- as.matrix(y0[observation_indices,])
  SigmaA_tr <- Sigma_A0
  
  # ---- trainning ----
  relevant_features <- list()
  
  g_ome = cv.glasso_data(X_mid, SigmaA_tr,  0)
  Omega = g_ome$wi
  
  repo_rank <- Rank_nopre(X_mid, y_mid, SigmaA_tr, Omega = Omega  )
  relevant_features[[1]] <- repo_rank$S
  relevant_features[[2]] <- repo_rank$Splus

  
  #or use CVglasso(Z, tau = Sigma_A0)
  obj<- score.nodewiselasso(X_mid, wantTheta = TRUE,
                            verbose = FALSE,
                            lambdaseq = "quantile",
                            parallel = FALSE,
                            oldschool = FALSE,
                            lambdatuningfactor = 1,
                            cv.verbose = FALSE,
                            do.ZnZ = TRUE)
  Omega = obj$out
  #Omega = g_ome$wi
  
  repo <- DCocolasso(X_mid, y_mid, Sigma_A0, Omega = Omega   )
  relevant_features[[3]] <- repo$ldco
  
  

  g_ome = cv.glasso_data(X_mid, Sigma_A0)
  Omega = g_ome$wi
  # asd = CVglasso(X_mid, tau = Sigma_A0)
  # Omega= asd$Omega
  repo_dp <- DPknockoff_nopre(X_mid, y_mid, SigmaA_tr, Omega = Omega)

  
  relevant_features[[4]] <- repo_dp$S
  relevant_features[[5]] <- repo_dp$Splus
  
  # ---- use the trainned model and X_inde to predict y
  X_inde <- Z0[(index_prediction), ]
  y <-as.matrix(y0[(index_prediction)])
  
  #rank
  error_1 <- linear_prediction_loss(relevant_features[1:10], X_mid, y_mid, X_inde, y)
  #error_1
  #dcoco and DP
  error_2 <- linear_prediction_loss2(relevant_features[-(1:10)], X_mid, y_mid, X_inde, y, SigmaA_tr)
  #error_2
  error_t = c(error_1,error_2)
  # error_t
  # relevant_features
  
  error[[ii]] = error_t
  relevant[[ii]] = relevant_features
  
}

ftemp=NULL
for(i in 1:iter){
  dtemp<-paste('iter',i,sep = '')
  ftemp<-c(ftemp,dtemp)
}
dim1 = ftemp
dim2=c("rank",
       "rank+",
       "deco",
       "DP", 
       "DP+")

err_final = array(0, dim = c(iter, 5), dimnames = list(dim1,dim2))
rr_final = array(0, dim = c(iter, 5), dimnames = list(dim1,dim2))
set_final = array(0, dim = c(iter, 5), dimnames = list(dim1,dim2))

for(i in 1:iter){
  for(j in 1:5){
    temp <- error[[i]][[j]]
    err_final[i,j] <- temp[1]
    set_final[i,j] <- length(relevant[[i]][[j]])
    
  }
}

std <- function(x) sd(x)/sqrt(length(x))


colMeans(set_final)
colMeans(err_final/20)
apply(err_final/20, 2, std)
apply(set_final, 2, std)


#freqcency of selected variables
method_name= c("rank",
               "rank+",
               "deco",
               "DP", 
               "DP+")
var_name = colnames(Z0)
sele_all =  array(0, dim = c(5, 319), dimnames = list(method_name,var_name))
for(i in 1:iter){
  for(j in 1:5){
    temp <- relevant[[i]][[j]]
    sele_all[j, temp] = sele_all[j, temp] + 1
  }
}

for(i in 1:5){
  aaa = sort(sele_all[i,], decreasing = TRUE)
  print(aaa[1:5])
}



# prediction funtion and feature selection methods ----
linear_prediction_loss <- function(relevant_features, X_mid, y_mid, x, y) {
  
  # lm(y_mid~X_mid[, numeric(0)])
  result <- list()
  #prediction_error <- c()
  for (i in 1:length(relevant_features)) {
    if (length(relevant_features[[i]]) > 0) {
      #coef <- lm(y_mid ~ X_mid[, relevant_features[[i]]])$coef
      tmp = X_mid[, relevant_features[[i]]]
      coef1 <- solve(t(tmp)%*%tmp)%*%t(tmp)%*%y_mid
      intercept <- mean(y_mid) - mean(tmp%*%coef1)
      coef <- c(intercept, coef1)
      
    } else {
      coef <- lm(y_mid ~ 1)$coef
      #cat(coef, '--372\n')
    }
    
    SSE =  sum((c(   cbind(1, x[,relevant_features[[i]]] ) %*% coef ) -  y )^2)
    SS = sum((y - mean(y))^2)
    R2 = 1 - SSE/SS
    
    #RMSE = sqrt(SSE/length(y))
    result[[i]] <- c(SSE, R2)
    
    #prediction_error <- c(prediction_error, (c(c(1, x[relevant_features[[i]]]) %*% coef) - y)^(2) )
  }
  return(result)
}

linear_prediction_loss2 <- function(relevant_features, X_mid, y_mid, x, y, tau) {
  
  result <- list()
  #prediction_error <- c()
  for (i in 1:length(relevant_features)) {
    if (length(relevant_features[[i]]) > 0) {
      #coef <- lm(y_mid ~ X_mid[, relevant_features[[i]]])$coef
      nonzero = relevant_features[[i]]
      err <- tau[nonzero,nonzero]
      tmp = X_mid[, relevant_features[[i]]]
      gram = t(tmp)%*%tmp - err
      coef1 <- solve(gram)%*%t(tmp)%*%y_mid
      intercept <- mean(y_mid) - mean(tmp%*%coef1)
      coef <- c(intercept, coef1)
      
    } else {
      coef <- lm(y_mid ~ 1)$coef
      #cat(coef, '--372\n')
    }
    
    SSE =  sum((c(   cbind(1, x[,relevant_features[[i]]] ) %*% coef ) -  y )^2)
    SS = sum((y - mean(y))^2)
    R2 = 1 - SSE/SS
    
    #RMSE = sqrt(SSE/length(y))
    result[[i]] <- c(SSE, R2)
    
    #prediction_error <- c(prediction_error, (c(c(1, x[relevant_features[[i]]]) %*% coef) - y)^(2) )
  }
  return(result)
}



DCocolasso <- function(X, y, Sigma_A, Omega = NULL) {
  
  n = dim(X)[1]; p = pz = dim(X)[2]
  Z = X
  Theta = Omega
  
  
  reco=coco(Z,y,n,pz,tau=Sigma_A, noise="additive", block=F, center.Z = F, scale.Z = F, center.y = F, scale.y = F, penalty="lasso")
  beta_coco=reco$beta.opt
  beta_co = beta_coco
  tau = Sigma_A
  sigw = Sigma_A
  part2b = n^(-1)*t(Z)%*%y - (n^(-1)*t(Z)%*%Z - sigw)%*%beta_co
  bwh = beta_co + Theta%*%part2b
  
  #stanard w
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
  
  
  names(pval) = c(1:p)
  S_dco = list()
  q = 0.2
  res = bhq(pval, q)
  num_rej = res$Rejections
  ad.pva = res$Adjusted.pvalues
  
  if(num_rej != 0 ){
    sele = ad.pva[1:num_rej]
    num = rep(0, num_rej)
    for(j in 1:num_rej){
      num[j] = which(names(pval) == names(sele)[j])
    }
    S_dco = num
  }else{
    S_dco = which(1>2) #integer 0
  }
  
  return(list(ldco = S_dco))
}


Rank_nopre <- function(X, y, Sigma_A,  Omega = NULL) {
  n = dim(X)[1]; p = pz = dim(X)[2]
  Sigma_A = Sigma_A
  Z = X
  
  Znew=gknockoffX(Z,Omega)
  
  CV_la=cv.glmnet(Znew,  y,  alpha = 0.9,  nfolds = 10)
  beta_la=coef(CV_la, s = "lambda.min", exact=TRUE)[-1]
  beta2= as.numeric(beta_la)
  
  W2 = rep(0, pz)
  for (j in 1:pz) {
    W2[j] = abs(beta2[j]) - abs(beta2[j+pz])
  }
  
  q = 0.2
  S1 = S2 = list()
  
  t_rank=knockoff.threshold(W2, fdr  =q, offset = 0)
  S1=which(W2>=t_rank)
  
  t_rank2=knockoff.threshold(W2, fdr  =q, offset = 1)
  S2=which(W2>=t_rank2)
  
  
  
  return(list(S = S1, Splus = S2))
}



DPknockoff_nopre <- function(X, y, Sigma_A, Omega = NULL) {
  n = dim(X)[1]; p = pz = dim(X)[2]
  Z = X
  Znew=gknockoffX(Z,Omega)
  s= as.numeric(1/eigen(Omega)$value[1])
  GZ=t(Znew)%*%Znew/n
  subome=Omega
  jiu=(diag(pz)-subome%*%(s*diag(pz)))
  ER11=Sigma_A
  ER12=ER11%*%jiu
  ER21=t(jiu)%*%ER11
  ER22=t(jiu)%*%ER11%*%jiu
  ER1=cbind(ER11,ER12)
  ER2=cbind(ER21,ER22)
  ER=rbind(ER1,ER2)

  
  reco=coco(Znew,y,n,2*p,tau=ER, noise="additive", block=F, center.Z = F, scale.Z = F, center.y = F, scale.y = F, penalty="lasso")
  beta=reco$beta.opt
  
  W2 = c(rep(0, pz))
  for (j in 1:pz) {
    W2[j] = abs(beta[j]) - abs(beta[j+pz])
  }
  
  q = 0.2
  S1 = S2 = list()
  
  t_rank=knockoff.threshold(W2, fdr  =q, offset = 0)
  S1=which(W2>=t_rank)
  
  t_rank2=knockoff.threshold(W2, fdr  =q, offset = 1)
  S2=which(W2>=t_rank2)
  
  
  return(list(S = S1, Splus = S2))
}


cv.glasso_data=function(Z, Sigma_A, noise){
  
  n=nrow(Z)
  p=ncol(Z)
  ntau=20
  n.split=4
  split.ratio=1-1/log(n)
  n.train = ceiling(n*split.ratio)
  n.test = n-n.train
  tau.min = 0.5
  tau.max = 1
  tau.path = seq(tau.min,tau.max,length.out=ntau)*sqrt(log(p)/n)
  loss.re <- matrix(0, n.split, length(tau.path))
  if(noise == 1){
    for(i in 1:n.split){
      tr <- sample(1:n, round(n/n.split))
      S.tr <- t(Z[-tr,])%*%Z[-tr,]/(n-length(tr))- Sigma_A
      
      if(min(eigen(S.tr)$value) < 1e-3)
      {
        S.tr=maxproj.cov(mat=S.tr, epsilon=1e-3)
      }else
      {
        S.tr=S.tr
      }
      
      S.te <- t(Z[tr,])%*%Z[tr,]/length(tr)- Sigma_A
      
      for(j in 1:length(tau.path)){
        lambda <- tau.path[j]
        omega=glassoFast(S.tr,lambda)$wi
        sigma=S.te*(1-1/round(n/n.split))
        loss.re[i,j]= loss.re[i,j]+(sum(diag(sigma%*%omega)))-log(det(omega))
      }
    }
    loss.mean <- apply(loss.re, 2, mean)
    loss.sd <- apply(loss.re, 2, sd)
    tau.id=which.min(loss.mean)
    
    
    Sig=t(Z)%*%Z/n- Sigma_A
    if(min(eigen(Sig)$value) < 1e-3)
    {
      pdSig=maxproj.cov(mat=Sig, epsilon=1e-3)  
    }else
    {
      pdSig=Sig
    }
    
    g=glassoFast(pdSig,tau.path[tau.id])  
  }
  else{
    for(i in 1:n.split){
      tr <- sample(1:n, round(n/n.split))
      S.tr <- t(Z[-tr,])%*%Z[-tr,]/(n-length(tr))
      S.te <- t(Z[tr,])%*%Z[tr,]/length(tr)
      
      for(j in 1:length(tau.path)){
        lambda <- tau.path[j]
        omega=glassoFast(S.tr,lambda)$wi
        sigma=S.te*(1-1/round(n/n.split))
        loss.re[i,j]= loss.re[i,j]+(sum(diag(sigma%*%omega)))-log(det(omega))
      }
    }
    loss.mean <- apply(loss.re, 2, mean)
    loss.sd <- apply(loss.re, 2, sd)
    tau.id=which.min(loss.mean)
    Sig=t(Z)%*%Z/n
    
    g=glassoFast(Sig,tau.path[tau.id])  
  }
  return(g)
}



bhq <- function(u, alpha = 0.05) {
  n = length(u)
  r = rank(u, ties.method = "max")
  bh = max(c(r[u <= (r/n) * alpha], 0), na.rm = T)
  su <- sort(u)
  jj <- which(u == 1)
  if (length(jj) != 0) 
    pi0 <- 1
  else pi0 <- min((-1/n) * sum(log(1 - u)), 1)
  if (bh == 0) {
    FDR_BH <- 0
  }
  else {
    FDR_BH <- round((pi0 * su[bh])/(ecdf(u)(su[bh])), 
                    4)
  }
  ad.p = numeric(n)
  ad.p[n] <- sort(u)[n]
  for (i in (n - 1):1) {
    ad.p[i] <- min(sort(u)[i] * (n/i), ad.p[i + 1])
  }
  names(ad.p) = names(sort(u))
  return(c(list(Rejections = bh, FDR = min(FDR_BH, 1), 
                Adjusted.pvalues = sort(ad.p))))
}
