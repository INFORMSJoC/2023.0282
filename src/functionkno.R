bhq <- function(pval, St, q){
  adjusted_p_values <- p.adjust(pval, method = "BH")
  S_de <- which(adjusted_p_values <= q)
  if(length(S_de>0)){
    f2=length(setdiff(S_de,St))/length(S_de)
    t2=length(intersect(S_de,St))/length(St)
  }else{
    f2=0
    t2=0}
  
  return(list(f2 = f2, t2 = t2))
  
}


cv.glasso_clean=function(Z,noise){

  n=nrow(Z)
  p=ncol(Z)
  ntau=20
  n.split=5
  split.ratio=1-1/log(n)
  n.train = ceiling(n*split.ratio)
  n.test = n-n.train
  tau.min = 0.5
  tau.max = 1
  tau.path = seq(tau.min,tau.max,length.out=ntau)*sqrt(log(p)/n)
  loss.re <- matrix(0, n.split, length(tau.path))
  
  K=5
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  folds = sample(cut(seq(1,n),breaks=K,labels=FALSE))
  
    for(i in 1:n.split){
      tr<- which(folds==i, arr.ind= TRUE)
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
  
  return(g)
}


cv.glasso=function(Z,noise){
  
  n=nrow(Z)
  p=ncol(Z)
  ntau=20
  n.split=5
  split.ratio=1-1/log(n)
  n.train = ceiling(n*split.ratio)
  n.test = n-n.train
  tau.min = 0.1
  tau.max = 1
  tau.path = seq(tau.min,tau.max,length.out=ntau)*sqrt(log(p)/n)
  loss.re <- matrix(0, n.split, length(tau.path))
  for(i in 1:n.split){
    tr <- sample(1:n, round(n/n.split))
    #S.tr <- t(Z[-tr,])%*%Z[-tr,]/(n-length(tr))-tau^2*diag(p)
    S.tr <- t(Z[-tr,])%*%Z[-tr,]/(n-length(tr))- ER
    # 
    # if(min(eigen(S.tr)$value) < 1e-3)
    # {
    #   S.tr=maxproj.cov(mat=S.tr, epsilon=1e-3)
    # }else
    # {
    #   S.tr=S.tr
    # }
    
    S.tr=ADMM_proj(S.tr, mu=10, etol=1e-4)$mat
    
    S.te <- t(Z[tr,])%*%Z[tr,]/length(tr)-ER
    # if(min(eigen(S.te)$value) < 1e-3)
    # {
    #   S.te=maxproj.cov(mat=S.te, epsilon=1e-3)
    # }else
    # {
    #   S.te=S.te
    # }
    S.te=ADMM_proj(S.te, mu=10, etol=1e-4)$mat
    
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
  
  
  
  Sig=t(Z)%*%Z/n-ER
  # if(min(eigen(Sig)$value) < 1e-3)
  # {
  #   pdSig=maxproj.cov(mat=Sig, epsilon=1e-3)  #maxproj.cov 
  # }else
  # {
  #   pdSig=Sig
  # }
  pdSig=ADMM_proj(Sig, mu=10, etol=1e-4)$mat
  
  g=glassoFast(pdSig,tau.path[tau.id])
  #g=glasso(pdSig,tau.path[tau.id])  
  return(g)
}



gknockoffX <- function(X, Omega){
  n = nrow(X)
  p = ncol(X)
  # calculate minimum eigenvalue of Sigma = inv(Omega), i.e. 1/max eigenvalue of Omega
  # r = eigen(Sigma)
  # min.eig = r$values[p]
  # Omega=inv(Sigma)
  r=eigen(Omega)
  max.eig = r$values[1]
  max.eig=ifelse(is.numeric(max.eig),max.eig,as.numeric(max.eig))
  # diagonal matrix for construction of knockoff variables
  s = (1/max.eig)*diag(p)
  # min.eig=1/(r$values[1])
  # # diagonal matrix for construction of knockoff variables
  # s = min.eig*diag(p)
  obj = sqrtm(2*s - s%*%Omega%*%s)
  B = obj$B
  A = diag(p) - s%*%Omega
  # now construct knockoff variables conditional on X
  X_ko = X%*%A + matrix(rnorm(n*p),n,p)%*%B
  Xnew = cbind(X, X_ko)
  
  return(Xnew)
}


gen_sigma<-function(p){
  ini_sigma=matrix(0,p,p)
  rho <- 0.5
  for(j in 1:p){
    for(h in 1:p){
      if(j>=h){
        ini_sigma[j,h]<-rho^(abs(j-h))
      }else{
        ini_sigma[j,h]<-0
      }
    }
  }
  ini_diag=diag(ini_sigma)
  ini_sigma=(ini_sigma+t(ini_sigma))
  diag(ini_sigma)=ini_diag
  return(ini_sigma)
}


gen_sigma2<-function(p){
  ini_sigma=matrix(0,p,p)
  a=rep(1,p)
  b=rep(1,p)
  c=a%*%t(b)
  ini_sigma=c+ diag(p)
  return(ini_sigma)
}




glasso_error=function(Omega,Omega_0){
  me1=Omega
  me0=Omega_0
  err=me1 - me0
  err_2=norm(err,"2")
  err_F=norm(err,"F")
  err_max=max(abs(err))
  
  true_nonzero = which(abs(me0) >  1e-3 )
  true_zero = which( abs(me0) <=  1e-3)
  
  es_nonzero = which(abs(me1) >  1e-3 )
  es_zero = which( abs(me1) <=  1e-3)
  
  TP = length(intersect(es_nonzero,true_nonzero)) / length(true_nonzero)
  TN = length(intersect(es_zero,true_zero)) / length(true_zero)
  
  
  # St= true_nonzero
  # S= es_nonzero
  # f1=length(setdiff(S,St))/length(S)
  # t1=length(intersect(S,St))/length(St)
  
  re=c(err_2, err_F, err_max, TP, TN)
  
  return(re)
  #return(list(e2=err_2, ef=err_F, em=err_max, tp=TP, tn=TN,   f1 =f1, t1 = t1))
  
  
}

fdRnew<-function(beta, St, q){
  St=St
  beta=beta
  W = c(rep(0, p))
  for (j in 1:p) {
    W[j] = abs(beta[j]) - abs(beta[j+p])
  }
  
  t_ko=knockoff.threshold(W, fdr  =q, offset = 0)
  S=which(W>=t_ko)
  
  t_ko2=knockoff.threshold(W, fdr  =q, offset = 1)
  S_plus=which(W>=t_ko2)
  
  # set of discovered variables
  
  if(length(S>0)){
    f1=length(setdiff(S,St))/length(S)
    t1=length(intersect(S,St))/length(St)
  }else{
    f1=0
    t1=0
  }
  
  if(length(S_plus>0)){
    f2=length(setdiff(S_plus,St))/length(S_plus)
    t2=length(intersect(S_plus,St))/length(St)
  }else{
    f2=0
    t2=0
  }
  

  
  return(list(f1 =f1, t1 = t1, f2 = f2, t2 = t2))
}


maxproj.cov<-function(mat, epsilon=1e-4, mu=10, nitr.max=1e3, etol=1e-4){
  
  p<-nrow(mat)
  
  # Initialization
  R<-diag(mat)
  S<-matrix(0,p,p)
  L<-matrix(0,p,p)
  
  itr<-0
  while (itr<nitr.max) {
    Rp<-R
    Sp<-S
    
    # Subproblem I: R step
    W<-mat+S+mu*L
    W.eigdec<-eigen(W, symmetric=TRUE)	
    W.V<-W.eigdec$vectors
    W.D<-W.eigdec$values
    R<-W.V%*%diag(pmax(W.D,epsilon))%*%t(W.V)
    
    # Subproblem II: S step
    M<-R-mat-mu*L	
    S[lower.tri(S, diag = TRUE)]<-M[lower.tri(M, diag = TRUE)]-l1proj(v=M[lower.tri(M, diag = TRUE)],b=mu/2)	
    for (i in 2:p){
      for (j in 1:(i-1)){
        S[j,i]<-S[i,j]
      }
    }
    
    # L step: update the Lagrange parameter
    L<-L-(R-S-mat)/mu
    
    # Stopping Rule                        
    #cat("check the stopping criterion:",max(abs(R-S-mat)),"\n")
    if ((max(abs(R-Rp))<etol) && (max(abs(S-Sp))<etol) && (max(abs(R-S-mat))<etol)){
      itr<-nitr.max
    } else {
      itr<-itr+1
    }
    
    if (itr%%20==0) {
      mu<-mu/2
    }
  }
  
  return(R)
  
}

# Efficient projection onto L1 ball of specified radius (i.e. b), used by the admm algo
# Ref. Duchi et al. (2008). Efficient Projections onto the L1-Ball for Learning in High Dimensions, ICML
l1proj<-function(v, b){
  
  stopifnot(b>0)
  
  u <- sort(abs(v),decreasing=TRUE)
  sv <- cumsum(u)
  rho <- max(which(u>(sv-b)/1:length(u)))
  theta <- max(0, (sv[rho]-b)/rho)
  w <-sign(v) * pmax(abs(v)-theta,0)
  
  return(w)
}

#prediction-based selection#


mspbs<-function(X, y, betav){
  
  n=nrow(betav)
  PEvec=c(rep(0,n))
  for (i in 1:n)
  {
    PEvec[i]= t(y-X%*%t(t(betav[i,])))%*%(y-X%*%t(t(betav[i,])))
  }
  id=which.min(PEvec)
  
  return(id)
  
}

