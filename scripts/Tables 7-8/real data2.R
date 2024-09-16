###part 1: data prepocessing
library(dplyr)
rawdata <- read.csv("real data2.csv")
dim(rawdata)
df <- rawdata
na_count <- sapply(df, function(x) sum(is.na(x)))
print(na_count)
non_numeric_columns <- sapply(df, function(x) !is.numeric(x))
print(names(df)[non_numeric_columns])
na_count_per_row <- apply(df, 1, function(x) sum(is.na(x)))
print(na_count_per_row)
row_select <- which(na_count_per_row < 50)
df1 <- df[row_select,]
na_count <- sapply(df1, function(x) sum(is.na(x)))
na_count
na_count[which(na_count < 10)]
df2 <- df1[ , -which(na_count > 100)]
dim(df2)
dff <- df2
dim(dff)
na_count <- sapply(dff, function(x) sum(is.na(x)))
na_count
na_count[which(na_count > 0)]
dff2<- dff[, - which(na_count > 0)]
df<- read.csv("real data2.csv")
dim(df)
values_of_interest <- df[df$RID == 295, "MOCA"]
RID_set <-  dff2$RID 
averages <- df %>%
  filter(RID %in% RID_set) %>%
  group_by(RID) %>%
  summarise(Average_MOCA = mean(MOCA, na.rm = TRUE))
print(averages, n = 224)
dff2 <- merge(dff2, averages, by = "RID", all.x = TRUE)

counts_and_positions <- lapply(RID_set, function(rid) {
  positions <- which(df$RID == rid)
  list(count = length(positions), positions = positions)
})
positions_to_extract <- unlist(lapply(counts_and_positions, function(item) item$positions))
df_selected_rows <- df[positions_to_extract, ]
columns_in_df1 <- names(data[,-1])
df2 <- df_selected_rows[, columns_in_df1]
head(df2)
df_rep <- df2
na_rows_count <- sum(apply(df_rep, 1, function(row) any(is.na(row))))
df_cleaned <- na.omit(df_rep)
var <- scale(df_cleaned[,-1])
df_cleaned[,-1] <- var
df_cleaned <- as.data.frame(df_cleaned)
df_filtered <- df_cleaned %>%
  group_by(RID) %>%
  summarise(count = n(), .groups = 'drop') %>%
  filter(count > 5) %>%
  select(RID) %>%
  left_join(df_cleaned, by = "RID")
rids <- unique(df_filtered$RID)
length(rids)
df_cleaned <- df_filtered
predictors <- setdiff(names(df_cleaned), c("AGE",  "PTGENDER", "PTEDUCAT"))  
df_cleaned <- df_cleaned[ , predictors ]
rids <- unique(df_cleaned$RID)
length(rids)
p <- ncol(df_cleaned) - 1  
cov_matrix <- matrix(0, nrow = p, ncol = p)
total_degrees_of_freedom <- 0

for (rid in rids) {

  subset <- df_cleaned[df_cleaned$RID == rid, -1]
  means <- colMeans(subset)
  for (i in 1:nrow(subset)) {
    deviations <- as.matrix(subset[i, ] - means)
    cov_matrix <- cov_matrix + t(deviations) %*% deviations
  }
  total_degrees_of_freedom <- total_degrees_of_freedom + nrow(subset) - 1
}

sigma_a_hat <- cov_matrix / total_degrees_of_freedom
p <- nrow(sigma_a_hat) 
expanded_matrix <- matrix(0, nrow = p+3, ncol = p+3)
expanded_matrix[(4):(p+3), (4):(p+3)] <- sigma_a_hat
Sigma_A0 <- expanded_matrix
data <- dff2
data <- data %>%
  filter(RID %in% rids)
y0 <- data[,1]
Z0 <- data[,-c(1,2)]
Z0= scale(Z0)
y0 = scale(y0)


###part 2: prediction
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
library(hdi)
library(BDcocolasso)

Z0 <- as.matrix(Z0)
y0 <- y0[[1]]
y0 <- as.vector(y0)
Sigma_A0 <- as.matrix(Sigma_A0)

n <- dim(Z0)[1]
p <- dim(Z0)[2]


iter = 100
obs <- list()
for(i in 1:iter){
  squta = 1 : n
  observation_indices = sample(1:n, 60)
  obs[[i]] = observation_indices
}

relevant <- list()
error <- list()
for (ii in 1:iter) {
 
  observation_indices = obs[[ii]]
  squta = 1:n
  index_prediction = squta[-observation_indices]
  
  X_mid <- as.matrix(Z0[observation_indices,])
  y_mid <- as.vector(y0[observation_indices])
  SigmaA_tr <- Sigma_A0
  
  
  relevant_features <- list()
  tau = Sigma_A0
  # asd = CVglasso(X_mid, tau = 0)
  # Omega= asd$Omega

  asd =cv.glasso_clean(X_mid,0)
  Omega= asd$wi
  
  repo_rank <- Rank_nopre(X_mid, y_mid, SigmaA_tr, Omega = Omega  )
  relevant_features[[1]] <- repo_rank$S
  relevant_features[[2]] <- repo_rank$Splus

  
  y_mid = as.matrix(y_mid)
  ER = Sigma_A0
  
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
  repo <- DCocolasso(X_mid, y_mid, Sigma_A0, Omega = Omega   )
  #reoa <- Co_bhq(X_mid, y_mid)
  relevant_features[[3]] <- repo$ldco


  g_ome = cv.glasso_data(X_mid, Sigma_A0)
  Omega = g_ome$wi
  # asd = CVglasso(X_mid, tau = Sigma_A0)
  # Omega= asd$Omega
  repo_dp <- DPknockoff_nopre(X_mid, y_mid, SigmaA_tr, Omega = Omega)


  relevant_features[[4]] <- repo_dp$S
  relevant_features[[5]] <- repo_dp$Splus


  X_inde <- Z0[(index_prediction), ]
  y <-as.matrix(y0[(index_prediction)])
  
  #rank
  error_1 <- linear_prediction_loss(relevant_features[1:2], X_mid, y_mid, X_inde, y)
  #error_1
  #dcoco and DP
  error_2 <- linear_prediction_loss2(relevant_features[-(1:2)], X_mid, y_mid, X_inde, y, Sigma_A0)

  error_t = c(error_1,error_2)
  # error_t
  # relevant_features
  
  error[[ii]] = error_t
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
colMeans(err_final/19)
apply(err_final/19, 2, std)
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
      
    }
    
    SSE =  sum((c(   cbind(1, x[,relevant_features[[i]]] ) %*% coef ) -  y )^2)
    
    
    
    result[[i]] <- SSE
    
  }
  return(result)
}

linear_prediction_loss2 <- function(relevant_features, X_mid, y_mid, x, y, tau) {
  
  result <- list()
  #prediction_error <- c()
  for (i in 1:length(relevant_features)) {
    print(i)
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
      
    }
    
    SSE =  sum((c(   cbind(1, x[,relevant_features[[i]]] ) %*% coef ) -  y )^2)
    
    result[[i]] <- SSE
    
  }
  return(result)
}



Co_bhq <- function(X, y){
  pval = lasso.proj(X,y)$pval
  bh_adjusted <- p.adjust(pval, method = "BH")
  selected_variables <- which(bh_adjusted < 0.2)
  selected_variables
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
  CV_la=cv.glmnet(Znew,  y,  alpha = 0.9,  nfolds = 5)
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
  
  
  reco1=coco(Znew,y,n,2*p,tau=ER, noise="additive", block=F, K=5, center.Z = F, scale.Z = F, center.y = F, scale.y = F, penalty="lasso")
  beta_co1=reco1$beta.opt
  
  asd = cv.glasso_data(Znew, ER)
  Theta = asd$wi
  
  sigw <- ER
  part2b = n^(-1)*t(Znew)%*%y - (n^(-1)*t(Znew)%*%Znew - sigw)%*%beta_co1
  bwh1 = beta_co1 + Theta%*%part2b
  beta <- bwh1
  
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
  n.split=5
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
        omega=glasso(S.tr,lambda)$wi
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
    
    g=glasso(pdSig,tau.path[tau.id])  
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
