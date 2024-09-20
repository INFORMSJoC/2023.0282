##compare the M-glasso and glasso

data <- readRDS("simu1_data_p400.rds") 
#"simu1_data_p500.rds" "simu1_data_p600.rds"
data_B <- data$beta
data_Z <- data$Z
data_Y <- data$Y
p = 400 #500 600


##simulation example 1
ini_sigma=gen_sigma(p)
Omega0 = inv(ini_sigma)

iter = 100

result_mg = matrix(0, iter, 5)
result_g = matrix(0, iter, 5)
colnames(result_mg) <- c("L2", "LF", "Lmax", "TP", "TN")
colnames(result_g) <- colnames(result_mg)

for(ii in 1:iter){
  Z =  data_Z[[ii]]
  #or use cv.glasso(Z)
  asd1 = CVglasso(Z, tau = 0.25)
  Omega1= asd1$Omega 
  
  #or use cv.glasso_clean(Z)
  asd2 = CVglasso(Z, tau = 0)
  Omega2= asd2$Omega 
  
  re1= glasso_error(Omega1,Omega0)
  re2= glasso_error(Omega2,Omega0)
  

  result_mg[ii,]=re1
  result_g[ii,]=re2
  
  print(result_mg[ii,])
  print(result_g[ii,])
  
}

apply(result_mg, 2, mean)
apply(result_g, 2, mean)
apply(result_mg, 2, sd)
apply(result_g, 2, sd)

